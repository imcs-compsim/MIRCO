#include "mirco_nonlinearsolver_kokkos.h"

#include <KokkosLapack_gesv.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <vector>

ViewVector_h MIRCO::NonLinearSolver::Solve(
    const ViewMatrix_h& matrix, const ViewVector_d& b0_d, const ViewVector_h& p0, ViewVector_h& w_out)
{
  const std::string thisFctName =
      "MIRCO::NonLinearSolver::Solve";  // we should do this anywhere we want to use timers or
                                        // kokkos parallel annotation/naming?



  //- input H - Influence coefficient matrix (Discrete version of Green Function); symmetric
  // positive definite
  //- p* - optimal contact force vector (p_{i,j});
  //- p - initial guess?
  //- u* - optimal normal displacement vector (p_{i,j}); u - initial guess?
  //- input b0_d - compenetration vector with correction, i.e. b0_d(i) = Delta + w_el - (zmax -
  // topology(ri, ci));
  //- w - := vector w_{i,j}, := u - \overbar{u}

  //- \epsilon > 0 - tolerance
  //- K_max - max iterations



  // TODO: make these into input params? at least the important ones
  // TODO: also, make some comments here to describe what exactly they mean (or just in @brief when
  // they are params)
  double nnlstol = 1.0e-08;
  
  
  
  /* IDEA: (nvm)
  * we should use a different tolerance for the $w \ge -epsilon$ and $s_I \ge -epsilon$, because w and s_I have different units (length and pressure respectively). This is not shown in Bemporad 2015, however I think it is shown in the textbook "Soving Least Squares Problems" by Lawson and Hanson, as referenced also in the Bemporad 2015 paper.
  (nvm because the tolerance is also "for" floating point representation error, not "for" physical modelling or solution error (hence no relation to units))
  - strictly speaking, the tolerance check is not needed, but it is there in Bemporad 2015, so we use it
  */
  
  
  
  double maxiter = 10000;
  ///double eps = 2.2204e-16;  // # machine precision of double, i.e. smallest number such that 1.0 + eps > 1.0
  constexpr double eps = 1.0 / (1ULL << 52); // 2^{-52}; double machine precision, i.e. smallest number such that 1.0 + eps > 1.0
  constexpr double infty = 0x7FF0000000000000; // positive infinity for double
  double alphai = 0;
  double alpha = infty;
  int iter = 0;
  bool init = false;
  const int n0 = b0_d.extent(0);



  // Initialize active set
  /*old version with activeSetI.extent(0) = counter
  ViewVectorInt_d tempIndices_d("tempIndices_d", n0);  // Overallocate
  ViewScalarInt_d counter_d("counter_d");              // Holds final count

  Kokkos::parallel_scan(
      "build initial active set", n0, KOKKOS_LAMBDA(const int i, int& update, const bool final) {
        if (p0(i) > 0.0)
        {
          if (final) tempIndices_d(update) = i;
          update++;
        }
        if (final && i == n0 - 1) counter_d = update;
      });

  // Shrink to actual size
  int counter = Kokkos::create_mirror_view_and_copy(MemorySpace_Host_t, counter_d)(0);
  ViewVectorInt_d activeSetI("activeSetI", counter);
  Kokkos::deep_copy(activeSetI, Kokkos::subview(tempIndices_d, std::make_pair(0, counter)));*/

  // new version with activeSetI.extent(0) = n0
  ViewScalarInt_d counter_init_d("counter_init_d");
  ViewVectorInt_d activeSetI("activeSetI", n0);
  Kokkos::parallel_scan(
      "build initial active set", n0, KOKKOS_LAMBDA(const int i, int& update, const bool final) {
        if (p0(i) > 0.0)
        {
          if (final) activeSetI(update) = i;
          update++;
        }
        if (final && i == n0 - 1) counter_d_init = update;
      });

  
  int counter = Kokkos::create_mirror_view_and_copy(MemorySpace_Host_t, counter_d_init)(0);
  
  
  
  
  /* IDEAs for overall kokkos for these data structures which kind of change size:
  * we keep full-sized views on both host and device, and then we simply copy the subviews
  * at least for some things like activeSetI
  
  
What I said to chatGPT:
  
wait ok unrealted but back to the activeSetI thing. Now I have activeSetI a viewvector of length n0 on device but it only has counter valid elements, in general unordered, which say which indices are active. Now, I have my other vectors/matrices like p or w and H, which also need to be like just holding the active ones.

I suppose I have two options for how to treat them:
- I could just basically be keeping eg p as just the full p vector of length n0 and then just fill it completely by using p(activeSetI(i)) when traversing from 0 to counter-1. But I do have one operation which probably wont work with this raw as it is-- KokkosLapack::gesv(), which I mean yeah would need the compacted rhs and matrix (H). But, good thing is I basically never have to move elemnts in the vectors until the end where I sort it all
- I could also make a compacted version which basically matches activeSetI exactly. Then, activeSetI is not used to like give p the positions, but I would instead just call p(i) from 0 to counter-1. The gesv() also requires no change. But, then I have to do exactly the same operations which I do on activeSetI but just also on p, w, and H, when I erase or add an element
  
  
  */


  ViewVector_d w_d("w_d", n0);

  if (counter == 0)
  {
    // add StandardTimer to compare
    StandardTimer timerK("KOKKOS __firstParallel");
    timerK.start();
    Kokkos::parallel_for(
        thisFctName + "__firstParallel", Kokkos::RangePolicy<ExecSpace_Default_t>(0, n0),
        KOKKOS_LAMBDA(const int x) { w_d(x, 0) = -b0_d(x); }); //# change to more succinct synax
    // end StandardTimer
    timerK.stop();

    init = false;
  }
  else
  {
    init = true;
  }

  /// P = copy(positions)



  ViewVector_d p_star("p_star", n0);

  bool aux1 = true;
  while (1)  /// while (aux1)
  {


    double min_w_i;
    int argmin_w_i;
    Kokkos::parallel_reduce(
        "min_w", w.extent(0),
        KOKKOS_LAMBDA(const int i, Kokkos::MinLoc<double, int>::value_type& ml) {
          if (w(i) < ml.val)
          {
            ml.val = w(i);
            ml.loc = i;
          }
        },
        Kokkos::MinLoc<double, int>(min_w_i, argmin_w_i));

    if (((min_w_i >= -epsilon || counter == n0) && init) || iter >= maxiter) //# algo3 line(4)
    {
      break;
    }

    
    
    
    if (init) {
      double min_w_i_fromInactiveSet;
      int argmin_w_i_fromInactiveSet; // in algo3 line(5), we need the argmin where i \in {1, ..., n}\I
      Kokkos::parallel_reduce(
          "min_w", w.extent(0),
          KOKKOS_LAMBDA(const int i, Kokkos::MinLoc<double, int>::value_type& ml) {
            if (w(i) < ml.val)
            {
              ml.val = w(i);
              ml.loc = i;
            }
          },
          Kokkos::MinLoc<double, int>(min_w_i_fromInactiveSet, argmin_w_i_fromInactiveSet));
          
      Kokkos::deep_copy(Kokkos::subview(activeSetI, counter), argmin_w_i_fromInactiveSet);
      ++counter;
    }

    while (2)
    {
      ++iter;



      ViewVector_d P("P", n0);

      /// ViewVector_d s0("s0", n0);//# completely unnecessary

      // first the rhs, then the solution of the linear system matrix*x = b; or, really, H_I*s_I =
      // \overbar{u}_I
      ViewVector_d vector_bx("bx", counter);

      ViewMatrix_d solverMatrix("A", counter);



      vector_b.size(counter);
      solverMatrix.shape(counter);

      Kokkos::parallel_for(
          "assemble reduced system", counter, KOKKOS_LAMBDA(const int x) {
            vector_bx(x) = b0_d(activeSetI(x));
            for (int z = 0; z < counter; z++)
            {
              double val = matrix(activeSetI(x), activeSetI(z));
              solverMatrix(x, z) = val;
            }
          });

      // Solving solverMatrix*vector_x=vector_b
      KokkosLapack::gesv(solverMatrix, vector_bx);

      bool allBigger = true;
#pragma omp parallel for schedule(static, 16)  // Always same workload -> Static!
      for (int x = 0; x < counter;
          x++)  // # opposite of $s_i \ge -\epsilon \forall i \in activeSetI$
      {
        if (vector_bx(x) < -nnlstol)
        {
          allBigger = false;
        }
      }

      if (allBigger == true)
      {
#pragma omp parallel for schedule(guided, 16)
        for (int x = 0; x < counter; x++)
        {
          p_star(P[x]) = vector_bx(x); //# Algo3 line(7)
        }
        w.scale(0.0);
#pragma omp parallel for schedule(dynamic, 16)
        for (int a = 0; a < matrix.numRows(); a++)
        {
          w(a, 0) = 0;
          for (int b = 0; b < counter; b++)
          {
            w(a, 0) += (matrix(a, P[b]) * p_star(P[b]));
          }
          w(a, 0) -= b0_d[a];
        }
        
        break;
      }
      else//# Algo3 line(8)-(11)
      {
        // Searching for minimum value with index position
        double min_alpha_i = infty;
        int argmin_alpha_i = -1;
        
        Kokkos::parallel_reduce(
          "min_alpha", counter,
          KOKKOS_LAMBDA(const int i, Kokkos::MinLoc<double, int>::value_type& ml) {
            if(s0(activeSetI(i)) <= 0) {
              const alphai = p_star(activeSetI(i)) / (eps + p_star(activeSetI(i)) - s0(activeSetI(i)));
              if (alphai < ml.val)
              {
                ml.val = alphai;
                ml.loc = i;
              }
            }
          },
          Kokkos::MinLoc<double, int>(min_alpha_i, argmin_alpha_i));
#pragma omp parallel
        {
          int jP = 0;
          double alphaP = alpha;
#pragma omp parallel for schedule(guided, 16)  // Even, guided seems fitting
          for (int i = 0; i < counter; i++)
          {
            if (s0(activeSetI(i)) < nnlstol) //# Algo3 line(8), i.e. argmin is jP, and min the actual min is alphaP
            {
              alphai = p_star(activeSetI(i)) / (eps + p_star(activeSetI(i)) - s0(activeSetI(i)));
              if (alphai < alphaP)
              {
                alphaP = alphai;
                jP = i;
              }
            }
          }
#pragma omp critical
          {
            alpha = alphaP;
            j = jP;
          }
        }

        // TODO: WHAT BELONGS TO THIS LOOP?????????????????????
#pragma omp parallel for schedule(guided, 16)
        for (int a = 0; a < counter; a++)
          p_star(P[a]) = p_star(P[a]) + alpha * (s0(P[a], 0) - p_star(P[a])); //# Algo3 line(9)
        
        if (j > 0) //# why do we check the last j value? It should be erase I_h where h such that p_h = 0, and that also means all such h, not necessarily only one--though I suppose we should use eps here then too
        {
          //# TODO: think about whether there is a chance of counter ever = 0 in here. Probably not, but if so we need a check ofc. Yeah I think not, but make sure
          // jth entry in P leaves active set
          ///s0(P[j], 0) = 0; //# in fact unnecessary. We can leave the old value there
          Kokkos::deep_copy(Kokkos::subview(activeSetI, j), Kokkos::subview(activeSetI, counter-1));
          ///Kokkos::deep_copy(Kokkos::subview(activeSetI, counter-1), 0); //# this is actually strictly speaking also not needed. We can leave the old value there
          --counter;
        }
      }
    }
  }

  // p_star = p;
  // u_star = w + \overbar{u}



  return p_star;
}
