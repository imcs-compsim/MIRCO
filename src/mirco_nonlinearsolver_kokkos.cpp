#include "mirco_nonlinearsolver_kokkos.h"

#include <KokkosLapack_gesv.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <vector>

ViewVector_h MIRCO::NonLinearSolver::Solve(
    const ViewMatrix_d& matrix, const ViewVector_d& b0_d, const ViewVector_h& p0, ViewVector_h& w_out)
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
    // Hp - b0 for p = 0
    StandardTimer timerK("KOKKOS __firstParallel");
    timerK.start();
    Kokkos::parallel_for(
        thisFctName + "__firstParallel", Kokkos::RangePolicy<ExecSpace_Default_t>(0, n0),
        KOKKOS_LAMBDA(const int x) { w_d(x, 0) = -b0_d(x); }); //# change to more succinct synax
    // end StandardTimer
    timerK.stop();

    init = true; //# note: in the original implementation, init was the opposite of how it is defined/used in Bemporad 2015
  }
  else {
    
  }

  /// P = copy(positions)



  ViewVector_d p_d("p_d", n0);

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

    
    
    
    if (init) { // I empty or iter1 > 0 // if init, we extend the active set by something
      double min_w_i_fromInactiveSet;
      int argmin_w_i_fromInactiveSet; // in algo3 line(5), we need the argmin where i \in {1, ..., n}\I
      Kokkos::parallel_reduce(
          "min_w", Kokkos::RangePolicy<ExecSpace_Default_t>(counter, n0),
          KOKKOS_LAMBDA(const int i, Kokkos::MinLoc<double, int>::value_type& ml) {
            auto wi = w(activeSetI(i))
            if (wi < ml.val)
            {
              ml.val = wi;
              ml.loc = i;
            }
          },
          Kokkos::MinLoc<double, int>(min_w_i_fromInactiveSet, argmin_w_i_fromInactiveSet));
          
      Kokkos::deep_copy(Kokkos::subview(activeSetI, counter), argmin_w_i_fromInactiveSet);
      ++counter;
      
      // separate thing but also fits perfectly into this init
    }
    else {
      
      
      init = true;
    }

    while (2)
    {
      ++iter;
          
          
      ViewMatrix_d H_compact_d("H_compact", counter, counter);
      ViewVector_d b0s_compact_d("b0_compact", counter);

      Kokkos::parallel_for("define compact views", counter, KOKKOS_LAMBDA(const int i) {
        int row = activeSetI(i);
        b0s_compact_d(i) = b0_d(row);
        for (int j = 0; j < counter; ++j) {
          int col = activeSetI(j);
          H_compact_d(i, j) = H(row, col);
        }
      });


      // Solving H_compact_d b0s_compact_d [as output] = b0s_compact_d [as input]
      KokkosLapack::gesv(H_compact_d, b0s_compact_d);

      bool allGreater = true;
      Kokkos::parallel_reduce(
        "allGreater", counter,
        KOKKOS_LAMBDA(const int i, bool& local_allGreater) {
          if (b0s_compact_d(i) < -nnlstol)
          {
            local_allGreater = false;
          }
        },
        allGreater);

      if (allGreater)
      {
        Kokkos::parallel_for("p = s_compact", counter, KOKKOS_LAMBDA(const int i) {
            p_d(activeSet(i)) = b0s_compact_d(i);
          });
          
        // Algo3, line (3) after loop: TWO OPTIONS: (TODO: test if one is generally superior or when each is superior)
        /*constexpr bool option1ElseOption2 = false; //#tmp
        
        // first option: uses 2D range policy but needs two separate loops
        if(option1ElseOption2){
          Kokkos::parallel_for("init w = -b0", n0, KOKKOS_LAMBDA(const int i) {
              const int Ii = activeSetI(i);
              w(i) = -b0(i);
            });
            
          Kokkos::parallel_for(
            "compute H*p", 
            Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {n0, counter}),
            KOKKOS_LAMBDA(const int i, const int j) {
              const int col = activeSetI(j);
              Kokkos::atomic_add(&w(i), H_d(i, activeSetI(col)) * p_compact(activeSetI(col)));
            });
          
        }
          
        // second option: one loop with 1D range policy (probably better option)
        else {*/
          Kokkos::parallel_for(
            "w = H p - b0", n0, KOKKOS_LAMBDA(const int i) {
              double sum = 0.0;
              for(int j=0; j<counter; ++j)
                sum += matrix(i, j) * b0s_compact_d(j);
              w(i) = sum - b0(i);
            });
          
        ///}
        
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
              const alphai = p_d(activeSetI(i)) / (eps + p_d(activeSetI(i)) - s0(activeSetI(i)));
              if (alphai < ml.val)
              {
                ml.val = alphai;
                ml.loc = i;
              }
            }
          },
          Kokkos::MinLoc<double, int>(min_alpha_i, argmin_alpha_i));
          
        Kokkos::parallel_for(
          "assemble reduced system", counter, KOKKOS_LAMBDA(const int i) {
            p_d(activeSetI(i)) = p_d(activeSetI(i)) + alpha * (s0(activeSetI(i)) - p_d(activeSetI(i)))
          });
          
        if (j > -1) //# why do we check the last j value? It should be erase I_h where h such that p_h = 0, and that also means all such h, not necessarily only one--though I suppose we should use eps here then too
        {
          //# TODO: think about whether there is a chance of counter ever = 0 in here. Probably not, but if so we need a check ofc. Yeah I think not, but make sure
          // jth entry in P leaves active set
          p_d(P[j], 0) = 0; //# in fact unnecessary. We can leave the old value there; no, the old value may be 0 only within a tolerance
          Kokkos::deep_copy(Kokkos::subview(activeSetI, j), Kokkos::subview(activeSetI, counter-1));
          Kokkos::deep_copy(Kokkos::subview(activeSetI, j), Kokkos::subview(activeSetI, counter-1));
          ///Kokkos::deep_copy(Kokkos::subview(activeSetI, counter-1), 0); //# this is actually strictly speaking also not needed. We can leave the old value there
          --counter;
        }
      }
    }
  }

  // p_d = p;
  // u_star = w + \overbar{u}



  return p_d;
}
