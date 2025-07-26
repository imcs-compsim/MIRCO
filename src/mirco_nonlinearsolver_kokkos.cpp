#include "mirco_nonlinearsolver_kokkos.h"

#include <KokkosLapack_gesv.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <vector>

ViewVector_h MIRCO::NonLinearSolver::Solve(
    const ViewMatrix_h& matrix, const ViewVector_h& b0, const ViewVector_h& p0, ViewVector_h& w)
{
  const std::string thisFctName =
      "MIRCO::NonLinearSolver::Solve";  // we should do this anywhere we want to use timers or
                                        // kokkos parallel annotation/naming?



  //- input H - Influence coefficient matrix (Discrete version of Green Function); symmetric
  // positive definite
  //- p* - optimal contact force vector (p_{i,j});
  //- p - initial guess?
  //- u* - optimal normal displacement vector (p_{i,j}); u - initial guess?
  //- input b0 - compenetration vector with correction, i.e. b0(i) = Delta + w_el - (zmax -
  // topology(ri, ci));
  //- w - := vector w_{i,j}, := u - \overbar{u}

  //- \epsilon > 0 - tolerance
  //- K_max - max iterations



  // TODO: make these into input params? at least the important ones
  // TODO: also, make some comments here to describe what exactly they mean (or just in @brief when
  // they are params)
  double nnlstol = 1.0000e-08;
  double maxiter = 10000;
  double eps = 2.2204e-16;  // #what is this
  double alphai = 0;
  double alpha = 100000000;
  int iter = 0;
  bool init = false;
  const int n0 = b0.extent(0);



  // Initialize active set
  ViewVectorInt_d tempIndices_d("tempIndices_d", n0);  // Overallocate
  ViewScalarInt_d counter_d("counter_d");              // Holds final count

  Kokkos::parallel_scan(
      "build active set", n0, KOKKOS_LAMBDA(const int i, int& update, const bool final) {
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
  Kokkos::deep_copy(activeSetI, Kokkos::subview(tempIndices_d, std::make_pair(0, counter)));



  w_h = ViewVector_h("w", n0);

  if (counter == 0)
  {
    // # TODO: change this
    auto w_kokkos = toKokkos(w);
    auto b0_kokkos = toKokkos(b0);
    // add StandardTimer to compare
    StandardTimer timerK("KOKKOS __firstParallel");
    timerK.start();
    Kokkos::parallel_for(
        thisFctName + "__firstParallel", Kokkos::RangePolicy<ExecSpace_Default_t>(0, n0),
        KOKKOS_LAMBDA(const int x) { w_kokkos(x, 0) = -b0_kokkos(x); });
    // end StandardTimer
    timerK.stop();
    w = kokkosMatrixToTeuchos(w_kokkos);

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
    w.scale(0.0);
#pragma omp parallel for schedule(dynamic, 16)
    for (int a = 0; a < matrix.numRows(); a++)
    {
      w(a, 0) = 0;
      for (int b = 0; b < counter; b++)
      {
        w(a, 0) += (matrix(a, P[b]) * p_star(P[b]));
      }
      w(a, 0) -= b0[a];
    }



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

    // }

    if (((min_w_i >= -epsilon || counter == n0) && init) || iter >= maxiter)
    {
      break;
    }

    if (init) int j = 0;

    aux2 = true;
    while (aux2)
    {
      iter++;



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
            vector_bx(x) = b0(activeSetI(x));
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
        aux2 = false;
#pragma omp parallel for schedule(guided, 16)
        for (int x = 0; x < counter; x++)
        {
          p_star(P[x]) = vector_bx(x);
        }
      }
      else
      {
        j = 0;
        alpha = 1.0e8;  // #wtf

        // Searching for minimum value with index position
#pragma omp parallel
        {
          int jP = 0;
          double alphaP = alpha;
#pragma omp parallel for schedule(guided, 16)  // Even, guided seems fitting
          for (int i = 0; i < counter; i++)
          {
            if (s0(P[i], 0) < nnlstol)
            {
              alphai = p_star(P[i]) / (eps + p_star(P[i]) - s0(P[i], 0));
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
          p_star(P[a]) = p_star(P[a]) + alpha * (s0(P[a], 0) - p_star(P[a]));
        if (j > 0)
        {
          // jth entry in P leaves active set
          s0(P[j], 0) = 0;
          P.erase(P.begin() + j);
#pragma omp atomic  // Necessary?
          counter -= 1;
        }
      }
    }
  }

  // p_star = p;
  // u_star = w + \overbar{u}



  return p_star;
}
