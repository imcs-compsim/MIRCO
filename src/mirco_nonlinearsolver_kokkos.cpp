#include "mirco_nonlinearsolver_kokkos.h"

#include <KokkosLapack_gesv.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <vector>

ViewVector_h MIRCO::NonLinearSolver::Solve(
    const ViewMatrix_h& matrix, const ViewVector_h& b0, const ViewVector_h& y0, ViewMatrix_h& w)
{
  const std::string thisFctName =
      "MIRCO::NonLinearSolver::Solve";  // we should do this anywhere we want to use timers or
                                        // kokkos parallel annotation/naming?


  /*
  input H - Influence coefficient matrix (Discrete version of Green Function); symmetric positive
  definite p* - optimal contact force vector (p_{i,j}); p - initial guess? u* - optimal normal
  displacement vector (p_{i,j}); u - initial guess? input \overbar{u} - compenetration vector
  \overbar{u}_{i,j}, i.e. \Delta - (\xi_{max} - \xi_{i,j}) w - := vector w_{i,j}, := u - \overbar{u}

  \epsilon > 0 - tolerance
  K_max - max iterations






  */



  // TODO: make these into input params? at least the important ones
  // TODO: also, make some comments here to describe what exactly they mean
  double nnlstol = 1.0000e-08;
  double maxiter = 10000;
  double eps = 2.2204e-16;
  double alphai = 0;
  double alpha = 100000000;
  int iter = 0;
  bool init = false;
  const int n0 = b0.size();



  // Initialize active set
  std::vector<int> positions;
  Kokkos::View<double*, Device_Host_t> positions for (int i = 0; i < y0.extent(0); i++)
  {
    if (y0(i) >= nnlstol)  // # this should be >= -epsilon, no?
    {
      positions.push_back(i);
      counter++;
    }
  }



  w.reshape(n0, 1);  // # wtf
  w.putScalar(0.0);  // # why

  if (counter == 0)
  {
#if (!kokkosElseOpenMP)
    // add StandardTimer to compare
    StandardTimer timerO("OPENMP __firstParallel");
    timerO.start();
#pragma omp parallel for schedule(static, 16)  // Always same workload -> static
    for (int x = 0; x < n0; x++)
    {
      w(x, 0) = -b0[x];  // # if this is really necessary, use atomic
    }
    // end StandardTimer
    timerO.stop();
#else
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
#endif
    init = false;
  }
  else
  {
#pragma omp parallel for schedule(static, 16)  // Always same workload -> static
    for (int i = 0; i < counter; i++)
    {
      P[i] = positions[i];  // # what
    }

    init = true;
  }

  s0.shape(n0, 1);  // # what


  // New searching algorithm
  std::vector<double> values, newValues;
  std::vector<int> poss, newPositions;
  double minValue = w(0, 0);
  int minPosition = 0;

  bool aux1 = true;
  while (aux1)
  {
    // @{
    for (int i = 0; i < w.numRows(); i++)  // # wtf
    {
      values.push_back(w(i, 0));
      poss.push_back(i);
    }

    // Get all values bigger than initial one
    while (values.size() > 1)
    {
      for (unsigned long int i = 1; i < values.size(); i++)
      {
        if (values[i] < values[0])
        {
          newValues.push_back(values[i]);
          newPositions.push_back(poss[i]);
        }
      }

      if (newValues.size() == 0)
      {
        newValues.push_back(values[0]);
        newPositions.push_back(poss[0]);
      }

      values = newValues;
      poss = newPositions;
      newValues.clear();
      newPositions.clear();
    }
    minValue = values[0];
    minPosition = poss[0];
    values.clear();
    poss.clear();
    // }

    if (((counter == n0) || (minValue > -nnlstol) || (iter >= maxiter)) && (init == false))
    {
      aux1 = false;
    }
    else
    {
      if (init == false)
      {
        P[counter] = minPosition;
        counter += 1;
      }
      else
      {
        init = false;
      }
    }

    int j = 0;

    aux2 = true;
    while (aux2)
    {
      iter++;



      ViewVector_d P("P", n0);

      ViewVector_d s0("s0", n0);
      ViewVector_d vector_bx("bx", counter);
      ViewMatrix_d solverMatrix("A", counter);



      vector_b.size(counter);
      solverMatrix.shape(counter);

#pragma omp parallel for schedule(static, 16)  // Always same workload -> Static!
      for (int x = 0; x < counter; x++)
      {
        vector_b(x) = b0[P[x]];

        for (int z = 0; z < counter; z++)
        {
          if (x >= z)
            solverMatrix(x, z) = matrix(P[x], P[z]);
          else
            solverMatrix(z, x) = matrix(P[x], P[z]);
        }
      }
      // Solving solverMatrix*vector_x=vector_b
      KokkosLapack::gesv()

#pragma omp parallel for schedule(static, 16)  // Always same workload -> Static!
          for (int x = 0; x < counter; x++)
      {
        s0(P[x], 0) = vector_x(x);
      }

      bool allBigger = true;
#pragma omp parallel for schedule(static, 16)  // Always same workload -> Static!
      for (int x = 0; x < counter; x++)
      {
        if (s0(P[x], 0) < nnlstol)
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
          y(P[x]) = s0(P[x], 0);
        }
        w.scale(0.0);
#pragma omp parallel for schedule(dynamic, 16)
        for (int a = 0; a < matrix.numRows(); a++)
        {
          w(a, 0) = 0;
          for (int b = 0; b < counter; b++)
          {
            w(a, 0) += (matrix(a, P[b]) * y(P[b]));
          }
          w(a, 0) -= b0[a];
        }
      }
      else
      {
        j = 0;
        alpha = 1.0e8;

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
              alphai = y(P[i]) / (eps + y(P[i]) - s0(P[i], 0));
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
        for (int a = 0; a < counter; a++) y(P[a]) = y(P[a]) + alpha * (s0(P[a], 0) - y(P[a]));
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
  return y;
}
