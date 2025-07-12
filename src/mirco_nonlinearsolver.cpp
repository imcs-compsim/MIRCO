#include "mirco_nonlinearsolver.h"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialSymDenseMatrix.hpp>
#include <vector>

#include "mirco_linearsolver.h"

// tmp
#include "tmpHelpers/Timer.hpp"
#include "tmpHelpers/kokkosIntegration.hpp"

Teuchos::SerialDenseVector<int, double> MIRCO::NonLinearSolver::Solve(
    const Teuchos::SerialDenseMatrix<int, double>& matrix, const std::vector<double>& b0,
    const Teuchos::SerialDenseMatrix<int, double>& y0, Teuchos::SerialDenseMatrix<int, double>& w)
{
  double nnlstol = 1.0000e-08;
  double maxiter = 10000;
  double eps = 2.2204e-16;
  double alphai = 0;
  double alpha = 100000000;
  int iter = 0;
  bool init = false;
  const int n0 = b0.size();
  Teuchos::SerialDenseMatrix<int, double> s0;
  std::vector<int> P(n0);
  Teuchos::SerialDenseVector<int, double> vector_x, vector_b;
  Teuchos::SerialSymDenseMatrix<int, double> solverMatrix;

  // Initialize contact force solution
  Teuchos::SerialDenseVector<int, double> y;
  y.size(n0);
  y.putScalar(0.0);

  // Initialize active set
  std::vector<int> positions;
  for (int i = 0; i < y0.numRows(); i++)
  {
    if (y0(i, 0) >= nnlstol)
    {
      positions.push_back(i);
    }
  }

  int counter = 0;
  counter = positions.size();
  w.reshape(n0, 1);
  w.putScalar(0.0);

  if (counter == 0)
  {
#pragma omp parallel for schedule(static, 16)  // Always same workload -> static
    for (int x = 0; x < n0; x++)
    {
      w(x, 0) = -b0[x];
    }

    init = false;
  }
  else
  {
#pragma omp parallel for schedule(static, 16)  // Always same workload -> static
    for (int i = 0; i < counter; i++)
    {
      P[i] = positions[i];
    }

    init = true;
  }

  s0.shape(n0, 1);
  bool aux1 = true, aux2 = true;

  // New searching algorithm
  std::vector<double> values, newValues;
  std::vector<int> poss, newPositions;
  double minValue = w(0, 0);
  int minPosition = 0;

  while (aux1 == true)
  {
    // @{
    for (int i = 0; i < w.numRows(); i++)
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
    while (aux2 == true)
    {
      iter++;
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
      vector_x = MIRCO::LinearSolver::Solve(solverMatrix, vector_b);

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
