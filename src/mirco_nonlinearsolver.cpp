#include "mirco_nonlinearsolver.h"
#include <Epetra_SerialSpdDenseSolver.h>
#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>
#include "mirco_linearsolver.h"

void MIRCO::NonLinearSolver::NonlinearSolve(Epetra_SerialDenseMatrix& matrix,
    Epetra_SerialDenseMatrix& b0, Epetra_SerialDenseMatrix& y0, Epetra_SerialDenseMatrix& w,
    Epetra_SerialDenseMatrix& y)
{
  double nnlstol = 1.0000e-08;
  double maxiter = 10000;
  double eps = 2.2204e-16;
  double alphai = 0;
  double alpha = 100000000;
  int iter = 0;
  bool init = false;
  int n0 = b0.M();
  y.Shape(n0, 1);
  Epetra_SerialDenseMatrix s0;
  std::vector<int> P(n0);
  Epetra_SerialDenseMatrix vector_x, vector_b;
  Epetra_SerialSymDenseMatrix solverMatrix;

  // Initialize active set
  std::vector<int> positions;
  for (int i = 0; i < y0.M(); i++)
  {
    if (y0(i, 0) >= nnlstol)
    {
      positions.push_back(i);
    }
  }

  int counter = 0;
  counter = positions.size();
  w.Reshape(b0.M(), b0.N());

  if (counter == 0)
  {
#pragma omp parallel for schedule(static, 16)  // Always same workload -> static
    for (int x = 0; x < b0.M(); x++)
    {
      w(x, 0) = -b0(x, 0);
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

  s0.Shape(n0, 1);
  bool aux1 = true, aux2 = true;

  // New searching algorithm
  std::vector<double> values, newValues;
  std::vector<int> poss, newPositions;
  double minValue = w(0, 0);
  int minPosition = 0;

  while (aux1 == true)
  {
    // @{
    for (int i = 0; i < w.M(); i++)
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
      vector_x.Shape(counter, 1);
      vector_b.Shape(counter, 1);
      solverMatrix.Shape(counter, counter);

#pragma omp parallel for schedule(static, 16)  // Always same workload -> Static!
      for (int x = 0; x < counter; x++)
      {
        vector_b(x, 0) = b0(P[x], 0);

        for (int z = 0; z < counter; z++)
        {
          if (x >= z)
            solverMatrix(x, z) = matrix(P[x], P[z]);
          else
            solverMatrix(z, x) = matrix(P[x], P[z]);
        }
      }
      // Solving solverMatrix*vector_x=vector_b
      MIRCO::LinearSolver solution;
      solution.Solve(solverMatrix, vector_x, vector_b);

#pragma omp parallel for schedule(static, 16)  // Always same workload -> Static!
      for (int x = 0; x < counter; x++)
      {
        s0(P[x], 0) = vector_x(x, 0);
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
          y(P[x], 0) = s0(P[x], 0);
        }
        w.Scale(0.0);
#pragma omp parallel for schedule(dynamic, 16)
        for (int a = 0; a < matrix.M(); a++)
        {
          w(a, 0) = 0;
          for (int b = 0; b < counter; b++)
          {
            w(a, 0) += (matrix(a, P[b]) * y(P[b], 0));
          }
          w(a, 0) -= b0(a, 0);
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
              alphai = y(P[i], 0) / (eps + y(P[i], 0) - s0(P[i], 0));
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
          y(P[a], 0) = y(P[a], 0) + alpha * (s0(P[a], 0) - y(P[a], 0));
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
}
