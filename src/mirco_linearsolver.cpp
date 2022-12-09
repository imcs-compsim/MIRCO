#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_RCP.hpp>

#include "mirco_linearsolver.h"

void MIRCO::LinearSolver::Solve(Teuchos::SerialDenseMatrix<int,double>& matrix,
    Teuchos::SerialDenseMatrix<int,double>& vector_x, Teuchos::SerialDenseMatrix<int,double>& vector_b)
{
  Teuchos::SerialDenseSolver<int,double> solver;
  int err = solver.setMatrix(Teuchos::rcpFromRef(matrix));
  if (err != 0)
  {
    std::cout << "Error setting matrix for linear solver (1)";
  }

  err = solver.setVectors(Teuchos::rcpFromRef(vector_x), Teuchos::rcpFromRef(vector_b));
  if (err != 0)
  {
    std::cout << "Error setting vectors for linear solver (2)";
  }

  solver.factorWithEquilibration(true);
  err = solver.solve();
  if (err != 0)
  {
    std::cout << "Error setting up solver (3)";
  }
}
