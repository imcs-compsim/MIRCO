#include "mirco_linearsolver.h"

#include <Teuchos_RCP.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialSpdDenseSolver.hpp>
#include <Teuchos_SerialSymDenseMatrix.hpp>

Teuchos::SerialDenseVector<int, double> MIRCO::LinearSolver::Solve(
    Teuchos::SerialSymDenseMatrix<int, double>& matrix,
    Teuchos::SerialDenseVector<int, double>& vector_b)
{
  Teuchos::SerialSpdDenseSolver<int, double> solver;
  int err = solver.setMatrix(Teuchos::rcpFromRef(matrix));
  if (err != 0)
  {
    std::cout << "Error setting matrix for linear solver (1)";
  }

  Teuchos::SerialDenseVector<int, double> vector_x;
  vector_x.size(vector_b.length());

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

  return vector_x;
}
