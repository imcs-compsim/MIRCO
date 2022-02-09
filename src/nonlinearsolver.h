#include <Epetra_SerialSpdDenseSolver.h>
#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>

class NonLinearSolver {
public:
  void NonlinearSolve(Epetra_SerialDenseMatrix &matrix,
                      Epetra_SerialDenseMatrix &b0, std::vector<double> &y0,
                      Epetra_SerialDenseMatrix &w, Epetra_SerialDenseMatrix &y);
  NonLinearSolver() = default;
};