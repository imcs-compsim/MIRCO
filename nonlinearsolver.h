#include <vector>	// Seems obvious
#include <Epetra_SerialSpdDenseSolver.h>	// Seems obvious
#include <Epetra_SerialSymDenseMatrix.h>	// Seems obvious

class NonLinearSolver
{
public:
    void NonlinearSolve(Epetra_SerialDenseMatrix& matrix,
                    Epetra_SerialDenseMatrix& b0, std::vector<double>& y0,
                    Epetra_SerialDenseMatrix& w, Epetra_SerialDenseMatrix& y);
    NonLinearSolver()
    {
        
    }
};