#include <Epetra_SerialSpdDenseSolver.h>	// Seems obvious
#include <Epetra_SerialSymDenseMatrix.h>	// Seems obvious

class LinearSolver
{
public:
    void LinearSolve(Epetra_SerialSymDenseMatrix& matrix,
                 Epetra_SerialDenseMatrix& vector_x,
                 Epetra_SerialDenseMatrix& vector_b);
    LinearSolver()
    {
        
    }
};