#include <Epetra_SerialSpdDenseSolver.h>	// Seems obvious
#include <Epetra_SerialSymDenseMatrix.h>	// Seems obvious

#include "linearsolver.h"

void LinearSolver::LinearSolve(Epetra_SerialSymDenseMatrix& matrix,
                 Epetra_SerialDenseMatrix& vector_x,
                 Epetra_SerialDenseMatrix& vector_b) { 
	Epetra_SerialSpdDenseSolver solver; 
	int err = solver.SetMatrix(matrix);
	if (err != 0) {
		std::cout << "Error setting matrix for linear solver (1)";
	}

	err = solver.SetVectors(vector_x, vector_b);
	if (err != 0) {
		std::cout << "Error setting vectors for linear solver (2)";
	}
  
	err = solver.Solve();
	if (err != 0) {
		std::cout << "Error setting up solver (3)";
	}
}