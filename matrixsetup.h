#include <vector>	// Seems obvious
#include <Epetra_SerialSymDenseMatrix.h>	// Seems obvious

class MatrixGeneration
{
public:
    void SetUpMatrix(Epetra_SerialDenseMatrix& A, std::vector<double> xv0,
                 std::vector<double> yv0, double delta, double E,
                 int systemsize, int k);
    MatrixGeneration()
    {

    }
};