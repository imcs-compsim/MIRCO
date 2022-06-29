#include "warmstart.h"
#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>

using namespace std;

void MIRCO::Warmstart(Epetra_SerialDenseMatrix& x0, Epetra_SerialDenseMatrix xv0,
    Epetra_SerialDenseMatrix yv0, Epetra_SerialDenseMatrix& xvf, Epetra_SerialDenseMatrix& yvf,
    Epetra_SerialDenseMatrix& pf)
{
  // This function is used to determine the nodes which were in contact in the last iteration from
  // the current contact set. This helps in making an initial guess of the nodes in contact in the
  // current iteration and speeds up the computation.
  x0.Shape(xv0.N(), 1);

  for (int i = 0; i < xv0.N(); i++)
  {
    for (int j = 0; j < xvf.N(); j++)
    {
      if (xvf(0, j) == xv0(0, i) && yvf(0, j) == yv0(0, i))
      {
        x0(i, 0) = pf(0, j);
      }
    }
  }
}
