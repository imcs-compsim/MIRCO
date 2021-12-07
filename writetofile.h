#include <string>                        //std::to_string, std::stod
#include <vector>                        // Seems obvious
#include <Epetra_SerialSpdDenseSolver.h> // Seems obvious
#include <Epetra_SerialSymDenseMatrix.h> // Seems obvious
void writeForceToFile(Epetra_SerialDenseMatrix &y, std::string pathName);