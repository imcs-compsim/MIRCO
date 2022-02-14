#ifndef SRC_WRITETOFILE_H_
#define SRC_WRITETOFILE_H_

#include <Epetra_SerialSymDenseMatrix.h>
#include <string>
void writeForceToFile(Epetra_SerialDenseMatrix &y, std::string pathName);

#endif  // SRC_WRITETOFILE_H_
