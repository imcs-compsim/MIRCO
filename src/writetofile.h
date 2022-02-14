#ifndef SRC_WRITETOFILE_H_
#define SRC_WRITETOFILE_H_

#include <string>
#include <Epetra_SerialSymDenseMatrix.h>
void writeForceToFile(Epetra_SerialDenseMatrix &y, std::string pathName);

#endif //SRC_WRITETOFILE_H_