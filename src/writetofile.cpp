#include <Epetra_SerialSymDenseMatrix.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
#include "writetofile.h"

void writeForceToFile(Epetra_SerialDenseMatrix &y, string pathName)
{
  // The aim of this function is to write the output force on all the contact nodes to a file.
  std::size_t botDirPos = pathName.find_last_of("/") + 1;
  // get directory
  std::string dir = pathName.substr(0, botDirPos);
  // get file name
  std::string dat_file = pathName.substr(botDirPos, pathName.length());

  std::string file = dir + "result_force_" + dat_file;

  ofstream outfile;
  outfile.open(file, std::ofstream::trunc);

  for (int i = 0; i < y.M(); i++)
  {
    if (y(i, 0) != 0)
    {
      outfile << y(i, 0) << ";";
    }
  }
}
