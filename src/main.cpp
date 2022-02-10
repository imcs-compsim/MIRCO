// #include <fstream>
// #include <iostream>
#include <string>
#include "evaluate.h"

int main(int argc, char *argv[])
{
  // reading the input file name from the command line
  std::string inputFileName = argv[1];

  double force = 0.0;

  Evaluate(inputFileName, force);

}