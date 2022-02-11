#include "evaluate.h"
#include <Teuchos_TestForException.hpp>
#include <string>

int main(int argc, char *argv[]) {
  TEUCHOS_TEST_FOR_EXCEPTION(
      argc != 2, std::invalid_argument,
      "The code expects (only) an input file as argument");
  // reading the input file name from the command line
  std::string inputFileName = argv[1];

  double force = 0.0;

  Evaluate(inputFileName, force);
}