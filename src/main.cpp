#include <fstream>
#include <iostream>
#include <string>
using namespace std;
#include "evaluate.h"
#include <Teuchos_TestForException.hpp>

int main(int argc, char *argv[])
{
  TEUCHOS_TEST_FOR_EXCEPTION(argc!=2, std::invalid_argument, "Number of command line arguments does not equal 2.");
  
  string jsonFileName = argv[1]; // reading the json file name from the command line

  double force = 0;

  Evaluate(jsonFileName, force);
  
}