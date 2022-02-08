#include <fstream>
#include <iostream>
#include <string>
using namespace std;
#include "evaluate.h"

int main(int argc, char *argv[])
{
  string jsonFileName = argv[1]; // reading the json file name from the command line

  double force = 0;

  if (argc==2)
  {
    Evaluate(jsonFileName, force);
  }
  else
  {
    cout << "Invalid input arguments" << endl;
  }
  
}