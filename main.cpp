#include <fstream>                       //ifstream
#include <iostream>                      //ifstream
#include <string>                        //std::to_string, std::stod
using namespace std;
#include "evaluate.h"

int main(int argc, char *argv[])
{
  string jsonFileName = argv[1]; // reading the json file name from the command line

  double force = 0;
  

  Evaluate(jsonFileName, force);

  cout << "The force is: " << force << endl;
  
}