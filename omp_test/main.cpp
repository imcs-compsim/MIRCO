#include <omp.h>
#include <unistd.h>
#include <cmath>  //pow
#include <cstdio>
#include <fstream>   //ifstream
#include <iostream>  //ifstream
#include <string>    //std::to_string, std::stod
#include <vector>
#include "../include/Epetra_SerialSpdDenseSolver.h"
#include "../include/Epetra_SerialSymDenseMatrix.h"
#include <chrono>

using namespace std;


int main(int argc, char* argv[]) {
  omp_set_num_threads(1);

double time=10000;

for(int i=0; i<1000; i++)
{
  const int size = 262144;
  double sinTable[size];
  
auto start = std::chrono::high_resolution_clock::now();
  
  #pragma omp parallel for
  for(int n=0; n<size; ++n)
    sinTable[n] = std::sin(2 * M_PI * n / size);

auto finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed = finish - start;

if (elapsed.count()<time) time=elapsed.count();
}

std::cout << "Elapsed time: " << time << " s\n";

  #pragma omp parallel
{
  printf("Hello World from thread %d\n", omp_get_thread_num());
  if (omp_get_thread_num() == 0) printf("Number of threads is %d\n", omp_get_num_threads());
}

}
