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

/*
void TestingCode() { // Contains former test code
	double time=10000;

	  for(int i=0; i<1000; i++) {
		  const int size = 262144;
		  double sinTable[size];
		  auto start = std::chrono::high_resolution_clock::now();
	  
	#pragma omp parallel for
		  for(int n=0; n<size; ++n) {
			  sinTable[n] = std::sin(2 * M_PI * n / size);
		  }

		  auto finish = std::chrono::high_resolution_clock::now();
		  std::chrono::duration<double> elapsed = finish - start;

		  if (elapsed.count()<time) { time=elapsed.count(); }
	  }

	  std::cout << "Elapsed time: " << time << " s\n";
		
	  	  #pragma omp parallel 
	  { 
		  printf("Hello World from thread %d\n", omp_get_thread_num());
		  if (omp_get_thread_num() == 0) { printf("Number of threads is %d\n", omp_get_num_threads()); 
	  }
	  }
}
*/

void CreateTopology(int systemsize, Epetra_SerialDenseMatrix& topology,
                    string filePath) {
  // Readin for amount of lines -> dimension of matrix
  ifstream reader(filePath);
  string blaLine;
  int dimension = 0;
  while (getline(reader, blaLine)) {
    dimension += 1;
  }
  reader.close();
  topology.Shape(dimension, dimension);
  int lineCounter = 0;
  float elements[264];
  int position = 0;
  ifstream stream(filePath);
  string line;
  while (getline(stream, line)) {
    // Split up Values into Double-Array
    int separatorPosition = 0;
    lineCounter += 1;  // Has to happen here, since baseline value is 0.

    for (int i = 0; i < dimension; i++) {  // prevent duplication of values!
      separatorPosition = line.find_first_of(';');
      string container = line.substr(0, separatorPosition);
      line = line.substr(separatorPosition + 1, line.length());

      double value = stod(container);

      topology(lineCounter - 1, i) = value;

      position += 1;
    }
  }
  stream.close();
}

// Generating necessary constants
void SetParameters(double& E1, double& E2, int& csteps, int& flagwarm,
                   double& lato, double& zref, double& ampface, double& nu1,
                   double& nu2, double& G1, double& G2, double& E, double& G,
                   double& nu, double& alpha, double& H, double& rnd,
                   double& k_el, double& delta, double& nnodi, double& errf,
                   double& to1) {
  E1 = 1;
  E2 = 1;
  nu1 = 0.3;
  nu2 = 0.3;
  G1 = E1 / (2 * (1 + nu1));
  G2 = E2 / (2 * (1 + nu2));
  E = 1 / ((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2) / E2));
  G = 1 / ((2 - nu1) / (4 * G1) + (2 - nu2) / (4 * G2));
  nu = E / (2 * G) - 1;

  vector<double> alpha_con{0.778958541513360, 0.805513388666376,
                           0.826126871395416, 0.841369158110513,
                           0.851733020725652, 0.858342234203154,
                           0.862368243479785, 0.864741597831785};
  int nn = 2;  // Matrix sent has the parameter nn=2!
  alpha = alpha_con[nn - 1];
  csteps = 1;
  ampface = 1;
  flagwarm = 0;
  lato = 1000;  // Lateral side of the surface [micrometers]
  H = 0.1;      // Hurst Exponent (D = 3 - H)
  rnd = 95.0129;
  zref = 50;  // Reference for the Scaling, former value = 25
  k_el = lato * E / alpha;
  delta = lato / (pow(2, nn) + 1);
  nnodi = pow(pow(2, nn + 1), 2);

  errf = 100000000;
  to1 = 0.01;
}

// Generating matrix A
void SetUpMatrix(Epetra_SerialDenseMatrix& A, std::vector<double> xv0,
                 std::vector<double> yv0, double delta, double E,
                 int systemsize, double& time) {
  double r;
  double pi = atan(1) * 4;
  double raggio = delta / 2;
  double C = 1 / (E * pi * raggio);

  // MULTITHREADING
  // Reading in time
  auto start = std::chrono::high_resolution_clock::now();
  
#pragma omp parallel for
  for (int i = 0; i < systemsize; i++) {
    A(i, i) = 1 * C;
  }

  // Writing time to console
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  time = elapsed.count();
  // std::cout << "\n	Calculating time for setting up A. \n";
  // std::cout << "Elapsed time: " << time << " s\n"; 
  // MULTITHREADING END
  		  
  for (int i = 0; i < systemsize; i++) {
    for (int j = 0; j < i; j++) {
      r = sqrt(pow((xv0[j] - xv0[i]), 2) + pow((yv0[j] - yv0[i]), 2));
      A(i, j) = C * asin(raggio / r);
      A(j, i) = C * asin(raggio / r);
    }
  }
}

void calculateTimes(double& elapsedTime1, double& elapsedTime2) {
	// Setup constants
	int csteps, flagwarm;
	double nu1, nu2, G1, G2, E, G, nu, alpha, H, rnd, k_el, delta, nnodi, to1, E1,
	      E2, lato, zref, ampface, errf;
	double Delta = 50;  // TODO only used for debugging
	string randomPath = "sup2.dat";
	
	// std::cout << "Test file for generating data is: " + randomPath + ".\n";

	SetParameters(E1, E2, csteps, flagwarm, lato, zref, ampface, nu1, nu2, G1, G2,
	                E, G, nu, alpha, H, rnd, k_el, delta, nnodi, errf, to1);
	vector<double> x;
	for (double i = delta / 2; i < lato; i = i + delta) {
		x.push_back(i);
	}

	// Setup Topology
	Epetra_SerialDenseMatrix topology, y;
	CreateTopology(topology.N(), topology, randomPath);
	
	double zmax = 0;
	double zmean = 0;
	
	for (int i = 0; i < topology.M(); i++) {
		for (int j = 0; j < topology.N(); j++) {
				zmean += topology(i, j);
				if (topology(i, j) > zmax) zmax = topology(i, j);
		}
	}
	zmean = zmean / (topology.N() * topology.M());
	
	vector<double> force0, area0;
	double w_el = 0, area = 0, force = 0;
	int k = 0;
	int n0;
	std::vector<double> xv0, yv0, b0, x0, xvf, yvf, pf;  // x0: initialized in Warmstart!
	double nf, xvfaux, yvfaux, pfaux;
	Epetra_SerialDenseMatrix A;
	// At this point would have been a while-loop, but not necessary since only one calculation
	// is necessary
	
	// while(...)
	vector<int> col, row;
	double value = zmax - Delta - w_el;
	
	for (int i = 0; i < topology.N(); i++) {
		//      cout << "x= " << x[i] << endl;
		for (int j = 0; j < topology.N(); j++) {
			if (topology(i, j) >= value) {
				row.push_back(i);
				col.push_back(j);
			}
		}
	}
	n0 = col.size();
	
	for (int b = 0; b < n0; b++) {
		xv0.push_back(x[col[b]]);
	}
	for (int b = 0; b < n0; b++) {
		yv0.push_back(x[row[b]]);
	}
	
	for (int b = 0; b < n0; b++) {
		b0.push_back(Delta + w_el - (zmax - topology(row[b], col[b])));
	}
	
	int err = A.Shape(xv0.size(), xv0.size());
	
	// Construction of the Matrix H = A
	SetUpMatrix(A, xv0, yv0, delta, E, n0, elapsedTime1);
	
	// MULTITHREADING
	y.Shape(A.N(), 0);
	auto start = std::chrono::high_resolution_clock::now();
	
#pragma omp parallel for
	// y = A * b (Matrix * Vector)
	for (int a = 0; a < A.M(); a++) {
		for (int b = 0; b < A.N(); b++) {
			y(a, 0) += (A(b, a) * b0[b]);
		}
	}
	
	// Writing time to console
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	elapsedTime2 = elapsed.count();
	// std::cout << "\n	Calculating time for y = A * b0. \n";
	// std::cout << "Elapsed time: " << elapsedTime2 << " s\n"; 
	// MULTITHREADING END
}

double generateMean(vector<double> values){
	double sum = 0; int elements = values.size();
	for (int i = 0; i < elements; i++){
		sum += values[i];
	}
	sum = sum / elements;
	return sum;
}

void sortValues(vector<double>& values){ // Impossible to parallelize
	int elements = values.size();
	double min = values[0]; int position = 0;
	vector<double> sortedValues;
	for (int a = 0; a < elements; a++){
		// Find minimum, then insert into new vector and erase from old vector.
		min = values[0]; position = 0;
		for (int i = 1; i < values.size(); i++){
			if (min > values[i]) { min = values[i]; position = i; }
		}
		sortedValues.push_back(min);
		values.erase(values.begin() + position);
	}
	values = sortedValues;
}

int main(int argc, char* argv[]) {
	// Setup Thread Amount
	int threadAmount = 1;
	omp_set_num_threads(threadAmount);
	double time1 = 0, time2 = 0;
	vector<double> times1, times2;
	// Generate different time-data
	for (int i = 0; i < 100; i++){
		calculateTimes(time1, time2);
		times1.push_back(time1);
		times2.push_back(time2);
	}
	
	std::cout << "Mean time for setting up A is:" + to_string(generateMean(times1)) + "s. \n";
	std::cout << "Mean time for calculating y = A * b0 is: " + to_string(generateMean(times2)) + "s. \n";
	
	sortValues(times1);
	sortValues(times2);
	
	// Insert method to write data to files.
}