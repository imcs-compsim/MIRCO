#include <omp.h> // omp
#include <unistd.h>
#include <cmath>  //pow
#include <cstdio>
#include <fstream>   //ifstream, ofstream
#include <iostream>  //ifstream, ofstream
#include <string>    //std::to_string, std::stod
#include <vector>
#include "../include/Epetra_SerialSpdDenseSolver.h"
#include "../include/Epetra_SerialSymDenseMatrix.h"
#include <chrono> // time stuff

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

void SetUpMatrix_Static(Epetra_SerialDenseMatrix& A, std::vector<double> xv0,
                 std::vector<double> yv0, double delta, double E,
                 int systemsize, double& time, int cachesize) {
  double r;
  double pi = atan(1) * 4;
  double raggio = delta / 2;
  double C = 1 / (E * pi * raggio);

  // MULTITHREADING
  // Reading in time
  auto start = std::chrono::high_resolution_clock::now();
  
#pragma omp parallel for schedule(static, cachesize)
  for (int i = 0; i < systemsize; i++) {
    A(i, i) = 1 * C;
  }

  // Writing time to console
  auto finish = std::chrono::high_resolution_clock::now();
  time = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
  // Hint: Micro is 10^-6
  
  // std::chrono::duration<double> elapsed = finish - start;
  // time = elapsed.count();
  
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

void SetUpMatrix_Dynamic(Epetra_SerialDenseMatrix& A, std::vector<double> xv0,
                 std::vector<double> yv0, double delta, double E,
                 int systemsize, double& time, int cachesize) {
  double r;
  double pi = atan(1) * 4;
  double raggio = delta / 2;
  double C = 1 / (E * pi * raggio);

  // MULTITHREADING
  // Reading in time
  auto start = std::chrono::high_resolution_clock::now();
  
#pragma omp parallel for schedule(dynamic, cachesize)
  for (int i = 0; i < systemsize; i++) {
    A(i, i) = 1 * C;
  }

  // Writing time to console
  auto finish = std::chrono::high_resolution_clock::now();
  time = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
  // Hint: Micro is 10^-6
  
  // std::chrono::duration<double> elapsed = finish - start;
  // time = elapsed.count();
  
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

void SetUpMatrix_Guided(Epetra_SerialDenseMatrix& A, std::vector<double> xv0,
                 std::vector<double> yv0, double delta, double E,
                 int systemsize, double& time, int cachesize) {
  double r;
  double pi = atan(1) * 4;
  double raggio = delta / 2;
  double C = 1 / (E * pi * raggio);

  // MULTITHREADING
  // Reading in time
  auto start = std::chrono::high_resolution_clock::now();
  
#pragma omp parallel for schedule(guided, cachesize)
  for (int i = 0; i < systemsize; i++) {
    A(i, i) = 1 * C;
  }

  // Writing time to console
  auto finish = std::chrono::high_resolution_clock::now();
  time = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
  // Hint: Micro is 10^-6
  
  // std::chrono::duration<double> elapsed = finish - start;
  // time = elapsed.count();
  
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

void calculateTimes_Static(double& elapsedTime1, double& elapsedTime2, int cachesize, string randomPath) {
	// Setup constants
	int csteps, flagwarm;
	double nu1, nu2, G1, G2, E, G, nu, alpha, H, rnd, k_el, delta, nnodi, to1, E1,
	      E2, lato, zref, ampface, errf;
	double Delta = 50;  // TODO only used for debugging
	
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
	SetUpMatrix_Static(A, xv0, yv0, delta, E, n0, elapsedTime1, cachesize);
	
	// MULTITHREADING
	y.Shape(A.M(), 1);
	auto start = std::chrono::high_resolution_clock::now();
	
#pragma omp parallel for schedule(static, cachesize)
	// y = A * b (Matrix * Vector)
	for (int a = 0; a < A.M(); a++) {
		for (int b = 0; b < A.N(); b++) {
			y(a, 0) += (A(b, a) * b0[b]);
		}
	}
	
	// Writing time to console
	auto finish = std::chrono::high_resolution_clock::now();
	elapsedTime2 = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
	// Hint: Micro is 10^-6
	
	// std::chrono::duration<double> elapsed = finish - start;
	// elapsedTime2 = elapsed.count();
	
	// MULTITHREADING END
}

void calculateTimes_Dynamic(double& elapsedTime1, double& elapsedTime2, int cachesize, string randomPath) {
	// Setup constants
	int csteps, flagwarm;
	double nu1, nu2, G1, G2, E, G, nu, alpha, H, rnd, k_el, delta, nnodi, to1, E1,
	      E2, lato, zref, ampface, errf;
	double Delta = 50;  // TODO only used for debugging
	
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
	SetUpMatrix_Dynamic(A, xv0, yv0, delta, E, n0, elapsedTime1, cachesize);
	
	// MULTITHREADING
	y.Shape(A.M(), 1);
	auto start = std::chrono::high_resolution_clock::now();
	
#pragma omp parallel for schedule(dynamic, cachesize)
	// y = A * b (Matrix * Vector)
	for (int a = 0; a < A.M(); a++) {
		for (int b = 0; b < A.N(); b++) {
			y(a, 0) += (A(b, a) * b0[b]);
		}
	}
	
	// Writing time to console
	auto finish = std::chrono::high_resolution_clock::now();
	elapsedTime2 = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
	// Hint: Micro is 10^-6
	
	// std::chrono::duration<double> elapsed = finish - start;
	// elapsedTime2 = elapsed.count();
	
	// MULTITHREADING END
}

void calculateTimes_Guided(double& elapsedTime1, double& elapsedTime2, int cachesize, string randomPath) {
	// Setup constants
	int csteps, flagwarm;
	double nu1, nu2, G1, G2, E, G, nu, alpha, H, rnd, k_el, delta, nnodi, to1, E1,
	      E2, lato, zref, ampface, errf;
	double Delta = 50;  // TODO only used for debugging
	
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
	SetUpMatrix_Guided(A, xv0, yv0, delta, E, n0, elapsedTime1, cachesize);
	
	// MULTITHREADING
	y.Shape(A.M(), 1);
	auto start = std::chrono::high_resolution_clock::now();
	
#pragma omp parallel for schedule(guided, cachesize)
	// y = A * b (Matrix * Vector)
	for (int a = 0; a < A.M(); a++) {
		for (int b = 0; b < A.N(); b++) {
			y(a, 0) += (A(b, a) * b0[b]);
		}
	}
	
	// Writing time to console
	auto finish = std::chrono::high_resolution_clock::now();
	elapsedTime2 = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
	// Hint: Micro is 10^-6
	
	// std::chrono::duration<double> elapsed = finish - start;
	// elapsedTime2 = elapsed.count();
	
	// MULTITHREADING END
}

double generateMinimum(vector<double> values){
	int min = values[0];
	for (int i = 1; i < values.size(); i++){
		if (min > values[i]) { min = values[i]; }
	}
	return min;
}

void writeToFile(string filepath, vector<double> values, string name){
	ofstream outfile; outfile.open(filepath);
	
	outfile << "Writing " + to_string(values.size()) + " data elements from " + 
			name + " to file." << endl;
	
	for (int i = 0; i < values.size(); i++){
		outfile << to_string(values[i]) << endl;
	}
	
	outfile << "End of data." << endl;
	
	outfile.close();
}

// @overload
void writeToFile(string filepath, Epetra_SerialDenseMatrix values, int dim1, int dim2){
	// Write file in MatLab style to increase data flexibility and handling
	ofstream outfile; outfile.open(filepath);
	outfile << std::scientific;
	for (int y = 0; y < dim2; y++){
		for (int x = 0; x < dim1; x++){
			if (x == (dim1 - 1)){
				outfile << to_string(values(x, y)) << endl;
			} else {
				outfile << to_string(values(x, y)) + ",";
			}
		}
	}
	outfile.close();
}

int main(int argc, char* argv[]) {
	// Setup Thread Amount and Cache_Size
	int maxThreads = 12, maxCache = 18;
	double time1 = 0, time2 = 0, min1 = 0, min2 = 0;
	vector<double> times1, times2, mins1, mins2;
	Epetra_SerialDenseMatrix matrix1, matrix2;
	string filePath = "sup8.dat";
	matrix1.Shape(maxCache, maxThreads); matrix2.Shape(maxCache, maxThreads);
	
	std::cout << "Starting Static Runtime" << endl;
	
	for (int cachesize = 1; cachesize < (maxCache + 1); cachesize++){
		for (int threadAmount = 1; threadAmount < (maxThreads + 1); threadAmount++){
			// Generate runtime-data
			omp_set_num_threads(threadAmount);
			for (int i = 0; i < 100; i++){ // Should be sufficient
				calculateTimes_Static(time1, time2, cachesize, filePath);
				times1.push_back(time1);
				times2.push_back(time2);
			}
			min1 = generateMinimum(times1);
			min2 = generateMinimum(times2);
			// std::cout << to_string(min1) << endl;
			matrix1(cachesize, threadAmount) = min1;
			matrix2(cachesize, threadAmount) = min2;
			min1 = 0; times1.clear();
			min2 = 0; times2.clear();
			std::cout << "Thread " + to_string(threadAmount) + " done." << endl;
		}
		std::cout << "Cache " + to_string(cachesize) + " done." << endl;
	}
	
	writeToFile("datatimes1_static_" + filePath, matrix1, maxCache, maxThreads);
	writeToFile("datatimes2_static_" + filePath, matrix2, maxCache, maxThreads);
	
	std::cout << "Starting Dynamic Runtime" << endl;
		
	for (int cachesize = 1; cachesize < (maxCache + 1); cachesize++){
		for (int threadAmount = 1; threadAmount < (maxThreads + 1); threadAmount++){
			// Generate runtime-data
			omp_set_num_threads(threadAmount);
			for (int i = 0; i < 100; i++){ // Should be sufficient
				calculateTimes_Dynamic(time1, time2, cachesize, filePath);
				times1.push_back(time1);
				times2.push_back(time2);
			}
			min1 = generateMinimum(times1);
			min2 = generateMinimum(times2);
			// std::cout << to_string(min1) << endl;
			matrix1(cachesize, threadAmount) = min1;
			matrix2(cachesize, threadAmount) = min2;
			min1 = 0; times1.clear();
			min2 = 0; times2.clear();
			std::cout << "Thread " + to_string(threadAmount) + " done." << endl;
		}
		std::cout << "Cache " + to_string(cachesize) + " done." << endl;
	}
		
	writeToFile("datatimes1_dynamic_" + filePath, matrix1, maxCache, maxThreads);
	writeToFile("datatimes2_dynamic_" + filePath, matrix2, maxCache, maxThreads);
		
	std::cout << "Starting Guided Runtime" << endl;
			
	for (int cachesize = 1; cachesize < (maxCache + 1); cachesize++){
		for (int threadAmount = 1; threadAmount < (maxThreads + 1); threadAmount++){
			// Generate runtime-data
			omp_set_num_threads(threadAmount);
			for (int i = 0; i < 100; i++){ // Should be sufficient
				calculateTimes_Static(time1, time2, cachesize, filePath);
				times1.push_back(time1);
				times2.push_back(time2);
			}
			min1 = generateMinimum(times1);
			min2 = generateMinimum(times2);
			// std::cout << to_string(min1) << endl;
			matrix1(cachesize, threadAmount) = min1;
			matrix2(cachesize, threadAmount) = min2;
			min1 = 0; times1.clear();
			min2 = 0; times2.clear();
			std::cout << "Thread " + to_string(threadAmount) + " done." << endl;
		}
		std::cout << "Cache " + to_string(cachesize) + " done." << endl;
	}
	
	writeToFile("datatimes1_guided_" + filePath, matrix1, maxCache, maxThreads);
	writeToFile("datatimes2_guided_" + filePath, matrix2, maxCache, maxThreads);
	
	std::cout << "All jobs done!" << endl;
}