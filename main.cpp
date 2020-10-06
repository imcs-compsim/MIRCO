#include <omp.h>	// Seems obvious
#include <unistd.h>		// Linux stuff
#include <cmath>  //pow
#include <cstdio>	// IO stuff
#include <fstream>   //ifstream
#include <iostream>  //ifstream
#include <string>    //std::to_string, std::stod
#include <vector>	// Seems obvious
#include "include/Epetra_SerialSpdDenseSolver.h"	// Seems obvious
#include "include/Epetra_SerialSymDenseMatrix.h"	// Seems obvious
#include <chrono> // time stuff
#include <ctime>
using namespace std;

// Declaration for std::vector<int> reduction in parallel loops.
#pragma omp declare reduction(mergeI:std::vector<int>:omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

// Declaration for std::vector<double> reduction in parallel loops.
#pragma omp declare reduction(mergeD:std::vector<double>:omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

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
  int nn = 8;  // CHANGE THIS WHEN CHANGING FILES
  alpha = alpha_con[nn - 1];
  csteps = 1;
  ampface = 1;
  flagwarm = 1;
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

/*------------------------------------------*/
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
  float elements[514];
  int position = 0, separatorPosition, lineCounter = 0;
  ifstream stream(filePath);
  string line, container;
  double value;
  while (getline(stream, line)) {
    separatorPosition = 0;
    lineCounter += 1;
    for (int i = 0; i < dimension; i++) {    // Parallel not possible, since no synchronization points possible
    	separatorPosition = line.find_first_of(';');
    	container = line.substr(0, separatorPosition);
    	position += 1;
    	line = line.substr(separatorPosition + 1, line.length());
    	value = stod(container);
    	topology(lineCounter - 1, i) = value;
    }
  }
  stream.close();
}

/*------------------------------------------*/
void SetUpMatrix(Epetra_SerialDenseMatrix& A, std::vector<double> xv0,
                 std::vector<double> yv0, double delta, double E,
                 int systemsize, int k) {
	double r, pi = atan(1) * 4, raggio = delta / 2;
	double C = 1 / (E * pi * raggio);

#pragma omp parallel for schedule (static, 16) // Always same workload -> static
	for (int i = 0; i < systemsize; i++) {
		A(i, i) = 1 * C;
	}
	
#pragma omp parallel for schedule (static, 16) private(r) // Always same workload -> static
	// Every iteration needs to have a different r! -> private(r)
	for (int i = 0; i < systemsize; i++) {
		for (int j = 0; j < i; j++) {
			r = sqrt(pow((xv0[j] - xv0[i]), 2) + pow((yv0[j] - yv0[i]), 2));
			A(i, j) = C * asin(raggio / r);
			A(j, i) = C * asin(raggio / r);
		}
	}
}

void writeToFile(string time, string fileName, int threadAmount){
	ofstream outfile;
	outfile.open("study_data_" + fileName, std::ofstream::app);
	outfile << "(" + to_string(threadAmount) + ") Time for " + fileName + " is: " + time + " seconds." << endl;
	outfile.close();
}

Epetra_SerialDenseMatrix Warmstart(Epetra_SerialDenseMatrix xv0, Epetra_SerialDenseMatrix yv0, Epetra_SerialDenseMatrix& xvf,
				Epetra_SerialDenseMatrix& yvf, Epetra_SerialDenseMatrix& pf) {
    Epetra_SerialDenseMatrix x0; x0.Shape(xv0.N(), 1);
    Epetra_SerialDenseMatrix combinedMatrix;
    combinedMatrix.Shape(2, xv0.N());
    // matfin = [xv0, yv0]
#pragma omp parallel for schedule (static, 16) // Always same workload -> Static
    for (int i = 0; i < xv0.N(); i++) {
    	combinedMatrix(0, i) = xv0(0, i);
    }
    
#pragma omp parallel for schedule (static, 16)
    for (int i = 0; i < yv0.N(); i++){
    	combinedMatrix(1, i) = yv0(0, i);
    }

    vector<int> index;
#pragma omp parallel for schedule (dynamic, 16) // Workload can differ vastly -> Dynamic
    for (int i = 0; i < pf.N(); i++) {
        // ind=find(matfin(:,1)==xvf(i) & matfin(:,2)==yvf(i));
        for (int j = 0; j < xvf.N(); j++) {
            if ((combinedMatrix(j, 0) == xvf(i, 0)) && (combinedMatrix(j, 1) == yvf(i, 0))) {
                index.push_back(j);
            }
        }

        // x0(ind,1)=pf(i);
#pragma omp parallel for schedule(static, 16) // Always same workload -> Static
        for (int y = 0; y < index.size(); y++) {
            x0(y, 0) = pf(y, 0);
        }
        index.clear();
    }
    return x0;
}

/*------------------------------------------*/
// Parallel not possible, since multiple accesses on the same object.
// In addition, parallel programming not supported for LinearSolve.

void LinearSolve(Epetra_SerialSymDenseMatrix& matrix,
                 Epetra_SerialDenseMatrix& vector_x,
                 Epetra_SerialDenseMatrix& vector_b) { 
	Epetra_SerialSpdDenseSolver solver; 
	int err = solver.SetMatrix(matrix);
	if (err != 0) {
		std::cout << "Error setting matrix for linear solver (1)";
	}

	err = solver.SetVectors(vector_x, vector_b);
	if (err != 0) {
		std::cout << "Error setting vectors for linear solver (2)";
	}
  
	err = solver.Solve();
	if (err != 0) {
		std::cout << "Error setting up solver (3)";
	}
}

/*------------------------------------------*/
double NonlinearSolve(Epetra_SerialDenseMatrix& matrix, string filename, 
                    Epetra_SerialDenseMatrix& b0, std::vector<double>& y0,
                    Epetra_SerialDenseMatrix& w, Epetra_SerialDenseMatrix& y) {
	// matrix -> A, b0 -> b, y0 -> y0 , y -> y, w-> w; nnstol, iter, maxiter ->
	// unused
	double nnlstol = 1.0000e-08; double maxiter = 10000;
	double eps = 2.2204e-16; double alphai = 0;
	double alpha = 100000000; int iter = 0;
	bool init = false; int n0 = b0.M();
	double elapsedTime = 0;
	y.Shape(n0, 1); elapsedTime = 0;
	Epetra_SerialDenseMatrix s0;
	vector<int> P(n0);
	Epetra_SerialDenseMatrix vector_x, vector_b;
	Epetra_SerialSymDenseMatrix solverMatrix;

	// Initialize active set
	vector<int> positions; 
#pragma omp parallel for schedule(guided, 16) reduction(mergeI:positions)
	for (int i = 0; i < y0.size(); i++){
		if (y0[i] >= nnlstol) {
			positions.push_back(i);
		}
	}
  
	int counter = 0;
	counter = positions.size();
	w.Reshape(b0.M(), b0.N());

	if (counter == 0) {
#pragma omp parallel for schedule (static, 16) // Always same workload -> static
	for (int x = 0; x < b0.M(); x++) {
		w(x, 0) = -b0(x, 0);
	}
    
    init = false;
	} else {
#pragma omp parallel for schedule (static, 16) // Always same workload -> static
		for (int i = 0; i < counter; i++) {
			P[i] = positions[i];
		}
		init = true;
	}

	s0.Shape(n0, 1);  // Replacement for s
	bool aux1 = true, aux2 = true;
  
	// New searching algorithm
	vector<double> values, newValues;
	vector<int> poss, newPositions;
	double minValue = w(0, 0); int minPosition = 0;

	while (aux1 == true) {
		if (filename == "sup8.dat"){
#pragma omp parallel
					{
						double minVP = w(0,0); int minPosP = 0;
#pragma omp parallel for schedule (guided, 16) // Guided and static seem to be even but guided makes more sense
						for (int i = 0; i < w.N(); i++){
							if (minVP > w(i,0)) {
								minVP = w(i, 0); minPosP = i;
							}
						}
#pragma omp critical
						{
							minValue = minVP; minPosition = minPosP;
						}
					}
				} else {
					// This is slightly slower than the optimal one. So far at least. Should have a bit better scaling.
					// @{
#pragma omp parallel for schedule(static, 16) reduction(mergeI:poss) reduction(mergeD:values)
					for(int i = 0; i < w.M(); i++){
						values.push_back(w(i, 0)); poss.push_back(i);
					}
				  
					// Get all values smaller than initial one
					while(values.size() > 1) {
#pragma omp parallel for schedule(dynamic, 16) reduction(mergeI:newPositions) reduction(mergeD:newValues)
						for (int i = 1; i < values.size(); i++){
							if (values[i] < values[0]){
								newValues.push_back(values[i]); newPositions.push_back(poss[i]);
							}
						}
					  
						if (newValues.size() == 0){
							newValues.push_back(values[0]); newPositions.push_back(poss[0]);
						}	
					  
						values = newValues; poss = newPositions;
						newValues.clear(); newPositions.clear();
					}
					minValue = values[0]; minPosition = poss[0];
					values.clear(); poss.clear();
				}
				
		// }
	  
		if (((counter == n0) || (minValue > -nnlstol) || (iter >= maxiter)) &&
				(init == false)) {
			aux1 = false;
		} else {
			if (init == false) {
				P[counter] = minPosition;
				counter += 1;
			} else {
				init = false;
			}
		}

		int j = 0;
		aux2 = true;
		while (aux2 == true) {
			iter++;
			vector_x.Shape(counter, 1); vector_b.Shape(counter, 1);
			solverMatrix.Shape(counter, counter);
      
#pragma omp parallel for schedule (static, 16) // Always same workload -> Static!
			for (int x = 0; x < counter; x++) {
				vector_b(x, 0) = b0(P[x], 0);

				for (int z = 0; z < counter; z++) {
					if (x >= z)
						solverMatrix(x, z) = matrix(P[x], P[z]);
					else
						solverMatrix(z, x) = matrix(P[x], P[z]);
				}
			}

			// Measure time LinearSolve needs
			auto start = std::chrono::high_resolution_clock::now();
			
			LinearSolve(solverMatrix, vector_x, vector_b);
			
			auto finish = std::chrono::high_resolution_clock::now();
			elapsedTime += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();

#pragma omp parallel for schedule (static, 16) // Always same workload -> Static!
			for (int x = 0; x < counter; x++) {
				s0(P[x], 0) = vector_x(x, 0);
			}

			bool allBigger = true;
#pragma omp parallel for schedule (static, 16) // Always same workload -> Static!
			for (int x = 0; x < counter; x++) {
				if (s0(P[x], 0) < nnlstol) {
					allBigger = false;
					// break; // -> Terminate Loop!
				}
			}

			if (allBigger == true) {
				aux2 = false;
#pragma omp parallel for schedule(guided, 16)
				for (int x = 0; x < counter; x++) {
					y(P[x], 0) = s0(P[x], 0);
				}

				// w=A(:,P(1:nP))*y(P(1:nP))-b;
				w.Scale(0.0);
#pragma omp parallel for schedule (dynamic, 16)
				for (int a = 0; a < matrix.M(); a++) {
					w(a, 0) = 0;
					for (int b = 0; b < counter; b++) { 
						w(a, 0) += (matrix(a, P[b]) * y(P[b], 0));
					}
					w(a, 0) -= b0(a, 0);
				}
			} else {
				j = 0;
				alpha = 1.0e8;
        
				// Searching for minimum value with index position
#pragma omp parallel
				{
					int jP = 0;
					double alphaP = alpha;
#pragma omp parallel for schedule(guided, 16) // Even, guided seems fitting
					for (int i = 0; i < counter; i++) {
						if (s0(P[i], 0) < nnlstol) {
							alphai = y(P[i], 0) / (eps + y(P[i], 0) - s0(P[i], 0));
							if (alphai < alphaP) {
								alphaP = alphai;
								jP = i;
							}
						}
					}
#pragma omp critical
					{
						alpha = alphaP;
						j = jP;
					}
				}
        
#pragma omp parallel for schedule (guided, 16)
				for (int a = 0; a < counter; a++)
					y(P[a], 0) = y(P[a], 0) + alpha * (s0(P[a], 0) - y(P[a], 0));
				if (j > 0) {
					// jth entry in P leaves active set
					s0(P[j], 0) = 0;
					P.erase(P.begin() + j);
#pragma omp atomic // Necessary?
					counter -= 1;
				}
			}
		}
	}
	return elapsedTime;
}

/*------------------------------------------*/

int main(int argc, char* argv[]) {
	// Commenting might not cause wrong thread amount on server
	// omp_set_num_threads(6); // 6 seems to be optimal
  
	auto start2 = std::chrono::high_resolution_clock::now(); // Timer for geneeralized time
	int csteps, flagwarm;
	double nu1, nu2, G1, G2, E, G, nu, alpha, H, rnd, k_el, delta, nnodi, to1, E1,
		E2, lato, zref, ampface, errf, sum = 0;
	double Delta = 50;  // only used for debugging
	double elapsedTimeM = 0;
	string randomPath = "sup8.dat";
	SetParameters(E1, E2, csteps, flagwarm, lato, zref, ampface, nu1, nu2, G1, G2,
                E, G, nu, alpha, H, rnd, k_el, delta, nnodi, errf, to1);

	std::cout << "File is " + randomPath << endl;
  
	time_t now = time(0);
	tm *ltm = localtime(&now);
	std::cout << "Time is: " << ltm->tm_hour << ":";
	std::cout << 1 + ltm->tm_min << endl;

	// Meshgrid-Command
	// Identical Vectors/Matricies, therefore only created one here.
	int iter = int (ceil((lato - (delta / 2)) / delta)); // Replacement for "for (double i = delta / 2; i < lato; i = i + delta)"
	vector<double> x(iter);
  
#pragma omp parallel for schedule(static, 16) // Same amount of work -> static
	for (int i = 0; i < iter; i++) {
		x[i] = (delta / 2) + i * delta;
	}
  
	// Setup Topology
	Epetra_SerialDenseMatrix topology, y;
	CreateTopology(topology.N(), topology, randomPath);

	double zmax = 0;
	double zmean = 0;
	int cont = 0;
  
#pragma omp parallel for schedule (guided, 16) reduction(+:zmean) reduction(max:zmax) 
	// Static and Guided seem even but Guided makes more sense
	for (int i = 0; i < topology.M(); i++){
		for (int j = 0; j < topology.N(); j++) {
			zmean += topology(i, j);
			if (topology(i, j) > zmax){
				zmax = topology(i, j);
			}
		}
	}
  
	zmean = zmean / (topology.N() * topology.M());

	vector<double> force0, area0;
	double w_el = 0, area = 0, force = 0;
	int k = 0, n0 = 0;
	std::vector<double> xv0, yv0, b0, x0, xvf, 
	yvf, pf; vector<int> col, row;
	double nf, xvfaux, yvfaux, pfaux;
	Epetra_SerialDenseMatrix A;
	int nf2 = floor(nf);
	
	while (errf > to1 && k < 100) {
		// First predictor for contact set
		// All points, for which gap is bigger than the displacement of the rigid
		// indenter, cannot be in contact and thus are not checked in nonlinear solve
		// @{
		
		// [ind1,ind2]=find(z>=(zmax-(Delta(s)+w_el(k))));
		double value = zmax - Delta - w_el;
		row.clear(); col.clear();
    
#pragma omp parallel for schedule(guided, 16) reduction(mergeI:row) reduction(mergeI:col)
		// Data is even, guided makes more sense
		for (int i = 0; i < topology.N(); i++){
			for (int j = 0; j < topology.N(); j++){
				if (topology(i, j) >= value){
					row.push_back(i); col.push_back(j);
				}
			}
		}
		
		n0 = col.size();
    
		// @{
		xv0.clear(); xv0.resize(n0);
		yv0.clear(); yv0.resize(n0);
		b0.clear(); b0.resize(n0);
		// @} Parallelizing slows down program here, so not parallel

#pragma omp for schedule (guided, 16) // Always same workload but testing might be good -> Guided?
		for (int b = 0; b < n0; b++){
			try{
				xv0[b] = x[col[b]];
			} catch (const std::exception& e){}
		}
    
#pragma omp parallel for schedule (guided, 16) // Same
		for (int b = 0; b < n0; b++) {
			try{
				yv0[b] = x[row[b]];
			} catch (const std::exception& e){}
		}
    
#pragma omp parallel for schedule (guided, 16) // Same
		for (int b = 0; b < n0; b++) {
			try{
				b0[b] = Delta + w_el - (zmax - topology(row[b], col[b]));
			} catch (const std::exception& e) {}
		}

		int err = A.Shape(xv0.size(), xv0.size());

		// Construction of the Matrix H = A
		SetUpMatrix(A, xv0, yv0, delta, E, n0, k);


        // Second predictor for contact set
        // @{
        
		Epetra_SerialDenseMatrix xv0t, yv0t, xvft, yvft, pft, x0temp; // Temporary variables for warmup
		if (flagwarm == 1 && k > 1) {
			// x0=warm_x(xv0(1:n0(k),k),yv0(1:n0(k),k),xvf(1:nf(k-1),k-1),yvf(1:nf(k-1),k-1),pf(1:nf(k-1),k-1));
			xv0t.Shape(1, n0); yv0t.Shape(1, n0); xvft.Shape(1, nf2);
			yvft.Shape(1, nf2); pft.Shape(1, nf2);
        
#pragma omp parallel for schedule (static, 16) // Always same workload -> Static
			for (int i = 0; i < n0; i++) {
				xv0t(0, i) = xv0[i];
				yv0t(0, i) = yv0[i];
			}
    	
#pragma omp parallel for schedule (static, 16) // Always same workload -> Static
			for (int j = 0; j < nf2; j++) {
				xvft(0, j) = xvf[j];
				yvft(0, j) = yvf[j];
				pft(0, j) = pf[j];
			}
    	
			x0temp = Warmstart(xv0t, yv0t, xvft, yvft, pft);
    	
#pragma omp parallel for schedule (static, 16) // Always same workload -> Static
			for (int i = 0; i < x0temp.N(); i++){
				x0[i] = x0temp(i, 1);
			}
		} else {
			if (b0.size() > 0) {
				// x0.Shape(b0.N(), 1);
				x0.resize(b0.size());
			}
		}
		// }
    
		// {
		x0.clear(); x0.resize(b0.size());
		Epetra_SerialDenseMatrix b0new;
		b0new.Shape(b0.size(), 1);
		// } Parallel region makes around this makes program slower
    
#pragma omp parallel for schedule (static, 16) // Always same workload -> Static!
		for (int i = 0; i < b0.size(); i++) {
			b0new(i, 0) = b0[i];
		}

		Epetra_SerialDenseMatrix w;
    
		elapsedTimeM += NonlinearSolve(A, b0new, x0, w, y) * pow(10, -3);  // y -> sol, w -> wsol; x0 -> y0
    
		// Compute residial
		// @{
		Epetra_SerialDenseMatrix res1;
		if (A.M() != y.N()) {
			std::runtime_error("Error 1: Matrix dimensions imcompatible");
		}
		res1.Shape(A.M(), A.M());
    
		// res1=A*sol-b0(:,k)-wsol;
#pragma omp parallel for schedule (static, 16) // Always same workload -> Static
		for (int x = 0; x < A.N(); x++) {
			for (int z = 0; z < y.M(); z++) {
				res1(x, 0) += A(x, z) * y(z, 0);
			}
			res1(x, 0) -= b0new(x, 0) + w(x, 0);  // [...]-b0(:,k) - wsol;
		}
		// }

		// Compute number of contact node
		// @{
    
		// @{
		xvf.clear(); xvf.resize(y.M());
		yvf.clear(); yvf.resize(y.M());
		pf.resize(y.M());
		cont = 0;
		// @} Parallelizing this slows down program, so removed it.
  
#pragma omp for schedule(guided, 16)
		for (int i = 0; i < y.M(); i++) {
			if (y(i, 0) != 0) {
#pragma omp critical
				{
					xvf[cont] = xv0[i];
					yvf[cont] = yv0[i];
					pf[cont] = y(i, 0);
					cont += 1;
				}
			}
		}
    
		nf = cont;
		// }

		// Compute contact force and contact area
		// @{
		force0.push_back(0);
    
		sum = 0; iter = ceil(nf);
#pragma omp parallel for schedule (static, 16) reduction(+: sum) // Always same workload -> Static!
		for (int i = 0; i < iter; i++){
			sum += pf[i];
		}
		force0[k] += sum;
		area0.push_back(nf * (pow(delta, 2) / pow(lato, 2)) * 100);
		w_el = force0[k] / k_el;
    
		// }

		// Compute error due to nonlinear correction
		// @{
		if (k > 0) {
			errf = abs(force0[k] - force0[k - 1]) / force0[k];

			// errw(k) = abs((w_el0(k+1)-w_el0(k))/w_el0(k+1));
			// It appears that this is only a debugging variable without any uses,
			// therefore im not gonna implement this here.
		}
		k += 1;
		// }
	}
	
	// @{

	force = force0[k - 1];
	area = area0[k - 1];

	// }
	// }

	double sigmaz = force / pow(lato, 2); // Mean pressure
	double pressz = sigmaz; // Pressure unit per depth
	cout << "k= " << k << " nf= " << nf << endl;
	cout << "Force= " << force << endl;
	cout << "area= " << area << endl;
	cout << "Mean pressure is:" + std::to_string(sigmaz) +
              " ; pressure unit per depth is:" + std::to_string(pressz) +
              " . \n";
	if (abs(sigmaz - 0.130720) > to1) {
		std::runtime_error("Differenz ist zu groß!");  // for nn=2
	}
	// if (abs(sigmaz - 0.246623) > to1)
	//   std::runtime_error("Differenz ist zu groß!");  // for nn=5
  
	auto finish = std::chrono::high_resolution_clock::now();
	double elapsedTimeG = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start2).count() * pow(10, -3);
	double elapsedTime = (elapsedTimeG - elapsedTimeM);
	
	std::cout << "LinearSolve: elapsed time is: " + to_string(elapsedTimeM) + "s." << endl;
	std::cout << "Parallel elapsed time is: " + to_string(elapsedTime) + "s." << endl;
	std::cout << "General elapsed time is: " + to_string(elapsedTimeG) + "s." << endl;
  
	// Get number of threads
	int thread_amount = 99; // error number
  
	// Thread amount is only correct if used in parallel region
#pragma omp parallel for schedule (static, 16) // Scheduling type irrelevant, only one iteration
	for (int i = 0; i < 1; i++){
		thread_amount = omp_get_num_threads();
	}
  
	writeToFile(to_string(elapsedTime), randomPath, thread_amount);
}
