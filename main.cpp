#include <omp.h>
#include <unistd.h>
#include <cmath>
#include <cstdio>	
#include <fstream> 
#include <iostream> 
#include <string>
#include <vector>
#include <Epetra_SerialSpdDenseSolver.h>
#include <Epetra_SerialSymDenseMatrix.h>
#include <chrono> 
#include <ctime>
#include <memory>
#include <jsoncpp/json/json.h>
using namespace std;

#include "topology.h"
#include "topologyfactory.h"

// Declaration for std::vector<int> reduction in parallel loops.
#pragma omp declare reduction(mergeI:std::vector<int>:omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

// Declaration for std::vector<double> reduction in parallel loops.
#pragma omp declare reduction(mergeD:std::vector<double>:omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

void SetParameters(double& E1, double& E2,
                   double& lato, double& nu1,
                   double& nu2, double& G1, double& G2, double& E,
                   double& alpha,
                   double& k_el, double& delta, double& nnodi, double& errf,
                   double& to1, double& Delta, string& zfilePath, int& n, string& jsonFileName, bool& rmg_flag, double& Hurst, bool& rand_seed_flag) {
  
  
  Json::Value parameterlist;   // will contain the root value after parsing.
  ifstream stream(jsonFileName, std::ifstream::binary);
  stream >> parameterlist; 

  rmg_flag = parameterlist["rmg_flag"].asBool();
  rand_seed_flag = parameterlist["rand_seed_flag"].asBool();
  zfilePath = parameterlist["z_file_path"].asString();
  E1 = parameterlist["parameters"]["material_parameters"]["E1"].asDouble();
  E2 = parameterlist["parameters"]["material_parameters"]["E2"].asDouble();
  nu1 = parameterlist["parameters"]["material_parameters"]["nu1"].asDouble();
  nu2 = parameterlist["parameters"]["material_parameters"]["nu2"].asDouble();
  G1 = E1 / (2 * (1 + nu1));
  G2 = E2 / (2 * (1 + nu2));
  E = 1 / ((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2) / E2));
  vector<double> alpha_con{0.778958541513360, 0.805513388666376,
                           0.826126871395416, 0.841369158110513,
                           0.851733020725652, 0.858342234203154,
                           0.862368243479785, 0.864741597831785};
  n = parameterlist["parameters"]["geometrical_parameters"]["n"].asInt();
  Hurst = parameterlist["parameters"]["geometrical_parameters"]["H"].asDouble(); // Hurst component
  alpha = alpha_con[n - 1];
  lato = parameterlist["parameters"]["geometrical_parameters"]["lato"].asDouble();  // Lateral side of the surface [micrometers]
  k_el = lato * E / alpha; 
  delta = lato / (pow(2, n) + 1);
  nnodi = pow(pow(2, n + 1), 2);
  errf = parameterlist["parameters"]["geometrical_parameters"]["errf"].asDouble();
  to1 = parameterlist["parameters"]["geometrical_parameters"]["tol"].asDouble();
  Delta = parameterlist["parameters"]["geometrical_parameters"]["Delta"].asDouble();
}


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

// Parallel not possible, since multiple accesses on the same object.
// In addition, parallel programming not supported for LinearSolve.

/*------------------------------------------*/
Epetra_SerialDenseMatrix Warmstart(Epetra_SerialDenseMatrix xv0, Epetra_SerialDenseMatrix yv0, Epetra_SerialDenseMatrix& xvf,
    Epetra_SerialDenseMatrix& yvf, Epetra_SerialDenseMatrix& pf) {
    Epetra_SerialDenseMatrix x0; x0.Shape(xv0.N(), 1);
    Epetra_SerialDenseMatrix combinedMatrix;
    combinedMatrix.Shape(2 * xv0.N(), xv0.N());
    // matfin = [xv0, yv0]
#pragma omp parallel for schedule (static, 16) // Always same workload -> Static
    for (int i = 0; i < xv0.N(); i++) {
        for (int j = 0; j < xv0.N(); j++) {
            combinedMatrix(i, j) = xv0(i, j);
        }
    }
    
#pragma omp parallel for schedule (static, 16) // Always same workload -> Static
    for (int i = xv0.N(); i < yv0.N() * 2; i++) {
        for (int j = 0; j < yv0.N(); j++) {
            combinedMatrix(i, j) = yv0(i - xv0.N(), j);
        }
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
void writeForceToFile(Epetra_SerialDenseMatrix& y, string pathName){

  int i;

  std::size_t botDirPos = pathName.find_last_of("/") + 1;
  // get directory
  std::string dir = pathName.substr(0, botDirPos);
  // get file name
  std::string dat_file = pathName.substr(botDirPos, pathName.length());

  std:string file = dir + "result_force_" + dat_file;

	ofstream outfile;
	outfile.open(file, std::ofstream::trunc);

  for (int i = 0; i < y.M(); i++) {
    if (y(i, 0) != 0) {
      outfile << y(i,0) << ";";
    }
  }
}
/*------------------------------------------*/
Epetra_SerialDenseMatrix Warmstart2(Epetra_SerialDenseMatrix xv0, Epetra_SerialDenseMatrix yv0, Epetra_SerialDenseMatrix& xvf,
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
void NonlinearSolve(Epetra_SerialDenseMatrix& matrix,
                    Epetra_SerialDenseMatrix& b0, std::vector<double>& y0,
                    Epetra_SerialDenseMatrix& w, Epetra_SerialDenseMatrix& y) {
  // matrix -> A, b0 -> b, y0 -> y0 , y -> y, w-> w; nnstol, iter, maxiter ->
  // unused
  double nnlstol = 1.0000e-08;
  double maxiter = 10000;
  double eps = 2.2204e-16;
  double alphai = 0;
  double alpha = 100000000;
  int iter = 0;
  bool init = false;
  int n0 = b0.M();
  y.Shape(n0, 1);
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
	  // [wi,i]=min(w);
	  // This is slightly slower than the optimal one. So far at least. Should have a bit better scaling.
	  // @{
#pragma omp parallel for schedule(static, 16) reduction(mergeI:poss) reduction(mergeD:values)
	  for(int i = 0; i < w.M(); i++){
		  values.push_back(w(i, 0)); poss.push_back(i);
	  }
	  
	  // Get all values bigger than initial one
	  while(values.size() > 1){
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

      LinearSolve(solverMatrix, vector_x, vector_b);

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
        
        // TODO: WHAT BELONGS TO THIS LOOP?????????????????????
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
}

/*------------------------------------------*/

int main(int argc, char* argv[]) {
	omp_set_num_threads(6); // 6 seems to be optimal
  
	auto start = std::chrono::high_resolution_clock::now();
	int csteps, flagwarm, n;
	double nu1, nu2, G1, G2, E, alpha, k_el, delta, nnodi, to1, E1,
		E2, lato, errf, sum = 0, Delta;
  bool rmg_flag;
  bool rand_seed_flag;
  double Hurst;
	string zfilePath;

  string jsonFileName = argv[1]; // reading the json file name from the command line

	SetParameters(E1, E2, lato, nu1, nu2, G1, G2,
                E, alpha, k_el, delta, nnodi, errf, to1, Delta, zfilePath, n, jsonFileName, rmg_flag, Hurst, rand_seed_flag);

  
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
  int N = pow(2,n);
  topology.Shape(N+1,N+1);

  std::shared_ptr<TopologyGeneration> surfacegenerator;

  CreateSurfaceObject(n,Hurst,rand_seed_flag,zfilePath,rmg_flag,surfacegenerator); // creating the correct surface object

  surfacegenerator->GetSurface(topology);

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
    	
    	x0temp = Warmstart2(xv0t, yv0t, xvft, yvft, pft);
    	
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
    NonlinearSolve(A, b0new, x0, w, y);  // y -> sol, w -> wsol; x0 -> y0

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

  // Mean pressure
  double sigmaz = force / pow(lato, 2);
  // Pressure unit per depth
  double pressz = sigmaz;
  cout << "k= " << k << " nf= " << nf << endl;
  cout << "Force= " << force << endl;
  cout << "area= " << area << endl;
  cout << "Mean pressure is:" + std::to_string(sigmaz) +
              " ; pressure unit per depth is:" + std::to_string(pressz) +
              " . \n";
 // if (abs(sigmaz - 0.130720) > to1)
 //   std::runtime_error("Differenz ist zu groß!");  // for nn=2
   if (abs(sigmaz - 0.246623) > to1)
     cout << "Differenz ist zu groß!" << std::endl;  // for nn=5
  
  
  auto finish = std::chrono::high_resolution_clock::now();
  double elapsedTime2 = std::chrono::duration_cast<std::chrono::seconds>(finish - start).count();
  std::cout << "Elapsed time is: " + to_string(elapsedTime2) + "s." << endl;
  
  // Get number of threads
  int thread_amount = 99; // error number
#pragma omp parallel for schedule (static, 16) // Scheduling type irrelevant, only one iteration
  for (int i = 0; i < 1; i++){
	  thread_amount = omp_get_num_threads();
  }
  
  writeForceToFile(y, zfilePath);
}
