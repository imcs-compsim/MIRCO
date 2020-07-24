#include <cmath> //pow
#include <cstdio>
#include <iostream> //ifstream
#include <fstream> //ifstream
#include <string> //std::to_string, std::stod
#include <vector>
#include "include/Epetra_SerialSpdDenseSolver.h"
#include "include/Epetra_SerialSymDenseMatrix.h"
#include <unistd.h>

using namespace std;

/**
* Sets up any given parameter.
* Cross-checked for functionality. Should work as intented.
*/
void SetParameters(int& E1, int& E2, int& csteps, int& flagwarm, int& lato, int& zref, int& ampface,
    double& nu1, double& nu2, double& G1, double& G2, double& E, double& G, double& nu, double& alpha,
    double& H, double& rnd, double& k_el, double& delta, double& nnodi, int& errf, float& to1) {
    E1 = 1; E2 = 1;
    nu1 = 0.3; nu2 = 0.3;
    G1 = E1 / (2 * (1 + nu1));
    G2 = E2 / (2 * (1 + nu2));
    E = 1/((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2) / E2));
    G = 1/((2 - nu1) / (4 * G1) + (2 - nu2) / (4 * G2));
    nu = E / (2 * G) - 1;

    vector<double> alpha_con { 0.778958541513360, 0.805513388666376, 0.826126871395416, 0.841369158110513,
    0.851733020725652, 0.858342234203154, 0.862368243479785, 0.864741597831785 };
    int nn = 2; // Matrix sent has the parameter nn=2!
    alpha = alpha_con[nn-1];
    csteps = 1;
    ampface = 1;
    flagwarm = 0;
    lato = 1000; // Lateral side of the surface [micrometers]
    H = 0.1; // Hurst Exponent (D = 3 - H)
    rnd = 95.0129;
    zref = 50; // Reference for the Scaling, former value = 25
    k_el = lato * E / alpha;
    delta = lato / (pow(2, nn)+1);
    nnodi = pow(pow(2, nn + 1), 2);

    errf = 100000000;
    to1 = 0.01;
}


/*------------------------------------------*/

/**
* Readin Matrix from a file. Elements have to be separated by a ';'.
* Cross-checked for functionality. Should work as intended.
*/
void CreateTopology(int systemsize, Epetra_SerialDenseMatrix& topology, string filePath) {
    // Readin for amount of lines -> dimension of matrix
    ifstream reader(filePath);
    string blaLine;
    int dimension = 0;
    while (getline(reader, blaLine)) { dimension += 1; }
    reader.close();
    topology.Shape(dimension, dimension);
    int lineCounter = 0;
    float elements[264];
    int position = 0;
    try {
        ifstream stream(filePath);
        string line;
        bool negative = false;
        bool lastValue = false;
        while (std::getline(stream, line)) {

            // Split up Values into Double-Array
            int separatorPosition = 0;
            lineCounter += 1; // Has to happen here, since baseline value is 0.

            int separatorAmount = count(line.begin(), line.end(), ';');

            for (int i = 0; i < separatorAmount; i--) { // prevent duplication of values!
                separatorPosition = line.find_first_of(';');
                string container = line.substr(0, separatorPosition);
                line = line.substr(separatorPosition + 1, line.length());
                if (line.length() < 2) { i = -1; } // exit condition to avoid error's, non-scientific!
                if (container == "" || container == ";") { i = -1; } // exit condition to avoid duplication of values!

                if (container.substr(0, 1) == "-") { // Substring Double-Value!
                    negative = true;
                    container = container.substr(1, container.length() - 1);
                }
                double value = stod(container);
                if (negative == true) { value = value * (-1); }

                topology(lineCounter, i + 1); // +1 has to happen, since baseline value is 0.
                // elements[position] = value; happened before

                // Reset Conditions
                negative = false;
                position += 1;
            }
        }
        stream.close();
    }
    catch (const std::exception& e) {
        std::cout << e.what(); //Fatal Error, catch it!
    }
}

/*------------------------------------------*/

void SetUpMatrix(Epetra_SerialSymDenseMatrix& A,std::vector<double> xv0, std::vector<double> yv0,
		double delta, double E, int systemsize, int k) {
	
	// DEBUG
	cout << "Init done. \n";
	
    int r;
    double pi = atan(1) * 4;
    double raggio = delta / 2;
    double C = 1 / (E * pi * raggio);

    cout << "systemsize= " << systemsize << endl;
    cout << "size of A = " << A.N() << endl;

    for (int i = 0; i < systemsize; i++) {
        A(i, i) = 1 * C;
    }
    
    // DEBUG
    cout << "Part 2 done. \n";

    for (int i = 0; i < systemsize; i++) {
        for (int j = 0; j < i; j++) {
            r = sqrt(pow((xv0[j] - xv0[i]),2) + pow((yv0[j] - yv0[i]),2));
            A(i, j) = C * asin(raggio / r);
        }
    }
    
    // DEBUG
    cout << "Part 3 done. \n";
}

/*------------------------------------------*/

void LinearSolve(Epetra_SerialSymDenseMatrix& matrix,
    Epetra_SerialDenseMatrix& vector_x,
    Epetra_SerialDenseMatrix& vector_b) {
    Epetra_SerialSpdDenseSolver solver;
    int err = solver.SetMatrix(matrix);
    if (err != 0) { std::cout << "Error setting up matrix solver (1)"; }

    err = solver.SetVectors(vector_x, vector_b);
    if (err != 0) { std::cout << "Error setting up maxtix solver (2)"; }

    err = solver.Solve();
    if (err != 0) { std::cout << "Error setting up matrix solver (3)"; }
    std::cout << vector_x << std::endl;
}

/*------------------------------------------*/

void NonlinearSolve(Epetra_SerialSymDenseMatrix& matrix, Epetra_SerialDenseMatrix& b0,
    std::vector<double> y0, Epetra_SerialDenseMatrix& w, Epetra_SerialDenseMatrix& y) {
    // matrix -> A, b0 -> b, y0 -> y0 , y -> y, w-> w; nnstol, iter, maxiter -> unused
    double nnlstol = 1.0000e-08;
    double maxiter = 10000;
    bool init = false;
    int n0 = b0.N() * b0.M();
    y.Shape(n0, 1);
    vector<int> P(n0);
    Epetra_SerialDenseMatrix vector_x, vector_b;
    Epetra_SerialSymDenseMatrix solverMatrix;
    
    // DEBUG
    cout << "Pre-Setup done. \n";

    // Initialize active set
    vector<int> positions;
    int counter = 0;
    for (int i = 0; i < y0.size(); i++) {
        if ((y0[i] == nnlstol) || (y0[i] > nnlstol)) {
            positions.push_back(i);
            counter += 1;
        }
    }
    
    // DEBUG
    cout << "Active set initialized. \n";
    
    // Avoid .N() = 0 or .M() = 0, completely wrecks code otherwise
    w.Reshape(b0.M(), b0.N());
    
    cout << "After shaping w\n";
    cout << "Counter= " << counter << endl;
    if (counter == 0) {
            for (int x = 0; x < b0.M(); x++) {
                w(x, 0) = -b0(x, 0);
            }
    } else {
        for (int i = 0; i < counter; i++) {
            P[i]=positions[i];
        }
        init = true;
    }

    cout << "After if-else \n";

    Epetra_SerialDenseMatrix s0; s0.Shape(counter+1, 1); // Replacement for s
    bool aux1 = true, aux2 = true;
    
    // DEBUG
    cout << "While-Loop starting. \n";
    
    while (aux1 == true) {
        // [wi,i]=min(w);
        int minValue = w(1, 1), minPosition = 0;
        for (int i = 0; i < w.M(); i++) {
            if (minValue > w(i, 0)) {
                minValue = w(i, 0);
                minPosition = i;
            }
        }

        if ((counter == n0) || ((minValue > -nnlstol) && (init == false))) {
            aux1 = false;
        } else {
            if (init == false) {
                // Index #i enter active index
                counter += 1;
                P.push_back(minPosition);
            }
        }
        
        // DEBUG
        cout << "Part 1 done. \n";
        
        int j = 0;
        double eps = 2.2204e-16; int alphai = 0, alpha = 100000000, a = 0;
        while (aux2 == true) {
        	
        	// DEBUG
        	cout << "Counter has value: " + to_string(counter) + " .\n";
        	
        	// Avoid errors
        	if (counter < 1){ counter = 1; } 
        	
            vector_x.Shape(counter, 1);
            vector_b.Shape(counter, 1);
            solverMatrix.Shape(counter, counter);
            
            cout << "Shape done. \n";
            
            // TODO: Why does this not work???

            for (int x = 0; x < counter; x++) {
            	try {
            		vector_b(x, 0) = b0(P[x], 0);
            	} catch (const std::exception& e) { // Catch errors if P[x - 1] is not initialized
            		vector_b(x, 0) = 0;	// MatLab standart value if not initialized
            	}
            	cout << "vector_b(1, 1) is: " + to_string(vector_b(1, 1)) +  " .\n";
            	// Cant assign values to vector_b
            	printf("After catch phase\n");
                for (int z = 0; z < counter; z++) {
                	try {
                		solverMatrix(x, z) = matrix(P[x], P[z]);
                	} catch (const std::exception& e){ // Catch errors if P[x - 1], P[z - 1] are not initialized
                		solverMatrix(x, z) = 0; // MatLab standart value if not initialized
                	}
                }
            }
            cout << "solverMatrix(1, 1) is: " + to_string(solverMatrix(1, 1)) + " \n";
            // Cant assign values to solverMatrix
            
            // DEBUG
            cout << "Part 2 done. \n";
            
            // Linear solve
            LinearSolve(solverMatrix, vector_x, vector_b);
            
            cout << "Linear solve done. \n";

            cout << "s0.M= " << s0.M() << " vector_b.M= " << vector_b.M() << " P[0]= " << P[0] << endl;
            for (int x = 0; x < counter; x++) {
                s0(P[x], 0) = vector_b(x, 0);
            }

            cout << "Filled s0 \n";

            bool allBigger = true;
            for (int x = 0; x < counter; x++) {
                if (s0(P[x], 0) < nnlstol) { allBigger = false; }
            }

            cout << "After allbigger \n";
            if (allBigger == true) {
                aux2 = false;

                if (matrix.M() != y.N()) { std::runtime_error("Fehler 2: Ungültige Matrixdimension! \n"); }

                cout << "Before Matrixproduct \n";

                // w=A(:,P(1:nP))*y(P(1:nP))-b;
                for (int a = 0; a < matrix.M(); a++) {
                    w(a, 0) = 0;
                    for (int b = 0; b < matrix.N(); b++) {
                        w(a, 0) += (matrix(a, P[b]) * y(P[b], 0)) - b0(a, 0);
                    }
                }
                
                cout << "After Matrixproduct \n";

                aux1 = true; // Exit condition
            } else {
                for (int i = 0; i < counter; i++) {
                    if (s0(P[i], 0) < nnlstol) {
                        alphai = y(P[i], 0) / (eps + y(P[i], 0) - s0(P[i], 0));
                        if (alphai < alpha) {
                            alpha = alphai;
                            j = 1;
                        }
                    }
                }
            }

            cout << "Before counter while \n";
            while (a < counter) {
                a += 1;
                y(P[a], 0) = y(P[a], 0) + alpha * (s0(P[a], 0) - y(P[a], 0));
            }

            if (j > 0) {
                // jth entry in P leaves active set
                s0(P[j], 0) = 0;
                vector<int> P2 = P;
                P2.push_back(0);
                for (int i = 0; i < j; i++) {
                    P2.push_back(P[i]);
                }
                for (int i = j; i < (counter - 1); i++) {
                    P2.push_back(P[i + 1]);
                }
                P = P2;
                P[counter] = 0;
                counter -= 1;
            }
        }
    }
}
/*------------------------------------------*/

int main(int argc, char* argv[]) {
    int E1, E2, csteps, flagwarm, lato, zref, ampface, errf;
    double nu1, nu2, G1, G2, E, G, nu, alpha, H, rnd, k_el, delta, nnodi;
    float to1;
    SetParameters(E1, E2, csteps, flagwarm, lato, zref, ampface, nu1, nu2, G1, G2, E, G, nu, alpha, H, rnd, k_el, delta, nnodi, errf, to1);
    
    // Meshgrid-Command
    // Identical Vectors/Matricies, therefore only created one here.
    vector<double> x;
    for (int i = delta / 2; i < (lato - delta / 2); i = i + delta) { x.push_back(i); }

    // Setup Topology
    string randomPath = "sup2.dat";
    Epetra_SerialDenseMatrix topology, y;
    CreateTopology(topology.N(), topology, randomPath);


    double zmax = 0;
    double zmean = 0;
    zmean = topology.NormOne() / pow(topology.N(), 2);
    zmax = topology.NormInf();

    double Delta = 50; // TODO only used for debugging

    vector<double> force0, area0, w_el0;
    force0.push_back(0); w_el0.push_back(0);
    double w_el, area, force;
    int k = 0;
    int n0;
    std::vector<double> xv0, yv0, b0, x0,  xvf, yvf, pf; // x0: initialized in Warmstart!
    double nf, xvfaux, yvfaux, pfaux;
    
    Epetra_SerialSymDenseMatrix A;

    // DEBUG
    cout << "While-Loop started. \n";
    
    // So far so good
    
    while (errf > to1) {
        k += 1;
        
        // First predictor for contact set
        // All points, for which gap is bigger than the displacement of the rigid indenter, cannot be in contact and thus are not checked in
        // nonlinear solve
        // @{
        
        // [ind1,ind2]=find(z>=(zmax-(Delta(s)+w_el(k))));
        vector<int> col, row;
        double value = zmax - Delta + w_el;
        for (int i = 0; i < topology.N(); i++) {
        	for (int j = 0; j < topology.N(); j++) {
        		if (topology(i, j) >= value) {
        			row.push_back(i);
        			col.push_back(j);
        			// counter += 1;
        		}
        	}
        }
        n0 = col.size();
        

        // Works until this point
        

            for (int b = 0; b < n0; b++) {
                    xv0.push_back(x[row[b]]);
            }

            for (int b = 0; b < n0; b++) {
                    yv0.push_back(x[col[b]]);
            }

            for (int b = 0; b < n0; b++) {
            	b0.push_back( Delta + w_el - (zmax - topology(row[b], col[b])));
            }


        // }
        
        // DEBUG
        cout << "First predictor done. \n \n";
        
        cout << "n0 = " << n0 << " .\n";
        cout << "k = " << k << " .\n";
        cout << "xv0.M= " << xv0.size() << " und xv0.N()= " << xv0.size() << std::endl;
        cout << "Vor shape von A mit size = " << A.N() << endl;

        int err = A.Shape(xv0.size());

        cout << "Nach shape von A" << endl;

        // Construction of the Matrix H = A
		SetUpMatrix(A, xv0, yv0, delta, E, n0, k);
        
        // DEBUG
        cout << "Matrix set up. \n";

        // Second predictor for contact set
        // @{
        
        	
        	// DEBUG
        	cout << "Warmstart inactive. \n";
        	
        	    x0.clear();
        		x0.resize(b0.size());

        
        // DEBUG
        cout << "Second predictor done. \n";
        cout << "Size b0 = " << b0.size() << endl;

        Epetra_SerialDenseMatrix b0new;
        b0new.Shape(b0.size(), 1);

        cout << "After shaping b0new. \n";

        for (int i = 0; i < b0.size(); i++) {
            b0new(i, 0) = b0[i];
        }

        cout << "After Creation of b0new. \n";

        Epetra_SerialDenseMatrix w;
        NonlinearSolve(A, b0new, x0, w, y); // y -> sol, w -> wsol; x0 -> y0

        // DEBUG
        cout << "Nonlinear solve complete. \n"; // TODO: Error in nnls!

        // Compute residial
        // @{
        Epetra_SerialDenseMatrix res1;
        if (A.M() != y.N()) { std::runtime_error("Error 1: Matrix dimensions imcompatible"); }
        res1.Shape(A.N(), y.M());
        // res1=A*sol-b0(:,k)-wsol;
        // Should work now.
        int sum = 0;
        for (int x = 0; x < A.N(); x++) {
            for (int z = 0; z < y.M(); z++) {
                sum += A(x, z) * y(x, 1);
            }
            res1(x, 0) = sum - b0new(x, 0) - w(x, 0); // [...]-b0(:,k) - wsol;
            sum = 0;
        }
        // }

        // Compute number of contact nodes
        // @{
        int cont = 0;
        xvf.clear(); yvf.clear();
        xvf.resize(A.N()); yvf.resize(A.N()); pf.resize(A.N());
        for (int i = 0; i < A.N(); i++) {
            if (A(i, 1) != 0) {
                cont += 1;
                xvf[cont] = xv0[i];
                yvf[cont] = yv0[i];
                pf[cont] = A(i, 0);
            }
        }
        nf = cont;
        // }

        // Compute contact force and contact area
        // @{
        force0.push_back(0);
        for (int i = 0; i < nf; i++) {
            force0[k] = force0[k] + pf[i];
            area0.push_back(nf * (pow(delta, 2) / pow(lato, 2)) * 100);
        }
        w_el0.push_back(force0[k] / k_el);
        // }

        // Compute error because of nonlinear correction
        // @{
        if (k > 1) {
            // errf(k) = abs((force0(k)-force0(k-1))/force0(k));
            errf = (force0[k] - force0[k - 1]) / force0[k];
            if (errf < 0) { errf = (-1) * errf; }

            // errw(k) = abs((w_el0(k+1)-w_el0(k))/w_el0(k+1));
            // It appears that this is only a debugging variable without any uses, therefore im not gonna implement this here.
        }
        // }
    }

    // @{

    force = force0[k];
    area = area0[k];
    w_el = w_el0[k];

    // }
    // }

    // Mean pressure
    double sigmaz = force / pow(lato, 2);
    // Pressure unit per depth
    double pressz = sigmaz;
    cout << "Mean pressure is:" + std::to_string(sigmaz) + " ; pressure unit per depth is:" + std::to_string(pressz) + " . \n";
}
