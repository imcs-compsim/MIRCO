#include <cmath> //pow
#include <cstdio>
#include <iostream> //ifstream
#include <fstream> //ifstream
#include <string> //std::to_string, std::stod
#include <vector>
#include "include/Epetra_SerialSpdDenseSolver.h"
#include "include/Epetra_SerialSymDenseMatrix.h"

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
    E = pow(((1 - pow(nu1, 2)) / E1 + 1 - pow(nu2, 2) / E2), -1);
    G = pow(((2 - nu1) / (4 * G1) + (2 - nu2) / (4 * G2)), -1);
    nu = E / (2 * G) - 1;

    vector<double> alpha_con { 0.778958541513360, 0.805513388666376, 0.826126871395416, 0.841369158110513,
    0.851733020725652, 0.858342234203154, 0.862368243479785, 0.864741597831785 };
    int nn = 2; // Matrix sent has the parameter nn=2!
    alpha = alpha_con[nn];
    csteps = 50;
    ampface = 1;
    flagwarm = 1;
    lato = 1000; // Lateral side of the surface [micrometers]
    H = 0.1; // Hurst Exponent (D = 3 - H)
    rnd = 95.0129;
    zref = 50; // Reference for the Scaling, former value = 25
    k_el = lato * E / alpha;
    delta = lato / pow(2, nn + 1);
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

/* void createTopology(int systemsize, Epetra_SerialDenseMatrix& topology) {
    // Idea: GENERATING a topology, instead of I/O-Topology!
    // t0do: Insert ranmid2d_MP-Method here!

        double scalefactor = zref / (zmax - zmean);
        // z = scalefactor * z;
        for (int i = 0; i < topology.N(); i++) {
            for (int j = 0; j < topology.N(); j++) {
                topology(i, j) = z * topology(i, j);
            }
        }

        // setting minimum heigth to zero
        zmin = zmax;
        for (int i = 0; i < topology.N(); i++) {
            for (int j = 0; j < topology.N(); j++) {
                if (zmin > topology(i, j)) {
                    zmin = topology(i, j);
                }
            }
        }

        // z = z - min(min(z))
        for (int i = 0; i < topology.N(); i++) {
            for (int j = 0; j < topology.N(); j++) {
                topology(i, j) = topology(i, j) - zmin;
            }
        }
        // recalculate mean, max, min
        zmax = 0;
        zmin = topology(0, 0);
        zmean = 0;
        for (int i = 0; i < topology.N(); i++) {
            for (int j = 0; j < topology.N(); j++) {
                zmean = zmean + topology(i, j);
                if (zmax < topology(i, j)) {
                    zmax = topology(i, j);
                }
                if (zmin > topology(i, j)) {
                    zmin = toplogy(i, j);
                }
            }
        }
        zmean = zmean / pow(topology.N(), 2);
} */

Epetra_SerialSymDenseMatrix SetUpMatrix(Epetra_SerialDenseMatrix xv0, Epetra_SerialDenseMatrix yv0, double delta,
    double E, int systemsize, int k) {
    Epetra_SerialSymDenseMatrix A; int r;
    double pi = atan(1) * 4;
    double raggio = delta / 2;
    double C = 1 / (E * pi * raggio);
    A.Shape(systemsize);

    for (int i = 0; i < systemsize; i++) {
        A(i, i) = 1 * C;
    }

    for (int i = 0; i < systemsize; i++) {
        for (int j = 0; j < i; j++) {
            r = (xv0(j, k) - xv0(i, k)) * (yv0(j, k) - yv0(i, k));
            A(i, j) = C * asin(raggio / r);
        }
    }
    return A;
}

/*------------------------------------------*/
Epetra_SerialDenseMatrix Warmstart(Epetra_SerialDenseMatrix xv0, Epetra_SerialDenseMatrix yv0, Epetra_SerialDenseMatrix& xvf,
    Epetra_SerialDenseMatrix& yvf, Epetra_SerialDenseMatrix& pf) {
    Epetra_SerialDenseMatrix x0; x0.Shape(xv0.N(), 1);
    Epetra_SerialDenseMatrix combinedMatrix;
    combinedMatrix.Shape(2 * xv0.N(), xv0.N());
    // matfin = [xv0, yv0]
    for (int i = 0; i < xv0.N(); i++) {
        for (int j = 0; j < xv0.N(); j++) {
            combinedMatrix(i, j) = xv0(i, j);
        }
    }
    for (int i = xv0.N(); i < yv0.N() * 2; i++) {
        for (int j = 0; j < yv0.N(); j++) {
            combinedMatrix(i, j) = yv0(i - xv0.N(), j);
        }
    }

    // loop

    // TODO: vector.push_back()!
    vector<int> index;
    for (int i = 0; i < pf.N(); i++) {
        // ind=find(matfin(:,1)==xvf(i) & matfin(:,2)==yvf(i));
        for (int j = 0; j < xvf.N(); j++) {
            if ((combinedMatrix(j, 1) == xvf(i, 1)) && (combinedMatrix(j, 2) == yvf(i, 1))) {
                index.push_back(j);
            }
        }

        // x0(ind,1)=pf(i);
        for (int y = 0; y < index.size(); y++) {
            x0(y, 1) = pf(y, 1);
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
    if (err != 0) { std::cout << "Error setting up matrix solver (1)"; }

    err = solver.SetVectors(vector_x, vector_b);
    if (err != 0) { std::cout << "Error setting up maxtix solver (2)"; }

    err = solver.Solve();
    if (err != 0) { std::cout << "Error setting up matrix solver (3)"; }
    std::cout << vector_x << std::endl;
}

/*------------------------------------------*/

void NonlinearSolve(Epetra_SerialSymDenseMatrix& matrix, Epetra_SerialDenseMatrix& b0,
    Epetra_SerialDenseMatrix& y0, Epetra_SerialDenseMatrix& w, Epetra_SerialDenseMatrix& y) {
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
    for (int i = 0; i < y0.M(); i++) {
        if ((y0(i, 1) == nnlstol) || (y0(i, 1) > nnlstol)) {
            positions.push_back(i);
            counter += 1;
        }
    }
    
    // DEBUG
    cout << "Active set initialized. \n";
    
    // Avoid .N() = 0 or .M() = 0, completely wrecks code otherwise
    int b0_N = 1, b0_M = 1;
    if (b0.N() > 1){ b0_N = b0.N(); }
    if (b0.M() > 1){ b0_M = b0.N(); }
    w.Shape(b0_N, b0_M);
    
    if (counter == 0) {
        for (int x = 0; x < b0.N(); x++) {
            for (int z = 0; z < b0.M(); z++) {
                w(x, z) = -b0(x, z);
            }
        }
    } else {
        for (int i = 0; i < counter; i++) {
            P.push_back(positions[i]);
        }
        init = true;
    }
    
    Epetra_SerialDenseMatrix s0; s0.Shape(counter, 1); // Replacement for s
    bool aux1 = true, aux2 = true;
    
    // DEBUG
    cout << "While-Loop starting. \n";
    
    while (aux1 == true) {
        // [wi,i]=min(w);
        int minValue = w(1, 1), minPosition = 0;
        for (int i = 0; i < w.M(); i++) {
            if (minValue > w(1, i + 1)) {
                minValue = w(1, i + 1);
                minPosition = i;
            }
        }

        if ((counter == n0) || (minValue > -nnlstol) && (init == false)) {
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
        	
            vector_x.Shape(counter + 1, 1);
            vector_b.Shape(counter + 1, 1);
            solverMatrix.Shape(counter + 1, counter + 1);

            for (int x = 1; x < counter + 1; x++) {
            	try {
            		vector_b(x, 1) = b0(P[x], 1);
            	} catch (const std::runtime_error& e) {
            		vector_b(x, 1) = 0;	// MatLab standart value if not initialized
            	}
            	cout << "vector_b(1, 1) is: " + to_string(vector_b(1, 1)) +  " .\n";
                for (int z = 1; z < counter + 1; z++) {
                	try {
                		solverMatrix(x, z) = matrix(P[x], P[z]);
                	} catch (const std::runtime_error& e){
                		solverMatrix(x, z) = 0; // MatLab standart value if not initialized
                	}
                }
            }
            cout << "Matrix entry is: " + to_string(solverMatrix(1, 1)) + " \n";
            
            // DEBUG
            cout << "Part 2 done. \n";
            
            // Linear solve
            LinearSolve(solverMatrix, vector_x, vector_b);
            
            cout << "Linear solve done. \n";

            for (int x = 0; x < counter; x++) {
                s0(P[x], 1) = vector_b(x, 1);
            }

            bool allBigger = true;
            for (int x = 0; x < counter; x++) {
                if (s0(P[x], 1) < nnlstol) { allBigger = false; }
            }

            if (allBigger == true) {
                aux2 = false;

                if (matrix.M() != y.N()) { std::runtime_error("Fehler 2: Ungültige Matrixdimension! \n"); }

                // w=A(:,P(1:nP))*y(P(1:nP))-b;
                for (int a = 0; a < matrix.N(); a++) {
                    w(a, 1) = 0;
                    for (int b = 0; b < matrix.M(); b++) {
                        w(a, 1) += (matrix(a, P[b]) * y(P[b], 1)) - b0(a, 1);
                    }
                }
                
                aux1 = true; // Exit condition
            } else {
                for (int i = 0; i < counter; i++) {
                    if (s0(P[i], 1) < nnlstol) {
                        alphai = y(P[i], 1) / (eps + y(P[i], 1) - s0(P[i], 1));
                        if (alphai < alpha) {
                            alpha = alphai;
                            j = 1;
                        }
                    }
                }
            }

            while (a < counter) {
                a += 1;
                y(P[a], 1) = y(P[a], 1) + alpha * (s0(P[a], 1) - y(P[a], 1));
            }

            if (j > 0) {
                // jth entry in P leaves active set
                s0(P[j], 1) = 0;
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

    
    // DEBUG
    cout << "Parameters set. \n";
    
    // Meshgrid-Command
    // Identical Vectors/Matricies, therefore only created one here.
    vector<double> x;
    for (int i = delta / 2; i < (lato - delta / 2); i = i + delta) { x.push_back(i); }

    // Setup Topology
    string randomPath = "/home/bartsch/BEM/sub.dat"; // TODO: Change this before debugging!
    Epetra_SerialSymDenseMatrix topology, y;
    CreateTopology(topology.N(), topology, randomPath);
    // TODO: Remove 3rd argument when ranmid2d_MP is implemented!
    
    // DEBUG
    cout << "Topology created. \n";

    double zmax = 0;
    double zmean = 0;
    zmean = topology.NormOne() / pow(topology.N(), 2);

    // Can also use zmatrix = topology.NormInf() and
    // Can also use zmax = topology.NormInf() and
    // zmean = topology.NormOne()/pow(topology.N(), 2)
    for (int i = 1; i < topology.N() + 1; i++) {
        for (int j = 1; j < topology.N() + 1; j++) {
            if (zmax < topology(i, j)) {
                zmax = topology(i, j);
            }
        }
    }

    vector<double> nfaux(csteps);
    double Delta = 39.202067343593399;

    vector<double> force0, area0, w_el0;
    force0.push_back(0); w_el0.push_back(0);
    double w_el, area, force;
    int k = 0;
    vector<int> n0;
    Epetra_SerialDenseMatrix xv0, yv0, b0, x0, nf, xvfaux, yvfaux, pfaux, xvf, yvf, pf; // x0: initialized in Warmstart!
    nf.Shape(csteps, 1); xvfaux.Shape(csteps, 1); yvfaux.Shape(csteps, 1); pfaux.Shape(csteps, 1);
    
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
        			col.push_back(i);
        			row.push_back(j);
        			// counter += 1;
        		}
        	}
        }
        n0.push_back(col.size());
        
        if (k == 1) { // initialize vectors
            xv0.Shape(n0[0], 1); yv0.Shape(n0[0], 1); b0.Shape(n0[0], 1);
        } else { // Initialize matricies with already generated values and an empty row
            Epetra_SerialDenseMatrix a; a.Shape(n0[k], k);
            for (int b = 0; b < n0[k - 1]; b++) {
                for (int c = 0; c < (k - 1); c++) {
                    a(b, c) = xv0(b, c);
                }
            }
            xv0 = a;

            for (int b = 0; b < n0[k - 1]; b++) {
                for (int c = 0; c < (k - 1); c++) {
                    a(b, c) = yv0(b, c);
                }
            }
            yv0 = a;

            for (int b = 0; b < n0[k - 1]; b++) {
                for (int c = 0; c < (k - 1); c++) {
                    a(b, c) = b0(b, c);
                }
            }
            b0 = a;
        }

        for (int i = 0; i < n0[k - 1]; i++) {
            xv0(i, k) = x[row[i]];
            yv0(i, k) = y(row[i], 1);
            b0(i, k) = Delta + w_el - (zmax - topology(row[i], col[i]));
        }
        // }
        
        // DEBUG
        cout << "First predictor done. \n";

        // Construction of the Matrix H = A
        Epetra_SerialSymDenseMatrix A = SetUpMatrix(xv0, yv0, delta, E, n0[k], k);
        
        // DEBUG
        cout << "Matrix set up. \n";

        // Second predictor for contact set
        // @{
        
        Epetra_SerialDenseMatrix xv0t, yv0t, xvft, yvft, pft; // Temporary variables for warmup
        if (flagwarm == 1 && k > 1) {
        	// DEBUG
        	cout << "Warmstart active. \n";
        	
            // x0=warm_x(xv0(1:n0(k),k),yv0(1:n0(k),k),xvf(1:nf(k-1),k-1),yvf(1:nf(k-1),k-1),pf(1:nf(k-1),k-1));
            xv0t.Shape(1, n0[k]); yv0t.Shape(1, n0[k]); xvft.Shape(1, nf(k - 1, 1));
            yvft.Shape(1, nf(k - 1, 1)); pft.Shape(1, nf(k - 1, 1));
            for (int i = 0; i < n0[k]; i++) {
                xv0t(1, i) = xv0(i, k);
                yv0t(1, i) = yv0(i, k);
            }
            for (int i = 0; i < nf(k - 1, 1); i++) {
                xvft(1, i) = xvf(i, k - 1);
                yvft(1, i) = yvf(i, k - 1);
                pft(1, i) = pf(i, k - 1);
            }
            x0 = Warmstart(xv0t, yv0t, xvft, yvf, pft);
        } else {
        	// DEBUG
        	cout << "Warmstart inactive. \n";
        	
        	if (b0.N() > 0) {
        		x0.Shape(b0.N(), 1);
        	} else {
        		x0.Shape(1, 1);
        	}
        }
        // }
        
        // DEBUG
        cout << "Second predictor done. \n";

        Epetra_SerialDenseMatrix b0new; b0new.Shape(b0.M(), 1);
        for (int i = 0; i < b0.M(); i++) {
            b0new(i, 1) = b0(i, k);
        }

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
            res1(x, 1) = sum - b0new(x, 1) - w(x, 1); // [...]-b0(:,k) - wsol;
            sum = 0;
        }
        // }

        // Compute number of contact nodes
        // @{
        int cont = 0;
        xvf.Shape(A.N(), k); yvf.Shape(A.N(), k); pf.Shape(A.N(), k);
        for (int i = 0; i < A.N(); i++) {
            if (A(i, 1) != 0) {
                cont += 1;
                xvf(cont, k) = xv0(i, k);
                yvf(cont, k) = yv0(i, k);
                pf(cont, k) = A(i, 1);
            }
        }
        nf(k, 1) = cont;
        // }

        // Compute contact force and contact area
        // @{
        force0.push_back(0);
        for (int i = 0; i < nf(k, 1); i++) {
            force0[k] = force0[k] + pf(i, k);
            area0.push_back(nf(k, 1) * (pow(delta, 2) / pow(lato, 2)) * 100);
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
    /*
    for (int i = 0; i < nf(k, 1); i++) {
        xvfaux(i, 0) = xvf(i, k);
        yvfaux(i, 0) = yvf(i, k);
        pfaux(i, 0) = pf(i, k);
    }
    */

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