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
            double& H, double& rnd, double& k_el, double& delta, double& nnodi) {
        E1 = 1; E2 = 1;
        nu1 = 0.3; nu2 = 0.3;
        G1 = E1 / (2 * (1 + nu1));
        G2 = E2 / (2 * (1 + nu2));
        E = pow(((1 - pow(nu1, 2)) / E1 + 1 - pow(nu2, 2) / E2), -1);
        G = pow(((2 - nu1) / (4 * G1) + (2 - nu2) / (4 * G2)), -1);
        nu = E / (2 * G) - 1;

        vector<double> alpha_con{ 0.778958541513360, 0.805513388666376, 0.826126871395416, 0.841369158110513,
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
    while (getline(reader, blaLine)){dimension += 1;}
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
    catch (const std::exception & e) {
        //Fatal Error, catch it!
        std::cout << e.what(); // Opens a standart error message box, at least StackOverflow said that.
    }
}

void createTopology(int systemsize, Epetra_SerialDenseMatrix& topology) {
    // Idea: GENERATING a topology, instead of I/O-Topology!
    // TODO: Insert ranmid2d_MP-Method here!

     /* Only needed once ranmid2d_NP is implemented to create new topology? Put in CreateTopology
        once we are there.

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
    */
}

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

    for (int i = 0; i < systemsize; i++){
        for (int j = 0; j < i; j++) {
                r = (xv0(j, k) - xv0(i, k)) * (yv0(j ,k) - yv0(i, k));
                A(i, j) = C * asin(raggio / r);
        }
    }
    return A;
}

/*------------------------------------------*/
// pf, xvf, yvf -> wird dann in der Main implementiert werden müssen! (TODO)
Epetra_SerialDenseMatrix Warmstart(Epetra_SerialDenseMatrix xv0, Epetra_SerialDenseMatrix yv0, Epetra_SerialDenseMatrix &xvf,
        Epetra_SerialDenseMatrix &yvf, Epetra_SerialDenseMatrix& pf) {
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
    vector<int> index;
    int counter = 0;
    for (int i = 0; i < pf.N(); i++) {
        // ind=find(matfin(:,1)==xvf(i) & matfin(:,2)==yvf(i));
        for (int j = 0; j < xvf.N(); j++) {
            if ((combinedMatrix(j, 1) == xvf(i, 1)) && (combinedMatrix(j, 2) == yvf(i, 1))) {
                index[counter] = j;
            }
        }

        // x0(ind,1)=pf(i); TODO
        for (int y = 0; y < counter; y++) {
            x0(y, 1) = pf(y, 1);
        }
        index.clear();
        counter = 0;
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
    if (err != 0 ) { std::cout << "Error setting up maxtix solver (2)"; }

    err = solver.Solve();
    if (err != 0) { std::cout << "Error setting up matrix solver (3)"; }
    std::cout << vector_x << std::endl;
}

/*------------------------------------------*/

void NonlinearSolve(Epetra_SerialSymDenseMatrix& matrix, Epetra_SerialDenseMatrix& b0, 
                    Epetra_SerialDenseMatrix& x0, Epetra_SerialDenseMatrix& w, int iter, Epetra_SerialDenseMatrix& y) {
    Epetra_SerialDenseMatrix vector_x;
    Epetra_SerialDenseMatrix vector_b;
    Epetra_SerialSymDenseMatrix local_matrix;
    double nnlstol = 1.0000e-08;
    double maxiter = 10000;
    bool init = false;
    int n0 = b0.N() * b0.M();
    y.Shape(n0, 1);
    vector<int> P(n0), yf(n0);
    iter = 0;
    // Initialize active set, y0 -> changed to yf (done, namechange)

    // !!! This seems absolute garbage: !!!
    // inz=find(y0>=nnlstol);
    // nP = numel(inz);
    int counter = 0;
    vector<int> positions;
    for (int i = 0; i < yf.size(); i++) {
        if ((yf[i] == nnlstol) || (yf[i] > nnlstol)) {
            positions[counter] = i;
            counter += 1;
        }
    }

    if (counter == 0) {
        w.Shape(b0.N(), b0.M());
        for (int x = 0; x < b0.N(); x++) {
            for (int y = 0; y < b0.M(); y++) {
                w(x, y) = -b0(x, y);
            }
        }
    } else {
        for (int i = 0; i < counter; i++) {
            P[i] = positions[i];
        }
        w.Shape(counter, 1);
        init = true;
    }
    
    vector<double> s(n0);
    bool aux1 = true, aux2 = true;
    while (aux1 == true) {
        // [wi,i]=min(w);
        int minValue = w(0, 0), minPosition = 0;
        for (int i = 0; i < w.M(); i++){
            if (minValue > w(0, i)) {
                minValue = w(0, i);
                minPosition = i;
            }
        }

        if ((counter == n0) || ( minValue > -nnlstol) || ((iter == maxiter) || (iter > maxiter)) && (init == false)) {
            aux1 = false;
        } else {
            if (init == false) {
                // Index #i enter active index
                counter += 1;
                P[counter] = minPosition;
            }
        }
        int j = 0;
        double eps = 2.2204e-16, alphai = 0, alpha = 100000000, a = 0;
        while (aux2 == true) {
            iter += 1;

            // Bringe die Matrizen in die richtige Vektorform
             vector_x.Shape(counter, 1);
             vector_b.Shape(counter, 1);
             local_matrix.Shape(counter, counter);

             for (int x = 0; x < counter; x++) {
               vector_b(x, 1) = b0(P[x], 1);
               for (int x = 0; x < counter; x++) {
                 local_matrix(x, x) = matrix(P[x], P[x]);
               }
             }

             // Rufe den linearen Lösungsalgorithmus
             LinearSolve(matrix, vector_x, vector_b);

             for (int x = 0; x < counter; x++) {
               s[P[x]] = vector_b(x, 1);
             }

            bool allBigger = true;
            for (int x = 0; x < counter; x++) {
                if (s[P[x]] < nnlstol) { allBigger = false; }
            }
            if (allBigger == true) {
                aux2 = false;
                // Linear solve für w=A(:,P(1:nP))*y(P(1:nP))-b; ? TODO!!



                return; // TODO: Return in einer Void-Funktion?
            } else {
                for (int i = 0; i < counter; i++) {
                    if (s[P[i]] < nnlstol) {
                        alphai = y(P[i], 1) / (eps + y(P[i], 1) - s[P[i]]);
                        if (alphai < alpha) {
                            alpha = alphai;
                            j = 1;
                        }
                    }
                }
            }

            while (a < counter) {
                a += 1;
                y(P[a], 1) = y(P[a], 1) + alpha * (s[P[a]] - y(P[a], 1));
            }

            if (j > 0) {
                // jth entry in P leaves active set
                s[P[j]] = 0;
                vector<int> P2 = P;
                for (int i = j; i < (counter - 1); i++) {
                    P2[i + 1] = P[i];
                }
                for (int i = 0; i < j; i++) {
                    P2[i] = P[i];
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
    int E1, E2, csteps, flagwarm, lato, zref, ampface;
    double nu1, nu2, G1, G2, E, G, nu, alpha, H, rnd, k_el, delta, nnodi;
    SetParameters(E1, E2, csteps, flagwarm, lato, zref, ampface, nu1, nu2, G1, G2, E, G, nu, alpha, H, rnd, k_el, delta, nnodi);

    // Meshgrid-Command
    // Identical Vectors/Matricies, therefore only created one here.
    vector<double> x;
    int iterator = 0;
    for (int i = delta / 2; i < (lato - delta / 2); i = i + delta) {
        x[iterator] = i;
        iterator += 1;
    }
    string randomPath = "C:\..."; // TODO: Change this before debugging!
    Epetra_SerialSymDenseMatrix topology, y;
    CreateTopology(topology.N(), topology, randomPath); 
    // TODO: Remove 3rd argument when ranmid2d_MP is implemented!
    // hard-code file path for debugging? Nicer, hand into main as command line argument

    double zmax = 0;
    double zmean = 0;
    zmean = topology.NormOne() / pow(topology.N(), 2);

    // Can also use zmatrix = topology.NormInf() and
    // zmean = topology.NormOne()/pow(topology.N(), 2)
    // topology.N() should send dimension of matrix
    for (int i = 1; i < topology.N() + 1; i++) {
        for (int j = 1; j < topology.N() + 1; j++) {
            if (zmax < topology(i, j)) {
                zmax = topology(i, j);
            }
        }
    }

    vector<double> nfaux(csteps), Delta(csteps), force(csteps), area(csteps);

    // For (maybe) future loop: Here!
    // for (int s = 0; s < csteps; s++){
    // Delta(0) = ampface * (0.5) * zmean * 1 / csteps;
    // Delta(s) = ampface * (0.5) * zmean * s / csteps;
    Delta[0] = 39.202067343593399;
    int errf = 100000000;
    float to1 = 0.01;
    vector<double> force0(csteps), area0(csteps), w_el0(csteps);
    vector<int> w_el; // Change w_el to vector for loop
    w_el[0] = 0;
    int k = 0;
    vector<int> n0;
    Epetra_SerialDenseMatrix xv0, yv0, b0, x0, nf, xvfaux, yvfaux, pfaux, xvf, yvf, pf;
    nf.Shape(csteps, 1); xvfaux.Shape(csteps, 1); yvfaux.Shape(csteps, 1); pfaux.Shape(csteps, 1);
    while (errf > to1) {
        k += 1;

        // [ind1,ind2]=find(z>=(zmax-(Delta(s)+w_el(k))));
        vector<int> col, row;
        int counter = 0;
        double value = zmax - Delta[0] + w_el[k];
        for (int i = 0; i < topology.N(); i++) {
            for (int j = 0; j < topology.N(); j++) {
                if ((topology(i, j) > value) || (topology(i, j) == value)) {
                    col[counter] = i;
                    row[counter] = j;
                    counter += 1;
                }
            }
        }
        n0[k] = col.size();

        for (int i = 0; i < n0[k]; i++) {
            xv0(i, k) = x[row[i]];
            yv0(i, k) = y(row[i], 1);
            b0(i, k) = Delta[0] + w_el[k] - (zmax - topology(row[i], col[i]));
        }

        // Construction of the Matrix H = A
        Epetra_SerialSymDenseMatrix A = SetUpMatrix(xv0, yv0, delta, E, n0[k], k);

        // TODO: Fertig machen! // TODO: xvf, yvf, pf werden vor der Initialisierung verwendet!
        Epetra_SerialDenseMatrix xv0t, yv0t, xvft, yvft, pft, xvfauxt, yvfauxt, pfauxt; // Temporary variables for warmup
        if (flagwarm == 1 && k > 1) {
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
        } 
        /*
             else if (flagwarm == 1 && k == 1 && s > 1){
            // x0=warm_x(xv0(1:n0(k),k),yv0(1:n0(k),k),xvfaux(1:nfaux(s-1),s-1),yvfaux(1:nfaux(s-1),s-1),pfaux(1:nfaux(s-1),s-1));
            xv0t.Shape(1, n0[k]); yv0t.Shape(1, n0[k]); xvfauxt.Shape(1, nfaux[0]);
            yvfauxt.Shape(1, nfaux[0]); pfauxt.Shape(1, nfaux[0]); // Change this for loop!
            for (int i = 0; i < n0[k]; i++){
                xv0t(1, i) = xv0(i, k);
                yv0t(1, i) = yv0(i, k);
            }
            for (int i = 0; i < 1; i++) { // Change for loop!
                xvfauxt(1, i) = nfaux[i];
                yvfauxt(1, i) = nfaux[i];
                pfauxt(1, i) = nfaux[i];
            }
        }
        */
        else {
            x0.Shape(b0.N(), 1);
        }

        Epetra_SerialDenseMatrix b0new; b0new.Shape(b0.M(), 1);
        for (int i = 0; i < b0.M(); i++) {
            b0new(i, 1) = b0(i, k);
        }
        Epetra_SerialDenseMatrix w;
        int iter;
        NonlinearSolve(A, b0new, x0, w, iter, y); // y -> sol, w -> wsol
        
        Epetra_SerialDenseMatrix res1;
        if (A.M() != y.N()) { std::runtime_error("Error 1: Matrix dimensions imcompatible"); }
        res1.Shape(A.N(), y.M());
        // res1=A*sol-b0(:,k)-wsol;
        // For some weird reason, adding a vector to a matrix adds it to every coloum/row (???)
        int sum = 0;
        for (int x = 0; x < A.N(); x++) {
            for (int z = 0; z < y.M(); z++) {
                for (int i = 0; i < A.M(); i++) {
                    sum += A(x, i) * y(i, z);
                }
                res1(x, z) = sum - b0new(z, 1) - w(x, z); // [...]-b0(:,k) - wsol;
                sum = 0;
            }
        }

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
        force0[k] = 0;
        for (int i = 0; i < nf(k, 1); i++) {
            force0[k] = force0[k] + pf(i, k);
            area0[k] = nf(k, 1) * (pow(delta, 2) / pow(lato, 2)) * 100;
        }
        w_el0[k + 1] = force0[k] / k_el;
        if (k > 1) {
            // errf(k) = abs((force0(k)-force0(k-1))/force0(k));
            errf = (force0[k] - force0[k - 1]) / force0[k];
            if (errf < 0) { errf = (-1) * errf; }

            // errw(k) = abs((w_el0(k+1)-w_el0(k))/w_el0(k+1));
            // It appears that this is only a debugging variable without any uses, therefore im not gonna implement this here.
        }


    }

    // Change this for global loop!
    for (int i = 0; i < nf(k, 1); i++) {
        xvfaux(i, 0) = xvf(i, k);
        yvfaux(i, 0) = yvf(i, k);
        pfaux(i, 0) = pf(i, k);
    }
    // Change this for global loop aswell!
    force[0] = force0[k];
    area[0] = area0[k];
    w_el[0] = w_el0[k];
    // End of loop:
    // }

    // Mean pressure
    double sigmaz = force[0] / pow(lato, 2);
    // Pressure unit per depth
    double pressz = sigmaz; // ????

    // End of Code
}
