#include <cmath>
#include <cstdio>
#include <fstream> //ifstream
#include <string> //std::to_string, std::stod
#include "include/Epetra_SerialSpdDenseSolver.h"
#include "include/Epetra_SerialSymDenseMatrix.h"

using namespace std


/**
* Outdated function, implemented into Main-method for code structure.
*/
void SetParameters() {}

/*------------------------------------------*/

/**
* Readin Matrix from a file. Elements have to be separated by a ';'.
* TODO: Change this to spaces.
*/
Epetra_SerialSymDenseMatrix& CreateTopology(int systemsize, Epetra_SerialSymDenseMatrix& topology, string filePath) {
    // Readin for amount of lines -> dimension of matrix
    ifstream reader(filePath);
    string blaLine;
    int dimension = 0;
    while (getLine(reader, blaLine){dimension += 1;}
    reader.close();
    topology.shape(lineAmount);
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

    // Place elements into the matrix -> Should happen in the loop already!
    // int dimension = sqrt(position); // square matrix, allowed then.
    // topology.shape(dimension);
    // for (int i = 0; i < dimension; i++) {
    //    for (int j = 0; j < dimension; j++) {
    //        topology(i, j) = elements[i + j]; // Enters value for line i, colloumn j - just as read into the elements-array.
    //    }
    //}
    return topology;
  // 3. wollen wir vielleich die Funktion aus ranmid2d_MP nachimplementieren (Später)
}

// Used to be in "CreateTopology", now moved here for code structure.
/**
* Outdated version to create a Matrix. Use @CreateTopology to readin matrix.
*/
Epetra_SerialSymDenseMatrix& createSimpleMatrix(int systemSize){
    // Shapes a matrix to be a square Matrix. Main diagonal will be 2's, rest is 1's.
    topology.Shape(systemSize);
    for (int i = 0; i < systemSize; i++) {
        topology(i, i) = 2;
        for (int j = 0; j < i; j++) {
            topology(i, j) = 1;
            topology(j, i) = 1;
        }
    }
}

/*------------------------------------------*/

void SetUpMatrix() {
  // Hier wollen wir die Konstitutivematrix A aufstellen
}

/*------------------------------------------*/

void Warmstart() {
  // Das Programm ist effizienter, wenn wir einen warmstart machen. Das Ergebnis
  // sollte sich aber nicht ändern. Diese Funktion kann also auch erst ganz zum
  // Schluss implementiert werden.
}

/*------------------------------------------*/
// Die LinearSolve Funktion kann als blackbox betrachtet werden, in die die
// Matrix A, der Vektor x und der Vektor b der Gleichung Ax=b gegeben werden
void LinearSolve(Epetra_SerialSymDenseMatrix& matrix,
                 Epetra_SerialDenseMatrix& vector_x,
                 Epetra_SerialDenseMatrix& vector_b) {
  // Erstelle ein Objekt des Lösungsalgorithmus für Ax=b
  Epetra_SerialSpdDenseSolver solver;

  // Gebe dem Lösungsalgorithmus die zu lösende Matrix
  int err = solver.SetMatrix(matrix);

  // Gebe dem Lösungsalgorithmus den x vektor, in den die lösung reingeschrieben
  // werden soll und die rechte Seite b
  err = solver.SetVectors(vector_x, vector_b);

  // Löse das gegebene Gleichungssystem
  err = solver.Solve();

  // Gebe die Lösung in das Terminal aus (zum debuggen)
  std::cout << vector_x << std::endl;
}

/*------------------------------------------*/

void NonlinearSolve(int systemsize, Epetra_SerialSymDenseMatrix& matrix) {
  // Erstelle 2 Matrix objekte
  Epetra_SerialDenseMatrix vector_x;
  Epetra_SerialDenseMatrix vector_b;

  // todo das ist nicht die wirkliche Matrixgröße, wie sie in MATLAB
  // implementiert ist!

  // Bringe die matrizen in die richtige Vektorform
  vector_x.Shape(systemsize, 1);
  vector_b.Shape(systemsize, 1);

  // Befülle die rechte Seite b
  for (int i = 0; i < systemsize; i++) {
    vector_b(i, 0) = 1;
  }

  // Rufe den linearen Lösungsalgorithmus
  LinearSolve(matrix, vector_x, vector_b);

  // Hier soll die Funktion nnls rein
}
/*------------------------------------------*/

int main(int argc, char* argv[]) {
    int E1 = 1, E2 = 1;
    double nu1 = 0.3, nu2 = 0.3;
    double G1 = E1 / (2 * (1 + nu1)), G2 = E2 / (2 * (1 + nu2));
    double E = pow(((1 - nu1 ^ 2) / E1 + (1 - pow(nu2, 2)) / E2), -1);
    double G = pow(((2 - nu1) / (4 * G1) + (2 - nu2) / (4 * G2)), -1);
    double nu = E / (2 * G) - 1;

    vector<double> alpha_con = { 0.778958541513360, 0.805513388666376, 0.826126871395416, 0.841369158110513,
    0.851733020725652, 0.858342234203154, 0.862368243479785, 0.864741597831785 };
    int nn = 2; // Matrix sent has the parameter nn=2!
    double alpha = alpha_con[nn];
    int csteps = 50;
    double ampface = 1;
    int flagwarm = 1;
    int lato = 1000; // Lateral side of the surface [micrometers]
    double H = 0.1; // Hurst Exponent (D = 3 - H)
    double rnd = 95.0129;

    // vector<double> z = radmid2d_MP(nn, H, 1, 1, rnd); -> Noch zu definieren!

    int zref = 50; // Reference for the Scaling, former value = 25
    // double scalefactor = zref(max(max(z)), - mean(mean(z))) -> mean, max, min in anderem Projekt erstellt, dort überprüfen und anschließend hier einfügen!
    // z = scalefactor * z; -> Skalare Mutiplikation! Durch eine Loop laufen lassen, um dies in C++ zu schaffen!
    // WICHTIGER HINWEIS: z ist eine Matrix, kein Vektor!

    // int systemsize = pow(2, nn); //helper for dimension, outdated

    // TODO: Matrix, welche ich bekommen habe sollte nn=2 sein?

    // Definiere die Parameter
    double k_el = lato * E / alpha;
    double delta = lato / pow(2, nn + 1);
    double nnodi = pow(pow(2, nn + 1), 2);
    vector<vector<double>> x; //All Rows are identical
    vector<vector<double>> y; //All Colloums are identical
    // Access: x [Row][Colloum]!


    // Translation in Progress, Vector sufficient?
    //[x,y]=meshgrid(delta/2:delta:(lato-delta/2),delta/2:delta:(lato-delta/2));

    // Rufe die anderen Funktionen in der richtigen Reihenfolge

    Epetra_SerialSymDenseMatrix topology;

    topology = CreateTopology(systemsize, topology);

    NonlinearSolve(systemsize, topology);
}
