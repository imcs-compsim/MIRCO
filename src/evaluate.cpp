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
#include <jsoncpp/json/json.h>
#include <memory>
using namespace std;

#include "topology.h"
#include "topologyfactory.h"
#include "linearsolver.h"
#include "nonlinearsolver.h"
#include "matrixsetup.h"
#include "warmstart.h"
#include "setparameters.h"
#include "evaluate.h"

#include "evaluate.h"
#include "writetofile.h"

void Evaluate(std::string jsonFileName, double &force)
{
    omp_set_num_threads(6); // 6 seems to be optimal

    auto start = std::chrono::high_resolution_clock::now();
    bool flagwarm;
    int n;
    double nu1, nu2, G1, G2, E, alpha, k_el, delta, nnodi, to1, E1,
        E2, lato, errf, sum = 0, Delta;
    bool rmg_flag;
    bool rand_seed_flag;
    double Hurst;
    string zfilePath;

    SetParameters(E1, E2, lato, nu1, nu2, G1, G2,
                  E, alpha, k_el, delta, nnodi, errf, to1, Delta, zfilePath, n, jsonFileName, rmg_flag, Hurst, rand_seed_flag, flagwarm);

    time_t now = time(0);
    tm *ltm = localtime(&now);
    std::cout << "Time is: " << ltm->tm_hour << ":";
    std::cout << 1 + ltm->tm_min << endl;

    // Meshgrid-Command
    // Identical Vectors/Matricies, therefore only created one here.
    int iter = int(ceil((lato - (delta / 2)) / delta)); // Replacement for "for (double i = delta / 2; i < lato; i = i + delta)"
    vector<double> x(iter);

#pragma omp parallel for schedule(static, 16) // Same amount of work -> static
    for (int i = 0; i < iter; i++)
    {
        x[i] = (delta / 2) + i * delta;
    }

    // Setup Topology
    Epetra_SerialDenseMatrix topology, y;
    int N = pow(2, n);
    topology.Shape(N + 1, N + 1);

    std::shared_ptr<TopologyGeneration> surfacegenerator;

    CreateSurfaceObject(n, Hurst, rand_seed_flag, zfilePath, rmg_flag, surfacegenerator); // creating the correct surface object

    surfacegenerator->GetSurface(topology);

    double zmax = 0;
    double zmean = 0;
    int cont = 0;

#pragma omp parallel for schedule(guided, 16) reduction(+                      \
                                                        : zmean) reduction(max \
                                                                           : zmax)
    // Static and Guided seem even but Guided makes more sense
    for (int i = 0; i < topology.M(); i++)
    {
        for (int j = 0; j < topology.N(); j++)
        {
            zmean += topology(i, j);
            if (topology(i, j) > zmax)
            {
                zmax = topology(i, j);
            }
        }
    }

    zmean = zmean / (topology.N() * topology.M());

    vector<double> area0;
    vector<double> force0;
    double w_el = 0, area = 0;
    int k = 0, n0 = 0;
    std::vector<double> xv0, yv0, b0, x0, xvf,
        yvf, pf;
    vector<int> col, row;
    double nf=0;
    Epetra_SerialDenseMatrix A;
    int nf2 = floor(nf);

    while (errf > to1 && k < 100)
    {
        // First predictor for contact set
        // All points, for which gap is bigger than the displacement of the rigid
        // indenter, cannot be in contact and thus are not checked in nonlinear solve
        // @{
            
        // [ind1,ind2]=find(z>=(zmax-(Delta(s)+w_el(k))));
        double value = zmax - Delta - w_el;
        row.clear();
        col.clear();

        // Data is even, guided makes more sense
        for (int i = 0; i < topology.N(); i++)
        {
            for (int j = 0; j < topology.N(); j++)
            {
                if (topology(i, j) >= value)
                {
                    row.push_back(i);
                    col.push_back(j);
                }
            }
        }

        n0 = col.size();

        // @{
        xv0.clear();
        xv0.resize(n0);
        yv0.clear();
        yv0.resize(n0);
        b0.clear();
        b0.resize(n0);
        // @} Parallelizing slows down program here, so not parallel

#pragma omp for schedule(guided, 16) // Always same workload but testing might be good -> Guided?
        for (int b = 0; b < n0; b++)
        {
            try
            {
                xv0[b] = x[col[b]];
            }
            catch (const std::exception &e)
            {
            }
        }

#pragma omp parallel for schedule(guided, 16) // Same
        for (int b = 0; b < n0; b++)
        {
            try
            {
                yv0[b] = x[row[b]];
            }
            catch (const std::exception &e)
            {
            }
        }

#pragma omp parallel for schedule(guided, 16) // Same
        for (int b = 0; b < n0; b++)
        {
            try
            {
                b0[b] = Delta + w_el - (zmax - topology(row[b], col[b]));
            }
            catch (const std::exception &e)
            {
            }
        }

        A.Shape(xv0.size(), xv0.size());

        // Construction of the Matrix H = A
        MatrixGeneration matrix1;
        matrix1.SetUpMatrix(A, xv0, yv0, delta, E, n0);

        // Second predictor for contact set
        // @{

        Epetra_SerialDenseMatrix xv0t, yv0t, xvft, yvft, pft, x0temp; // Temporary variables for warmup
        if (flagwarm == 1 && k > 1)
        {
            // x0=warm_x(xv0(1:n0(k),k),yv0(1:n0(k),k),xvf(1:nf(k-1),k-1),yvf(1:nf(k-1),k-1),pf(1:nf(k-1),k-1));
            xv0t.Shape(1, n0);
            yv0t.Shape(1, n0);
            xvft.Shape(1, nf2);
            yvft.Shape(1, nf2);
            pft.Shape(1, nf2);

#pragma omp parallel for schedule(static, 16) // Always same workload -> Static
            for (int i = 0; i < n0; i++)
            {
                xv0t(0, i) = xv0[i];
                yv0t(0, i) = yv0[i];
            }

#pragma omp parallel for schedule(static, 16) // Always same workload -> Static
            for (int j = 0; j < nf2; j++)
            {
                xvft(0, j) = xvf[j];
                yvft(0, j) = yvf[j];
                pft(0, j) = pf[j];
            }

            Warmstarter warm1;
            x0temp = warm1.Warmstart2(xv0t, yv0t, xvft, yvft, pft);

#pragma omp parallel for schedule(static, 16) // Always same workload -> Static
            for (int i = 0; i < x0temp.N(); i++)
            {
                x0[i] = x0temp(i, 1);
            }
        }
        else
        {
            if (b0.size() > 0)
            {
                // x0.Shape(b0.N(), 1);
                x0.resize(b0.size());
            }
        }
        // }

        // {
        x0.clear();
        x0.resize(b0.size());
        Epetra_SerialDenseMatrix b0new;
        b0new.Shape(b0.size(), 1);
        // } Parallel region makes around this makes program slower

#pragma omp parallel for schedule(static, 16) // Always same workload -> Static!
        for (long unsigned int i = 0; i < b0.size(); i++)
        {
            b0new(i, 0) = b0[i];
        }

        Epetra_SerialDenseMatrix w;

        NonLinearSolver solution2;
        solution2.NonlinearSolve(A, b0new, x0, w, y); // y -> sol, w -> wsol; x0 -> y0

        // Compute residial
        // @{
        Epetra_SerialDenseMatrix res1;
        if (A.M() != y.N())
        {
            std::runtime_error("Error 1: Matrix dimensions imcompatible");
        }
        res1.Shape(A.M(), A.M());

        // res1=A*sol-b0(:,k)-wsol;
#pragma omp parallel for schedule(static, 16) // Always same workload -> Static
        for (int x = 0; x < A.N(); x++)
        {
            for (int z = 0; z < y.M(); z++)
            {
                res1(x, 0) += A(x, z) * y(z, 0);
            }
            res1(x, 0) -= b0new(x, 0) + w(x, 0); // [...]-b0(:,k) - wsol;
        }
        // }

        // Compute number of contact node
        // @{

        // @{
        xvf.clear();
        xvf.resize(y.M());
        yvf.clear();
        yvf.resize(y.M());
        pf.resize(y.M());
        cont = 0;
        // @} Parallelizing this slows down program, so removed it.

#pragma omp for schedule(guided, 16)
        for (int i = 0; i < y.M(); i++)
        {
            if (y(i, 0) != 0)
            {
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

        sum = 0;
        iter = ceil(nf);
#pragma omp parallel for schedule(static, 16) reduction(+ \
                                                        : sum) // Always same workload -> Static!
        for (int i = 0; i < iter; i++)
        {
            sum += pf[i];
        }
        force0[k] += sum;
        area0.push_back(nf * (pow(delta, 2) / pow(lato, 2)) * 100);
        w_el = force0[k] / k_el;

        // }

        // Compute error due to nonlinear correction
        // @{
        if (k > 0)
        {
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
        cout << "Differenz ist zu groß!" << std::endl; // for nn=5

    auto finish = std::chrono::high_resolution_clock::now();
    double elapsedTime2 = std::chrono::duration_cast<std::chrono::seconds>(finish - start).count();
    std::cout << "Elapsed time is: " + to_string(elapsedTime2) + "s." << endl;

    writeForceToFile(y, zfilePath);
}