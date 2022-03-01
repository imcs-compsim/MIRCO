#include <Epetra_SerialDenseMatrix.h>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
using namespace std;

#include "topology.h"

void ReadFile::GetSurface(Epetra_SerialDenseMatrix &z)
{
  ifstream reader(filepath);
  string blaLine;
  int dimension = 0;
  while (getline(reader, blaLine))
  {
    dimension += 1;
  }
  reader.close();
  z.Shape(dimension, dimension);
  int position = 0, separatorPosition, lineCounter = 0;
  ifstream stream(filepath);
  string line, container;
  double value;
  while (getline(stream, line))
  {
    separatorPosition = 0;
    lineCounter += 1;
    for (int i = 0; i < dimension; i++)
    {
      separatorPosition = line.find_first_of(';');
      container = line.substr(0, separatorPosition);
      position += 1;
      line = line.substr(separatorPosition + 1, line.length());
      value = stod(container);
      z(lineCounter - 1, i) = value;
    }
  }
  stream.close();
}

void Rmg::GetSurface(Epetra_SerialDenseMatrix &z)
{
  srand(time(NULL));

  int seed;

  if (rand_seed_flag)
  {
    seed = rand();
  }
  else
  {
    seed = rmg_seed;  // seed can be fixed to reproduce result
  }
  std::default_random_engine generate(seed);
  std::normal_distribution<double> distribution(
      0.0, 1.0);  // normal distribution: mean = 0.0, standard deviation = 1.0

  int N = pow(2, resolution);
  double alpha = 1 / sqrt(0.09);

  int D = N;
  int d = N / 2;

  for (int i = 0; i < resolution; i++)
  {
    alpha = alpha / (pow(sqrt(2), Hurst));

    for (int j = d; j < N - d + 1; j = j + D)
    {
      for (int k = d; k < N - d + 1; k = k + D)
      {
        z(j, k) = (z(j + d, k + d) + z(j + d, k - d) + z(j - d, k + d) + z(j - d, k - d)) / 4 +
                  alpha * distribution(generate);
      }
    }

    alpha = alpha / (pow(sqrt(2), Hurst));

    for (int j = d; j < N - d + 1; j = j + D)
    {
      z(j, 0) = (z(j + d, 0) + z(j - d, 0) + z(j, d)) / 3 + alpha * distribution(generate);
      z(j, N) = (z(j + d, N) + z(j - d, N) + z(j, N - d)) / 3 + alpha * distribution(generate);
      z(0, j) = (z(0, j + d) + z(0, j - d) + z(d, j)) / 3 + alpha * distribution(generate);
      z(N, j) = (z(N, j + d) + z(N, j - d) + z(N - d, j)) / 3 + alpha * distribution(generate);
    }

    for (int j = d; j < N - d + 1; j = j + D)
    {
      for (int k = D; k < N - d + 1; k = k + D)
      {
        z(j, k) = (z(j, k + d) + z(j, k - d) + z(j + d, k) + z(j - d, k)) / 4 +
                  alpha * distribution(generate);
      }
    }

    for (int j = D; j < N - d + 1; j = j + D)
    {
      for (int k = d; k < N - d + 1; k = k + D)
      {
        z(j, k) = (z(j, k + d) + z(j, k - d) + z(j + d, k) + z(j - d, k)) / 4 +
                  alpha * distribution(generate);
      }
    }

    D = D / 2;
    d = d / 2;
  }

  double zref = 50;
  double zmax = z(0, 0);

  double zmean;
  double sum = 0;

  for (int i = 0; i < N + 1; i++)
  {
    for (int j = 0; j < N + 1; j++)
    {
      if (zmax < z(i, j))
      {
        zmax = z(i, j);
      }
      sum = sum + z(i, j);
    }
  }
  zmean = sum / (pow((N + 1), 2));
  double scalefactor = zref / (zmax - zmean);
  for (int i = 0; i < N + 1; i++)
  {
    for (int j = 0; j < N + 1; j++)
    {
      z(i, j) = scalefactor * z(i, j);
    }
  }
  double zmin = z(0, 0);
  for (int i = 0; i < N + 1; i++)
  {
    for (int j = 0; j < N + 1; j++)
    {
      if (zmin > z(i, j))
      {
        zmin = z(i, j);
      }
    }
  }
  for (int i = 0; i < N + 1; i++)
  {
    for (int j = 0; j < N + 1; j++)
    {
      z(i, j) = z(i, j) - zmin;
    }
  }
}
