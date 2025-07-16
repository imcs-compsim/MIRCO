#include "mirco_topology.h"

// ##include <Teuchos_SerialDenseMatrix.hpp>
#include <cmath>
#include <ctime>
#include <fstream>
#include <random>
#include <string>

Teuchos::SerialDenseMatrix<int, double> MIRCO::CreateSurfaceFromFile(
    const std::string& filepath, int& N)
{
  std::ifstream reader(filepath);
  std::string blaLine;
  int dimension = 0;
  while (getline(reader, blaLine))
  {
    dimension += 1;
  }
  reader.close();
  Teuchos::SerialDenseMatrix<int, double> z(dimension, dimension);
  int separatorPosition, lineCounter = 0;
  std::ifstream stream(filepath);
  std::string line, container;
  double value;
  while (getline(stream, line))
  {
    separatorPosition = 0;
    lineCounter += 1;
    for (int i = 0; i < dimension; i++)
    {
      separatorPosition = line.find_first_of(';');
      container = line.substr(0, separatorPosition);
      line = line.substr(separatorPosition + 1, line.length());
      value = stod(container);
      z(lineCounter - 1, i) = value;
    }
  }
  stream.close();

  N = dimension;
  return z;
}

Teuchos::SerialDenseMatrix<int, double> MIRCO::CreateRmgSurface(int resolution,
    double InitialTopologyStdDeviation, double Hurst, bool RandomSeedFlag, int RandomGeneratorSeed)
{
  srand(time(NULL));

  int seed;

  if (RandomSeedFlag)
  {
    seed = rand();
  }
  else
  {
    seed = RandomGeneratorSeed;  // seed can be fixed to reproduce result
  }
  std::default_random_engine generate(seed);
  std::normal_distribution<double> distribution(
      0.0, 1.0);  // normal distribution: mean = 0.0, standard deviation = 1.0

  int N = pow(2, resolution) + 1;
  Teuchos::SerialDenseMatrix<int, double> z(N, N);

  const double scaling_factor = pow(sqrt(2), Hurst);
  double alpha = InitialTopologyStdDeviation * scaling_factor;

  const int D_0 = N - 1;
  int D = D_0;
  int d = D_0 / 2;

  for (int i = 0; i < resolution; i++)
  {
    alpha = alpha / scaling_factor;

    for (int j = d; j < D_0 - d + 1; j = j + D)
    {
      for (int k = d; k < D_0 - d + 1; k = k + D)
      {
        z(j, k) = (z(j + d, k + d) + z(j + d, k - d) + z(j - d, k + d) + z(j - d, k - d)) / 4 +
                  alpha * distribution(generate);
      }
    }

    alpha = alpha / scaling_factor;

    for (int j = d; j < D_0 - d + 1; j = j + D)
    {
      z(j, 0) = (z(j + d, 0) + z(j - d, 0) + z(j, d)) / 3 + alpha * distribution(generate);
      z(j, D_0) =
          (z(j + d, D_0) + z(j - d, D_0) + z(j, D_0 - d)) / 3 + alpha * distribution(generate);
      z(0, j) = (z(0, j + d) + z(0, j - d) + z(d, j)) / 3 + alpha * distribution(generate);
      z(D_0, j) =
          (z(D_0, j + d) + z(D_0, j - d) + z(D_0 - d, j)) / 3 + alpha * distribution(generate);
    }

    for (int j = d; j < D_0 - d + 1; j = j + D)
    {
      for (int k = D; k < D_0 - d + 1; k = k + D)
      {
        z(j, k) = (z(j, k + d) + z(j, k - d) + z(j + d, k) + z(j - d, k)) / 4 +
                  alpha * distribution(generate);
      }
    }

    for (int j = D; j < D_0 - d + 1; j = j + D)
    {
      for (int k = d; k < D_0 - d + 1; k = k + D)
      {
        z(j, k) = (z(j, k + d) + z(j, k - d) + z(j + d, k) + z(j - d, k)) / 4 +
                  alpha * distribution(generate);
      }
    }

    D = D / 2;
    d = d / 2;
  }

  // Finding minimum of topology
  double zmin = std::numeric_limits<double>::max();
  for (int i = 0; i < D_0 + 1; i++)
  {
    for (int j = 0; j < D_0 + 1; j++)
    {
      zmin = std::min(zmin, z(i, j));
    }
  }

  // Setting the minimum of topology to zero
  for (int i = 0; i < D_0 + 1; i++)
  {
    for (int j = 0; j < D_0 + 1; j++)
    {
      z(i, j) = z(i, j) - zmin;
    }
  }

  return z;
}
