#include <cmath>
#include <ctime>
#include <fstream>
#include <random>

#include "mirco_topology_kokkos.h"

namespace MIRCO
{
  ViewMatrix_h CreateSurfaceFromFile(const std::string& filepath, int& dimension)
  {
    dimension = 0;

    std::ifstream reader(filepath);
    std::string line;

    while (getline(reader, line))
    {
      ++dimension;
    }
    reader.clear();
    // Reuse reader
    reader.seekg(0, std::ios::beg);

    ViewMatrix_h z("CreateSurfaceFromFile(); z", dimension, dimension);
    std::ifstream stream(filepath);
    int lineCounter = 0;
    while (getline(reader, line))
    {
      ++lineCounter;
      for (int i = 0; i < dimension; i++)
      {
        const int separatorPosition = line.find_first_of(';');
        z(lineCounter - 1, i) = stod(line.substr(0, separatorPosition));
        line = line.substr(separatorPosition + 1, line.length());
      }
    }
    reader.close();

    return z;
  }

  ViewMatrix_h CreateRmgSurface(int resolution, double InitialTopologyStdDeviation, double Hurst,
      bool RandomSeedFlag, int RandomGeneratorSeed)
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

    int N = (1 << resolution) + 1;
    ViewMatrix_h z("CreateRmgSurface(); z", N, N);

    const double scaling_factor = pow(2.0, 0.5 * Hurst);
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

}  // namespace MIRCO
