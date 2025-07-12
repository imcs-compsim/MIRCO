#pragma once

/*
#include <chrono>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>
*/

// for switching between kokkos and openmp during integration process
#define kokkosElseOpenMP true

#include <Kokkos_Core.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>

// These can be adjusted later to use CudaSpace, etc.
using ExecSpace = Kokkos::DefaultExecutionSpace;
using HostSpace = Kokkos::HostSpace;

using ViewVector = Kokkos::View<double*, HostSpace>;
using ViewMatrix =
    Kokkos::View<double**, Kokkos::LayoutRight, HostSpace>;  // LayoutRight = row-major

inline ViewVector toKokkos(const Teuchos::SerialDenseVector<int, double>& vec)
{
  int n = vec.length();
  ViewVector kokkosVec("kokkosVec", n);
  for (int i = 0; i < n; ++i)
  {
    kokkosVec(i) = vec[i];
  }
  return kokkosVec;
}

inline Teuchos::SerialDenseVector<int, double> kokkosVectorToTeuchos(const ViewVector& kokkosVec)
{
  int n = kokkosVec.extent(0);
  Teuchos::SerialDenseVector<int, double> vec(n);
  for (int i = 0; i < n; ++i)
  {
    vec[i] = kokkosVec(i);
  }
  return vec;
}

inline ViewMatrix toKokkos(const Teuchos::SerialDenseMatrix<int, double>& mat)
{
  int numRows = mat.numRows();
  int numCols = mat.numCols();
  ViewMatrix kokkosMat("kokkosMat", numRows, numCols);
  for (int i = 0; i < numRows; ++i)
    for (int j = 0; j < numCols; ++j)
      kokkosMat(i, j) = mat(i, j);  // Access via (row, col) despite column-major internal layout
  return kokkosMat;
}

inline Teuchos::SerialDenseMatrix<int, double> kokkosMatrixToTeuchos(const ViewMatrix& kokkosMat)
{
  int numRows = kokkosMat.extent(0);
  int numCols = kokkosMat.extent(1);
  Teuchos::SerialDenseMatrix<int, double> mat(numRows, numCols);
  for (int i = 0; i < numRows; ++i)
    for (int j = 0; j < numCols; ++j) mat(i, j) = kokkosMat(i, j);
  return mat;
}

inline ViewVector toKokkos(const std::vector<double>& stdVec)
{
  ViewVector view("b0_kokkos", stdVec.size());
  for (std::size_t i = 0; i < stdVec.size(); ++i) view(i) = stdVec[i];
  return view;
}

inline std::vector<double> kokkosVectorToStdVector(const ViewVector& kokkosVec)
{
  int n = kokkosVec.extent(0);
  std::vector<double> stdVec(n);
  for (int i = 0; i < n; ++i)
  {
    stdVec[i] = kokkosVec(i);
  }
  return stdVec;
}
