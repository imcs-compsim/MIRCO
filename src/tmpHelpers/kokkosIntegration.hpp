#pragma once

/*
#include <chrono>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>
*/

#include <Kokkos_Core.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>

// These can be adjusted later to use CudaSpace, etc.
using ExecSpace = Kokkos::DefaultExecutionSpace;
using HostSpace = Kokkos::HostSpace;

using ViewVector = Kokkos::View<double*, HostSpace>;
using ViewMatrix =
    Kokkos::View<double**, Kokkos::LayoutRight, HostSpace>;  // LayoutRight = row-major

inline ViewVector teuchosVectorToKokkos(const Teuchos::SerialDenseVector<int, double>& vec)
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

inline ViewMatrix teuchosMatrixToKokkos(const Teuchos::SerialDenseMatrix<int, double>& mat)
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
