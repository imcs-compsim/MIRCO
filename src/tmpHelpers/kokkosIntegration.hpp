#pragma once


// #include <chrono>
// #include <mutex>
#include <string>
// #include <unordered_map>
#include <vector>


// for switching between kokkos and openmp during integration process
#define kokkosElseOpenMP true

#include <Kokkos_Core.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>

// These can be adjusted later to use CudaSpace, etc.
using ExecSpace_DefaultHost_t = Kokkos::DefaultHostExecutionSpace;
// using ExecSpace_Cuda_t = Kokkos::Cuda;
using ExecSpace_Default_t = Kokkos::DefaultExecutionSpace;

using MemorySpace_Host_t = Kokkos::HostSpace;
// using MemorySpace_Cuda_t = Kokkos::CudaSpace;
using MemorySpace_ofDefaultExec_t = ExecSpace_Default_t::memory_space;

using Device_Host_t = Kokkos::Device<ExecSpace_DefaultHost_t, MemorySpace_Host_t>;
/*using Device_Cuda_t = Kokkos::Device<
  ExecSpace_Default_t,
  MemorySpace_Host_t
  >;*/
using Device_Default_t = Kokkos::Device<ExecSpace_Default_t, MemorySpace_ofDefaultExec_t>;

// # or LayoutRight for some things? no idea
// # LayoutRight = row-major, which is what Teuchos uses
using ViewVector_h = Kokkos::View<double*, Device_Host_t>;
using ViewMatrix_h =
    Kokkos::View<double**, Kokkos::LayoutLeft, Device_Host_t>;  // LayoutRight = row-major

using ViewVectorInt_h = Kokkos::View<int*, Device_Host_t>;    // # or "array"?
using ViewVectorBool_h = Kokkos::View<bool*, Device_Host_t>;  // #

using ViewScalarInt_h = Kokkos::View<int, Device_Host_t>;
using ViewScalarDouble_h = Kokkos::View<double, Device_Host_t>;

/*using ViewVector_cuda = Kokkos::View<double*, Device_Cuda_t>;
using ViewMatrix_cuda =
    Kokkos::View<double**, Kokkos::LayoutLeft, Device_Cuda_t>;*/  // LayoutRight = row-major

using ViewVector_d = Kokkos::View<double*, Device_Default_t>;
using ViewMatrix_d =
    Kokkos::View<double**, Kokkos::LayoutLeft, Device_Default_t>;  // LayoutRight = row-major

using ViewVectorInt_d = Kokkos::View<int*, Device_Default_t>;
using ViewVectorBool_d = Kokkos::View<bool*, Device_Default_t>;

using ViewScalarInt_d = Kokkos::View<int, Device_Default_t>;
using ViewScalarDouble_d = Kokkos::View<double, Device_Default_t>;


#if (kokkosElseOpenMP)
constexpr bool HOSTONLY = std::is_same<MemorySpace_Host_t, MemorySpace_ofDefaultExec_t>();
#endif


inline ViewVector_d toKokkos(const Teuchos::SerialDenseVector<int, double>& vec)
{
  int n = vec.length();
  ViewVector_d kokkosVec("kokkosVec", n);

  // Always create host mirror (may just return kokkosVec itself if already in HostSpace)
  auto kokkosVec_h = Kokkos::create_mirror_view(MemorySpace_Host_t(), kokkosVec);

  // Fill the host mirror
  for (int i = 0; i < n; ++i)
  {
    kokkosVec_h(i) = vec[i];
  }

  // Deep copy into the original view (device if needed)
  Kokkos::deep_copy(kokkosVec, kokkosVec_h);

  return kokkosVec;
}

inline ViewVector_d toKokkos(const std::vector<double>& stdVec)
{
  int n = stdVec.size();
  ViewVector_d kokkosVec("kokkosVec", n);

  // Always create host mirror (may just return kokkosVec itself if already in HostSpace)
  auto kokkosVec_h = Kokkos::create_mirror_view(MemorySpace_Host_t(), kokkosVec);

  // Fill the host mirror
  for (int i = 0; i < n; ++i)
  {
    kokkosVec_h(i) = stdVec[i];
  }

  // Deep copy into the original view (device if needed)
  Kokkos::deep_copy(kokkosVec, kokkosVec_h);

  return kokkosVec;
}

inline ViewMatrix_d toKokkos(const Teuchos::SerialDenseMatrix<int, double>& mat)
{
  int numRows = mat.numRows();
  int numCols = mat.numCols();
  ViewMatrix_d kokkosMat("kokkosMat", numRows, numCols);

  auto kokkosMat_h = Kokkos::create_mirror_view(MemorySpace_Host_t(), kokkosMat);

  for (int i = 0; i < numRows; ++i)
    for (int j = 0; j < numCols; ++j)
      kokkosMat_h(i, j) = mat(i, j);  // Access via (row, col) despite column-major internal layout

  Kokkos::deep_copy(kokkosMat, kokkosMat_h);

  return kokkosMat;
}



inline Teuchos::SerialDenseVector<int, double> kokkosVectorToTeuchos(const ViewVector_d& kokkosVec)
{
  int n = kokkosVec.extent(0);

  auto kokkosVec_h = Kokkos::create_mirror_view_and_copy(MemorySpace_Host_t(), kokkosVec);

  Teuchos::SerialDenseVector<int, double> vec(n);
  for (int i = 0; i < n; ++i)
  {
    vec[i] = kokkosVec_h(i);
  }
  return vec;
}

inline std::vector<double> kokkosVectorToStdVector(const ViewVector_d& kokkosVec)
{
  int n = kokkosVec.extent(0);

  auto kokkosVec_h = Kokkos::create_mirror_view_and_copy(MemorySpace_Host_t(), kokkosVec);

  std::vector<double> stdVec(n);
  for (int i = 0; i < n; ++i)
  {
    stdVec[i] = kokkosVec_h(i);
  }
  return stdVec;
}

inline Teuchos::SerialDenseMatrix<int, double> kokkosMatrixToTeuchos(const ViewMatrix_d& kokkosMat)
{
  int numRows = kokkosMat.extent(0);
  int numCols = kokkosMat.extent(1);

  auto kokkosMat_h = Kokkos::create_mirror_view_and_copy(MemorySpace_Host_t(), kokkosMat);

  Teuchos::SerialDenseMatrix<int, double> mat(numRows, numCols);
  for (int i = 0; i < numRows; ++i)
    for (int j = 0; j < numCols; ++j)
    {
      mat(i, j) = kokkosMat_h(i, j);
    }
  return mat;
}



/*
inline ViewVector_Default toKokkos_old(const Teuchos::SerialDenseVector<int, double>& vec)
{
int n = vec.length();
ViewVector kokkosVec("kokkosVec", n);
for (int i = 0; i < n; ++i)
{
kokkosVec(i) = vec[i];
}
return kokkosVec;



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


}*/
