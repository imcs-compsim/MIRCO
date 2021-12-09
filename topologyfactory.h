#include <string>
#include <memory>
#include <Epetra_SerialSymDenseMatrix.h>

void SortSurf(int n, double Hurst, bool rand_seed_flag, std::string zfilePath, bool rmg_flag,std::shared_ptr<TopologyGeneration>& surfacegenerator);