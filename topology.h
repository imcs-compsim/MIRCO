#include <Epetra_SerialDenseMatrix.h> // Seems obvious

class TopologyGeneration
{
public:
int resolution; // resolution parameter
    virtual void GetSurface(Epetra_SerialDenseMatrix &z) = 0;
    TopologyGeneration(int nn)
    {
        resolution = nn;
    }
};

class ReadFile : public TopologyGeneration
{
public:
string filepath;
    void GetSurface(Epetra_SerialDenseMatrix &z) override;
    ReadFile(int nn, string ffilepath) : TopologyGeneration(nn)
    {
        filepath = ffilepath;
    }
};

class Rmg : public TopologyGeneration
{
public:
double Hurst; // Hurst component
bool rand_seed_flag;
    void GetSurface(Epetra_SerialDenseMatrix &z) override;
    Rmg(int nn, double HH, bool rsf) : TopologyGeneration(nn)
    {
        Hurst = HH;
        rand_seed_flag = rsf;
    }
};