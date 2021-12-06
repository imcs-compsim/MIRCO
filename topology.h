#include <Epetra_SerialDenseMatrix.h> // Seems obvious

class TopologyGeneration
{
public:
int n;
    virtual void GetSurface(Epetra_SerialDenseMatrix &z) = 0;
    TopologyGeneration(int nn)
    {
        n = nn;
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
double H;
    void GetSurface(Epetra_SerialDenseMatrix &z) override;
    Rmg(int nn, double HH) : TopologyGeneration(nn)
    {
        H = HH;
    }
};