#include "Test_CSR_2.h"

int main()
{
    int STlength = 10;
    int STdiags[STlength];
    for(int i=0;i<STlength;i++)
    {
        STdiags[i] = 2*i;
    }

    double STvals[STlength];
    int n=0;
    for(int i=-STlength/2;i<STlength/2;i++)
    {
        STvals[n] = 2*i;
        n++;
    }

    int STwidth = 40000;
    int STheight = 40000;
    Vectord x(STwidth);
    for(int i=0;i<STwidth;i++)
    {
        x.Vec[i]=i;
    }
    SparseToeplitz C(STheight, STwidth, STlength, STdiags, STvals);

    CSR C_CSR(C);

    Vectord a = C_CSR*x;

    return 0;
}
