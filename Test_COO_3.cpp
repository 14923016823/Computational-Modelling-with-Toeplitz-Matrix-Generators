#include "Test_COO_3.h"

int main()
{
    int STlength = 1000;
    int STdiags[STlength];
    for(int i=0;i<STlength;i++)
    {
        STdiags[i] = 2*i;
    }

    double STvals[STlength];
    int n=0;
    for(int i=-STlength/2;i<STlength/2;i++)
    {
        STvals[n] = (double)(2*i);
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

    COO C_COO(C);
    for(int n=0;n<1000;n++)
    Vectord a = C_COO*x;

    return 0;
}
