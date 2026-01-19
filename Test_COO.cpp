#include "Test_COO.h"

int main()
{
    int STlength = 3;
    int STdiags[STlength] = {-1,0,1};
    double STvals[STlength] = {-1,2,-1};
    int STwidth = 400;
    int STheight = 400;
    Vectord x(STwidth);
    for(int i=0;i<STwidth;i++)
    {
        x.Vec[i]=i;
    }
    SparseToeplitz C(STheight, STwidth, STlength, STdiags, STvals);

    CSR C_COO(C);
    Vectord a = C_COO*x;

    return 0;
}
