#include "Test_ST.h"

int main()
{
    int STlength = 3;
    int STdiags[STlength] = {-1,0,1};
    double STvals[STlength] = {-1,2,-1};
    int STwidth = 40000;
    int STheight = 40000;
    Vectord x(STwidth);
    for(int i=0;i<STwidth;i++)
    {
        x.Vec[i]=i;
    }
    SparseToeplitz C(STheight, STwidth, STlength, STdiags, STvals);

    Vectord a = C*x;

    return 0;
}
