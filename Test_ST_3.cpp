#include "Test_ST_2.h"

int main()
{
    int STlength = 100;
    int STdiags[STlength];
    for(int i=0;i<STlength;i++)
    {
        STdiags[i] = i+1;
    }

    double STvals[STlength];
    int n=0;
    for(int i=-STlength/2;i<STlength/2;i++)
    {
        STvals[n] = 2*i;
        n++;
    }

    int STwidth = 400;
    int STheight = 400;
    Vectord x(STwidth);
    for(int i=0;i<STwidth;i++)
    {
        x.Vec[i]=i;
    }
    SparseToeplitz C(STheight, STwidth, STlength, STdiags, STvals);

    Vectord a = C*x;

    return 0;
}
