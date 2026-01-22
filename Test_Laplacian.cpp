#include "Laplacian.h"



void Test()
{
    int xcells=1000;
    int ycells=1000;
    Laplacian2D_ToeplitzMatrix Periodic=Laplacian2D_ToeplitzMatrix(xcells,ycells);
    Vectord x=Vectord(xcells*ycells);
    for(int i=0;i<xcells*ycells;i++)
    {
        x.Vec[i]=1.0;
    }
    







}