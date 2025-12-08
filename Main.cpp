#include <vector>
#include <stdexcept>
#include <iostream>
//#include <tuple>
//#include "Vector.h"

#include "Main.h"

void TestKronecker() //This function tests the kornicker function by generating 2d laplacian with dirichlet boundary conditions
{
    int n=3;
    SparseToeplitz L1d=SparseToeplitz(n,n,3); //Laplacian1d
    L1d.Diags[0]=-1;
    L1d.Diags[1]=0;
    L1d.Diags[2]=1;
    L1d.Vals[0]=-1;
    L1d.Vals[1]=2;
    L1d.Vals[2]=-1;
    SparseToeplitz I=SparseToeplitz(n,n,1); //identity matrix
    I.Diags[0]=0;
    I.Vals[0]=1;
    Matrix* l2d=I.Kronecker(L1d);
 
    printf("%d",(L1d).rows());
    printf("%d",(I).rows());

    printf("%d",(*l2d).rows());
   

    Vectord input(n*n);
    input.Vec[0] = 1.0;
    input.PrintVector();
    Vectord result = (*l2d)*input;
    result.PrintVector();
    for(int i=1;i<(n*n);i++)
    {
        input.Vec[i-1]=0;
        input.Vec[i]=1.0;
        Vectord result = (*l2d)*input;
        result.PrintVector();
    }
};

int main() 
{
    TestKronecker(); 

    return 0;
}
