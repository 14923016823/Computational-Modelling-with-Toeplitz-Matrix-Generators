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
    printf("%f\n",L1d.Vals[0]);
    SparseToeplitz I=SparseToeplitz(n,n,1); //identity matrix
    I.Diags[0]=0;
    I.Vals[0]=1;
    printf("%f\n",I.Vals[0]);
    Matrix* l2d=I.Kronecker(L1d);
    printf("crashes here");
 
    //printf("%n",(L1d).rows());
    //printf("%d",(I).rows());

    //printf("%d\n",(*l2d).rows());
   

    Vectord input(n*n);
    input.Vec[0] = 1.0;
    input.PrintVector();
    fflush(stdout);
    Vectord result = (*l2d)*input;
    result.PrintVector();
    for(int i=1;i<(n*n);i++)
    {
        input.Vec[i-1]=0; 
        input.Vec[i]=1.0;
        result = (*l2d)*input;
        result.PrintVector();
    }
};

void TestMatmul()
{
    Vectord input(12);
    // default all zeros;

    
    //input.PrintVector();

    SparseToeplitz inc_1 = SparseToeplitz(4, 3, 2);
    SparseToeplitz inc0   = SparseToeplitz(4, 3, 2);
    SparseToeplitz inc1   = SparseToeplitz(4, 3, 2);

    inc_1.Vals[0] = 5.0;
    inc_1.Vals[1] = 5.0;
    inc0.Vals[0] = -1.0;
    inc0.Vals[1] = 1.0;
    inc1.Vals[0] = 2.0;
    inc1.Vals[1] = 4.0;
    inc_1.Diags[0]  = -1;
    inc_1.Diags[1]  = 0;
    inc0.Diags[0]  = -1;
    inc0.Diags[1]  = 0;
    inc1.Diags[0]  =-1;
    inc1.Diags[1]  = 0;
    


    BlockToeplitz inc = BlockToeplitz(16, 12, 3);
    inc.Vals[0] = &inc_1;
    inc.Vals[1] = &inc0;
    inc.Vals[2] = &inc1;
    inc.Diags[0] = -1;
    inc.Diags[1] = 0;
    inc.Diags[2] = 1;

    input.Vec[0] = 1.0;
    //input.Vec[1] = 1.0;
    input.PrintVector();
    fflush(stdout);
    Vectord result = inc*input;
    result.PrintVector();
    
    for(int i=1;i<12;i++)
    {
        input.Vec[i-1]=0;
        input.Vec[i]=1.0;
        Vectord result = inc*input;
        result.PrintVector();
    }

}

void TestSPMUL()
{
    Vectord input(3);
    SparseToeplitz inc_1 = SparseToeplitz(4, 3, 2);
    inc_1.Vals[0] =-5.0;
    inc_1.Vals[1] = 5.0;
    inc_1.Diags[0]  = -1;
    inc_1.Diags[1]  = 0;

    input.Vec[0]=1.0;
    Vectord output=inc_1*input;
    output.PrintVector();
    input.Vec[0]=0;
    input.Vec[1]=1.0;
    output=inc_1*input;
    output.PrintVector();
    input.Vec[1]=0;
    input.Vec[2]=1.0;
    output=inc_1*input;
    output.PrintVector();


}

int main() 
{
    //TestSPMUL();
    //TestMatmul();
    TestKronecker(); 
    return 0;
}
