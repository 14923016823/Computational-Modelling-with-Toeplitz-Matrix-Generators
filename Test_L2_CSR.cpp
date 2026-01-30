#include "Test_L2_CSR.h"

int main()
{  
    //Change the dimensions as needed
    int rows = 10000;
    int cols = 10000;

    Vectord vec2(cols*rows);
    for (int i = 0; i < rows * cols; ++i) {
        vec2.Vec[i] = i; // Example initialization
    }
    return 0;

    Laplacian2D_ToeplitzMatrix L2d=Laplacian2D_ToeplitzMatrix(rows,cols);
    Vectord a(rows*cols);
    CSR L2d_CSR = L2d.CSR_Laplacian();
    //Change the amount of loops as needed
    for(int i=0;i<10;i++)
        L2d_CSR.MatMul(vec2,a);
}