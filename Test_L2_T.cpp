#include "Test_L2_T.h"

int main()
{  
    //Change the dimensions as needed
    int rows = 3;
    int cols = 3;

    Vectord vec2(cols*rows);
    for (int i = 0; i < rows * cols; ++i) {
        vec2.Vec[i] = i; // Example initialization
    }

    Laplacian2D_ToeplitzMatrix L2d=Laplacian2D_ToeplitzMatrix(rows,cols);
    Vectord a(rows*cols);
    //Change the amount of loops as needed
    L2d.Laplacian2d(vec2,a);
    a.print();
    return 0;
}