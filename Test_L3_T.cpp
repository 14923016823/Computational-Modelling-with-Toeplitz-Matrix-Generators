#include "Test_L3_T.h"

int main()
{  
    //Change the dimensions as needed
    int rows = 100;
    int cols = 100;
    int arrays = 100;

    Vectord vec3(cols*rows*arrays);
    for (int i = 0; i < rows * cols * arrays; ++i) {
        vec3.Vec[i] = i; // Example initialization
    }

    Laplacian3D_ToeplitzMatrix L3d=Laplacian3D_ToeplitzMatrix(rows,cols,arrays);
    Vectord a(rows*cols);
    //Change the amount of loops as needed
    for(int i=0;i<10;i++)
        L3d.Laplacian3d(vec3,a);

    return 0;
}