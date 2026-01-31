#include "Test_L3_CSR.h"

int main()
{  
    //Change the dimensions as needed
    int rows = 1000;
    int cols = 1000;
    int arrays = 1000;

    Vectord vec3(cols*rows*arrays);
    for (int i = 0; i < rows * cols * arrays; ++i) {
        vec3.Vec[i] = i; // Example initialization
    }

    Laplacian3D_ToeplitzMatrix L3d=Laplacian3D_ToeplitzMatrix(rows,cols,arrays);
    std::cout << "Laplacian generated" << std::endl;
    Vectord a(rows*cols*arrays);
    CSR L3d_CSR = L3d.CSR_Laplacian();
    //Change the amount of loops as needed
    for(int i=0;i<10;i++)
        std::cout << "Iteration " << i << std::endl;
        L3d_CSR.MatMul(vec3,a);

    return 0;

}
