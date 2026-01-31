#include "Test_L3_COO.h"

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
    std::cout << "Laplacian generated" << std::endl;
    Vectord a(rows*cols*arrays);
    COO L3d_COO = L3d.COO_Laplacian();
    //Change the amount of loops as needed
    for(int i=0;i<10;i++)
        std::cout << "Iteration " << i << std::endl;
        L3d_COO.MatMul(vec3,a);
    return 0;

}
