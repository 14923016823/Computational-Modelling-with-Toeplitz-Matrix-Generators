#include "Test_L_CSR.h" 

int main()
{    
    int rows = 30000;
    int cols = 40000;
    int arrays = 1;

    Vectord vec(cols*rows*arrays);
    for (int i = 0; i < rows * cols * arrays; ++i) {
        vec.Vec[i] = 1.0; // Example initialization
    }
    
    //Laplacian2D_FullMatrix laplacian(rows, cols);

    Laplacian2D_CSRMatrix L2d=Laplacian2D_CSRMatrix(rows,cols);
    for(int i=0;i<1000;i++)
    {
    try
    {
        Vectord a = L2d*vec;
        //a.print();
    }
    catch(const char* msg)
    {
        std::cout << msg <<std::endl;
    }
    }
    return 0;
}