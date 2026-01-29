#include "BlockToeplitz.h"
#include "DiagonalMatrix.h"
#include "COO.h"
#include "CSR.h"

class Laplacian2D_ToeplitzMatrix{
    public:

    Laplacian2D_ToeplitzMatrix(const int rows, const int cols);

    ~Laplacian2D_ToeplitzMatrix();

    void generateMatrix(const int rows, const int cols);

    void Laplacian2d(const Vectord& x,Vectord& result);

    COO COO_Laplacian();
    CSR CSR_Laplacian();

    
    private:
    Matrix* Incidence; //the geometry of the problem
    Matrix* Diagonal;//This matrix encapsulates grid spacing and variable k
    Matrix* Incidence_T;//negative transpose of incidence matrix

    double k_func(double x, double y);    
    
    Vectord Tmp1;
    Vectord Tmp2;

    //Vectord generateMatrix(const int rows, const int cols, Vectord& b1);
};
