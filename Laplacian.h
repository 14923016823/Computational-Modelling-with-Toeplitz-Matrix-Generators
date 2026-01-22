#include "BlockToeplitz.h"
#include "DiagonalMatrix.h"

class Laplacian2D_ToeplitzMatrix {
    public:

    Laplacian2D_ToeplitzMatrix(const int rows, const int cols);

    void generateMatrix(const int rows, const int cols);

    Vectord Laplacian(Vectord& input);
    
    private:
    Matrix* Incidence; //the geometry of the problem
    Matrix* Diagonal;//This matrix encapsulates grid spacing and variable k
    Matrix* Incidence_T;//negative transpose of incidence matrix



    double k_func(double x, double y);    

    //Vectord generateMatrix(const int rows, const int cols, Vectord& b1);
};