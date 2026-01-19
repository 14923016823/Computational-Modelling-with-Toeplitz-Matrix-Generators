#include "BlockToeplitz.h"
#include "DiagonalMatrix.h"

class Laplacian2D_ToeplitzMatrix {
    public:
    Laplacian2D_ToeplitzMatrix(const int rows, const int cols, Vectord& b1);

    Vectord generateMatrix(const int rows, const int cols, Vectord& b1);
    
    private:
    double w(double k, double a);

    double k_func(double x, double y);    

    //Vectord generateMatrix(const int rows, const int cols, Vectord& b1);
};