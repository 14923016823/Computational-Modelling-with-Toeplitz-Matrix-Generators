#include "BlockToeplitz.h"
#include "DiagonalMatrix.h"

class Laplacian2D_ToeplitzMatrix {
    public:

    Laplacian2D_ToeplitzMatrix(const int rows, const int cols, Vectord& b1);

    void generateMatrix(const int rows, const int cols, Vectord& b1);

    Vectord Laplacian(Vectord& input);
    
    private:
    Matrix* Incidence; //the geometry of the problem
    Matrix* Diagonal;//This matrix encapsulates grid spacing and variable k
    Matrix* Incidence_T;//negative transpose of incidence matrix



    double k_func(double x, double y);    

    //Vectord generateMatrix(const int rows, const int cols, Vectord& b1);
};

class Laplacian3D_ToeplitzMatrix {
    public:

    Laplacian3D_ToeplitzMatrix(const int rows, const int cols, const int arrays, Vectord& b1);

    void generateMatrix(const int rows, const int cols, const int arrays, Vectord& b1);

    Vectord Laplacian(Vectord& input);
    
    private:
    Matrix* Incidence; //the geometry of the problem
    Matrix* Diagonal;//This matrix encapsulates grid spacing and variable k
    Matrix* Incidence_T;//negative transpose of incidence matrix



    double k_func(double x, double y, double z);    

    //Vectord generateMatrix(const int rows, const int cols, Vectord& b1);
};