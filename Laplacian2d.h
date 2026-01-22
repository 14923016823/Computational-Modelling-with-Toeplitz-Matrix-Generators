#include "BlockToeplitz.h"
#include "DiagonalMatrix.h"

class Laplacian2D_ToeplitzMatrix :Matrix{
    public:

    Laplacian2D_ToeplitzMatrix(const int rows, const int cols);

    void generateMatrix(const int rows, const int cols);

    Vectord Laplacian(Vectord& input);
    
    private:
    Matrix* Incidence; //the geometry of the problem
    Matrix* Diagonal;//This matrix encapsulates grid spacing and variable k
    Matrix* Incidence_T;//negative transpose of incidence matrix

    virtual Vectord operator*(Vectord& vec);
    double k_func(double x, double y);    


    void operator*=(double c) override;//all of these are not implemented
    Matrix* Kronecker(Matrix& B) override;
    Matrix* Clone(double c) override;
    Matrix* negativeTranspose() override;
    void printFullMatrix() override;
    double operator()(int i, int j) const override;

    

    //Vectord generateMatrix(const int rows, const int cols, Vectord& b1);
};
