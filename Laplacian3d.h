#pragma once

#include "BlockToeplitz.h"
#include "DiagonalMatrix.h"
#include "COO.h"
#include "CSR.h"

class Laplacian3D_ToeplitzMatrix :Matrix{
    public:

    Laplacian3D_ToeplitzMatrix(const int rows, const int cols, const int arrays);

    void generateMatrix(const int rows, const int cols, const int arrays);

    virtual Vectord operator*(Vectord& vec);
    
    private:
    Matrix* Incidence; //the geometry of the problem
    Matrix* Diagonal;//This matrix encapsulates grid spacing and variable k
    Matrix* Incidence_T;//negative transpose of incidence matrix

    double k_func(double x, double y, double z);    


    void operator*=(double c) override;//all of these are not implemented
    Matrix* Kronecker(Matrix& B) override;
    Matrix* Clone(double c) override;
    Matrix* negativeTranspose() override;
    void printFullMatrix() override;
    double operator()(int i, int j) const override;

    COO COO_Laplacian();
    CSR CSR_Laplacian();

    //Vectord generateMatrix(const int rows, const int cols, Vectord& b1);
};
