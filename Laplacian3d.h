#pragma once

#include "BlockToeplitz.h"
#include "DiagonalMatrix.h"
#include "COO.h"
#include "CSR.h"

class Laplacian3D_ToeplitzMatrix{
    public:

    Laplacian3D_ToeplitzMatrix(const int rows, const int cols, const int arrays);

    ~Laplacian3D_ToeplitzMatrix();

    void generateMatrix(const int rows, const int cols, const int arrays);

    void Laplacian3d(const Vectord& x,Vectord& result);

    COO COO_Laplacian();
    CSR CSR_Laplacian();
    
    //private:
    Matrix* Incidence; //the geometry of the problem
    Matrix* upperIncidenceMatrix_3D;
    Matrix* lowerIncidenceMatrix_3D;
    Matrix* Diagonal;//This matrix encapsulates grid spacing and variable k
    Matrix* Incidence_T;//negative transpose of incidence matrix
    Matrix* leftIncidenceT_Matrix_3D;
    Matrix* rightIncidenceT_Matrix_3D;

    double k_func(double x, double y, double z);    

    Vectord Tmp1;
    Vectord Tmp2;

    //Vectord generateMatrix(const int rows, const int cols, Vectord& b1);
};
