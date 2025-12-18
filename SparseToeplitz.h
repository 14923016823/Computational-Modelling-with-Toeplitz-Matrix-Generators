#pragma once

#include <vector>
#include <stdexcept>
#include <iostream>
#include "Matrix.h"
#include "VectorD.h"

class SparseToeplitz: public Matrix
{
public:
    int* Diags;//an int array with all nonzero diagonals (ordered smallest to largest)
    double* Vals;//a double array with the values at those diagonals (ordered the same way as the int array)
    int Num_Diags;
    SparseToeplitz(int nrows, int ncols, int ndiags);
    SparseToeplitz(SparseToeplitz& other, double c);
    ~SparseToeplitz();
   
    Vectord operator*(Vectord& vec) override;

    void operator*=(double c) override;
   
    Matrix* Kronecker(Matrix& B) override;
    
    Matrix* Clone(double c) override;

    Matrix* negativeTranspose() override;
    
    Matrix* printFullMatrix() override;
};




