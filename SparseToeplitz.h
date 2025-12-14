#pragma once

#include "VectorD.h"
#include "Matrix.h"

class SparseToeplitz: public Matrix
{
public:
    int Num_Diags;// the length of the Diags and Vals vectors
    int* Diags;//an int array with all nonzero diagonals (ordered smallest to largest)
    double* Vals;//a double array with the values at those diagonals (ordered the same way as the int array)
    SparseToeplitz(int num_rows, int num_cols, int num_diags);

    SparseToeplitz(int num_rows, int num_cols, int num_diags,int* diags,double* vals);

    SparseToeplitz(SparseToeplitz& other, double c);

    virtual Vectord operator*(Vectord& vec);

    void print();

    void operator*=(double c) override;
   
    Matrix* Kronecker(Matrix& B) override;

    Matrix* Clone(double c) override;
};