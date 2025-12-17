#ifndef SPARSE_TOEPLITZ_H
#define SPARSE_TOEPLITZ_H

#include "Matrix.h"

class SparseToeplitz: public Matrix
{
public:
    int* Diags;//an int array with all nonzero diagonals (ordered smallest to largest)
    double* Vals;//a double array with the values at those diagonals (ordered the same way as the int array)
    int Num_Diags;
    SparseToeplitz(int nrows, int ncols, int ndiags);
    SparseToeplitz(SparseToeplitz& other, double c);
    SparseToeplitz(int nrows, int ncols, int ndiags, int* diags, double* vals);

    ~SparseToeplitz();

    void print();
   
    Vectord operator*(Vectord& vec) override;

    void operator*=(double c) override;
   
    Matrix* Kronecker(Matrix& B) override;
    Matrix* Clone(double c) override;
    
};

#endif


