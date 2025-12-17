#ifndef BLOCKTOEPLITZ_H
#define BLOCKTOEPLITZ_H

#include "Matrix.h"

class BlockToeplitz: public Matrix
{
public:
    int Num_Diags;
    int* Diags;
    MatrixPointer* Vals;       //sub-Toeplitz matrices (only if recursive)

    // constructor
    BlockToeplitz(int nrows, int ncols, int ndiags);

    //copy and scale constructor
    BlockToeplitz(BlockToeplitz& other, double c);
    ~BlockToeplitz(); 

    Vectord operator*(Vectord& vec) override;
   

    void operator*=(double c) override;
    

    Matrix* Kronecker(Matrix& B) override;

    //Returns a cloned and scaled version of 
    Matrix* Clone(double c) override;
   
};

#endif
