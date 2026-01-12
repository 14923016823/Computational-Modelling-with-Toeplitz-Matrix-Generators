#pragma once

#include "Matrix.h"
#include "VectorD.h"
#include "BlockToeplitz.h"

class BlockCSR: public Matrix
{
public:
    int Num_Vals;
    int* Rows;
    int* Cols;
    MatrixPointer* Vals;       //sub-CSR matrices (only if recursive)

    // constructor
    BlockCSR(int nrows, int ncols, int nvals);

    //copy and scale constructor
    BlockCSR(BlockCSR& other, double c);
    ~BlockCSR(); 

    Vectord operator*(Vectord& vec) override;
   

    void operator*=(double c) override;
    

    Matrix* Kronecker(Matrix& B) override;

    //Returns a cloned and scaled version of 
    Matrix* Clone(double c) override;
   
    Matrix* negativeTranspose() override;
    
    Matrix* printFullMatrix() override;
    
    double operator()(int i, int j) const override;
};

