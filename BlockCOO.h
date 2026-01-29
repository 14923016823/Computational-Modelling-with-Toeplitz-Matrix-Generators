#pragma once

#include <tuple>

#include "Matrix.h"
#include "VectorD.h"
#include "BlockToeplitz.h"

typedef std::tuple<int,int,MatrixPointer> mtuple;

class BlockCOO: public Matrix
{
public:
    mtuple* Array;
    int Num_Vals;
    //MatrixPointer* Vals;       //sub-COO matrices (only if recursive)

    // constructor
    BlockCOO(int nrows, int ncols, int nvals);

    //copy and scale constructor
    BlockCOO(BlockCOO& other, double c);
    ~BlockCOO(); 

    Vectord operator*(Vectord& vec) override;
   

    void operator*=(double c) override;
    

    Matrix* Kronecker(Matrix& B) override;

    //Returns a cloned and scaled version of 
    Matrix* Clone(double c) override;
   
    Matrix* negativeTranspose() override;
    
    void printFullMatrix() override;
    
    double operator()(int i, int j) const override;

    void MatMulAdd(const Vectord& x,const int x_start,Vectord& result,const int result_start) override;
    void MatMul(const Vectord& x,Vectord& result) override;
};

