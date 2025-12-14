#pragma once

#include <vector>
#include <stdexcept>
#include <iostream>
#include <tuple>

#include "Matrix.h"
#include "VectorD.h"

typedef std::tuple<MatrixPointer,int,int> mtuple;

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
   
};

