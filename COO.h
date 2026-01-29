#pragma once

#include <tuple>

typedef std::tuple<int,int,double> tuple;

#include "Matrix.h"
#include "VectorD.h"
#include "SparseToeplitz.h"
#include "BlockCOO.h"
#include "BlockToeplitz.h"

class COO: public Matrix
{
public:
    tuple* Array;
    int Num_Vals;

    COO(int num_rows, int num_cols, int num_vals);

    COO(int num_rows, int num_cols, const std::initializer_list<tuple>& list);

    COO(COO& other, double c);

    Vectord operator*(Vectord& vect);


    void print();

    COO(SparseToeplitz& ST);

    COO(BlockToeplitz& ST);

    ~COO();

    void operator*=(double c) override;

    Matrix* Kronecker(Matrix& B) override;

    Matrix* Clone(double c) override;

    Matrix* negativeTranspose() override;
    
    void printFullMatrix() override;
    
    double operator()(int i, int j) const override;
};