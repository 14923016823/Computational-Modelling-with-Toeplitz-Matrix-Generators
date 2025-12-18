#pragma once

#include <tuple>

typedef std::tuple<double,int,int> tuple;

#include "Matrix.h"
#include "VectorD.h"
#include "SparseToeplitz.h"
#include "BlockCOO.h"

class COO: public Matrix
{
public:
    tuple* Array;
    int Num_Vals;

    COO(int length, int num_rows, int num_cols);

    COO(const std::initializer_list<tuple>& list, int num_rows, int num_cols);

    COO(COO& other, double c);

    Vectord operator*(Vectord& vect);

    void print();

    COO(SparseToeplitz& ST);

    ~COO();

    void operator*=(double c) override;

    Matrix* Kronecker(Matrix& B) override;

    Matrix* Clone(double c) override;
};