#pragma once

#include <iostream>
#include <exception>
#include <cmath>
#include <initializer_list>

#include "Matrix.h"
#include "VectorD.h"

class COO: public Matrix
{
public:
    tuple* Array;
    int Length;

    COO(int length, int num_rows, int num_cols);

    COO(const std::initializer_list<tuple>& list, int num_rows, int num_cols);

    virtual Vectord operator*(const Vectord vect) override;

    void print();

    COO(SparseToeplitz ST);

    ~CSR();
};