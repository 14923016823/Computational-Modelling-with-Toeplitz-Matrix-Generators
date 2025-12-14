#pragma once

#include <iostream>
#include <exception>
#include <cmath>
#include <initializer_list>

#include "Matrix.h"
#include "VectorD.h"
#include "SparseToeplitz.h"
#include "BlockCSR.h"


class CSR: public Matrix
{
public:
    double* Vals;
    int* Cols;
    int* Rows;
    int Num_Vals;

    //constructor
    CSR(double* vals, int* cols, int* rows, int num_vals, int num_rows, int num_cols);

    CSR(CSR& other, double c);

    //virtual Vectord operator*(Vectord& vect) override;

    Vectord operator*(Vectord& vec);

    CSR(SparseToeplitz& ST);

    void print();

    //destructor
    ~CSR();

    void operator*=(double c) override;

    Matrix* Kronecker(Matrix& B) override;

    Matrix* Clone(double c) override;
};