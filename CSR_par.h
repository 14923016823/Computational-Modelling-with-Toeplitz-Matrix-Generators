#pragma once

#include "Matrix.h"
#include "VectorD.h"
#include "SparseToeplitz.h"
#include "BlockCSR.h"
#include "COO.h"
#include "BlockToeplitz.h"


class CSR: public Matrix
{
public:
    double* Vals;
    int* Cols;
    int* Rows;
    int Num_Vals;

    //constructor
    void MatMulAdd(const Vectord& x,const int x_start,Vectord& result,const int result_start) override;
    void MatMul(const Vectord& x,Vectord& result) override;
    
    CSR(double* vals, int* cols, int* rows, int num_vals, int num_rows, int num_cols);

    CSR(int num_rows, int num_cols, int num_vals);

    CSR(CSR& other, double c);

    //virtual Vectord operator*(Vectord& vect) override;

    Vectord operator*(Vectord& vec);

    CSR(SparseToeplitz& ST);

    CSR(BlockToeplitz& ST);

    void print();

    //destructor
    ~CSR();

    void operator*=(double c) override;

    Matrix* Kronecker(Matrix& B) override;

    Matrix* Clone(double c) override;

    Matrix* negativeTranspose() override;
    
    void printFullMatrix() override;
    
    double operator()(int i, int j) const override;
};