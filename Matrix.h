#pragma once

#include "Vectord.h"

class Matrix 
{
protected:
    int Num_Rows;
    int Num_Cols;

public:
    int rows() const;
    int cols() const;

    Matrix();
    //Matrix(int cols, int rows);

    //virtual ~Matrix();

    // scalar multiplication
    virtual void operator*=(double scalar) = 0;

    // vector multiplication
    virtual Vectord operator*(Vectord& vec) = 0;

    virtual Matrix* Clone(double c) = 0;
 

    virtual Matrix* Kronecker(Matrix&) = 0;

    virtual Matrix* negativeTranspose() = 0;

    virtual Matrix* printFullMatrix() = 0;

    virtual double operator()(int i, int j) const { return 0.0; }
};

typedef  Matrix* MatrixPointer;