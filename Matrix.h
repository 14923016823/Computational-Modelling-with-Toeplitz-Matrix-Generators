#ifndef MATRIX_H
#define MATRIX_H

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
    Matrix(int cols, int rows);

    virtual ~Matrix()=default;

    // scalar multiplication
    virtual void operator*=(double scalar) = 0;

    // vector multiplication
    virtual Vectord operator*(const Vectord& vec) const;
    virtual void matvec(const Vectord& in, Vectord& out) const = 0;

    virtual Matrix* Clone(double c) = 0;
 

    virtual Matrix* Kronecker(Matrix&) = 0;
};

typedef Matrix* MatrixPointer;

#endif