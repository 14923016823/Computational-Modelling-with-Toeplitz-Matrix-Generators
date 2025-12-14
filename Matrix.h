#pragma once
#include "VectorD.h"

class Matrix
{
protected:
    int Num_Rows;
    int Num_Cols;
public:
    int rows() const;
    int cols() const;

    Matrix();
    
    virtual ~Matrix();

    virtual Vectord operator*(Vectord& vec);

    virtual void print();

    virtual Matrix* Clone(double c) = 0;
 
    virtual Matrix* Kronecker(Matrix&) = 0;

    virtual void operator*=(double scalar) = 0;
};

typedef  Matrix* MatrixPointer;