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

    virtual Matrix* negativeTranspose() = 0;

    virtual void printFullMatrix()=0;

    virtual double operator()(int i, int j) const { return 0.0; };

    virtual void MatMulAdd(const Vectord& x, const int x_start,Vectord& result,const int result_start) =0; //result+=Ax
    virtual void MatMul(const Vectord& x,Vectord& result) =0;//result=Ax
};

typedef  Matrix* MatrixPointer;