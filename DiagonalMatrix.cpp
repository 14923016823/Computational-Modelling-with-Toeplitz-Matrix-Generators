#include <vector>
#include <stdexcept>
#include <iostream>
//#include <tuple>
//#include "Vector.h"


#include "DiagonalMatrix.h"


DiagonalMatrix::DiagonalMatrix(int size)
{
    Num_Rows=size;
    Num_Cols=size;
    Diag_Vals=new double[Num_Cols];
}

//copy and scale constructor
DiagonalMatrix::DiagonalMatrix(DiagonalMatrix& other, double c)
{
    Num_Cols=other.Num_Cols;
    Num_Rows=other.Num_Rows;
    Diag_Vals=new double[Num_Cols];
    for(int i=0;i<Num_Cols;i++)
    {
        Diag_Vals[i]=other.Diag_Vals[i]*c;
    }
}

DiagonalMatrix::~DiagonalMatrix()
{
    delete Diag_Vals;
}

// scalar multiplication
void DiagonalMatrix::operator*=(double c) 
{
    for (int i=0;i<Num_Cols;i++)
    {
        Diag_Vals[i]*=c;
    }
}

// matrix-vector multiplication (fast)
Vectord DiagonalMatrix::operator*(Vectord& vec) 
{
    if (vec.len() != Num_Cols)
        throw std::runtime_error("Vector size mismatch");

    Vectord result=Vectord(vec.len());
    for (size_t i = 0; i < Num_Cols; i++) {
        result.Vec[i] = Diag_Vals[i] * vec.Vec[i];
    }
    return result;
}

//Returns a cloned and scaled version
Matrix* DiagonalMatrix::Clone(double c)
{
    DiagonalMatrix* copy = new DiagonalMatrix(*this, c);
    return copy;
}

Matrix* DiagonalMatrix::Kronecker(Matrix&) 
{
    throw std::runtime_error("Kronecker not implemented for DiagonalMatrix");
}

