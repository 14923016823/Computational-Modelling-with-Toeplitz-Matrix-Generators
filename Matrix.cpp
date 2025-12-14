#include "Matrix.h"

int Matrix::rows() const { return Num_Rows; }
int Matrix::cols() const { return Num_Cols; }

Matrix::Matrix()
{ 
    Num_Cols = 0;
    Num_Rows = 0;
}

Matrix::~Matrix(){}

Vectord Matrix::operator*(Vectord& vec) 
{
    return Vectord(0);
}

void Matrix::operator*=(double scalar) 
{ }

void Matrix::print(){}
