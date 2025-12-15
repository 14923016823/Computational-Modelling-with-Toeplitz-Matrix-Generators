#pragma once

#include <vector>
#include <stdexcept>
#include <iostream>
#include "Matrix.h"
#include "Vectord.h"

class DiagonalMatrix : public Matrix 
{
public:

    double* Diag_Vals;
    DiagonalMatrix(int size);
     ~DiagonalMatrix(); 
      
    //copy and scale constructor
    DiagonalMatrix(DiagonalMatrix& other,double c);

    // scalar multiplication
    void operator*=(double c) override;
    
    // matrix-vector multiplication (fast)
    Vectord operator*(Vectord& vec) override;
    
    //Returns a cloned and scaled version
    Matrix* Clone(double c) override;
   

    Matrix* Kronecker(Matrix&) override;

};


