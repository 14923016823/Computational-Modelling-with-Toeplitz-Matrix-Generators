#include <iostream>
#include <exception>
#include <cmath>
#include <initializer_list>

#include "Matrix.h"
#include "VectorD.h"

class CSR: public Matrix
{
public:
    double* Vals;
    int* Cols;
    int* Rows;
    int Length;

    //constructor
    CSR(double* vals, int* cols, int* rows, int length, int num_rows, int num_cols);

    virtual Vectord operator*(const Vectord vect) override;

    void print();

    //destructor
    ~CSR();
};