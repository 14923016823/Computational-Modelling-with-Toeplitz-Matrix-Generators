#ifndef DIAGONAL_MATRIX_H
#define DIAGONAL_MATRIX_H

#include "Matrix.h"
#include <vector>
#include <stdexcept>

// Templated diagonal matrix (header-only). Default element type is double.
template<typename T = double>
class DiagonalMatrix : public Matrix
{
public:
    std::vector<T> Diag_Vals;

    DiagonalMatrix(int size)
    {
        Num_Rows = size;
        Num_Cols = size;
        Diag_Vals.assign(Num_Cols, T(0));
    }

    // copy and scale constructor
    DiagonalMatrix(DiagonalMatrix& other, double c)
    {
        Num_Cols = other.Num_Cols;
        Num_Rows = other.Num_Rows;
        Diag_Vals.resize(Num_Cols);
        for (int i = 0; i < Num_Cols; ++i)
            Diag_Vals[i] = static_cast<T>(c * static_cast<double>(other.Diag_Vals[i]));
    }

    // scalar multiplication
    void operator*=(double c) override
    {
        for (int i = 0; i < Num_Cols; ++i)
            Diag_Vals[i] = static_cast<T>(static_cast<double>(Diag_Vals[i]) * c);
    }

    // matrix-vector multiplication (fast)
    void matvec(const Vectord& in, Vectord& out) const override
    {
        if (in.len() != Num_Cols || out.len() != Num_Rows)
            throw std::runtime_error("Vector size mismatch");

        for (int i = 0; i < Num_Cols; ++i) {
            out[i] = static_cast<double>(Diag_Vals[i]) * in[i];
        }
    }

    //Returns a cloned and scaled version
    Matrix* Clone(double c) override
    {
        DiagonalMatrix* copy = new DiagonalMatrix(*this, c);
        return copy;
    }

    Matrix* Kronecker(Matrix&) override
    {
        throw std::runtime_error("Kronecker not implemented for DiagonalMatrix");
    }
};

using DiagonalMatrixd = DiagonalMatrix<double>;

#endif