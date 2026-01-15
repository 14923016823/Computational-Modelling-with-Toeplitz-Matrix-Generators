#ifndef SPARSE_TOEPLITZ_H
#define SPARSE_TOEPLITZ_H

#include "Matrix.h"
#include "BlockToeplitz.h"
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>

// Header-only templated SparseToeplitz. Default scalar type is double.
template<typename T = double>
class SparseToeplitz: public Matrix
{
public:
    std::vector<int> Diags;
    std::vector<T> Vals;
    int Num_Diags;

    SparseToeplitz(int nrows, int ncols, int ndiags)
    {
        Num_Rows = nrows;
        Num_Cols = ncols;

        if (ndiags > nrows + ncols - 1)
        {
            throw std::invalid_argument("There are too many diagonals for Toeplitz matrix of this size");
        }
        Num_Diags = ndiags;
        Diags.assign(Num_Diags, 0);
        Vals.assign(Num_Diags, T(0));
    }

    SparseToeplitz(SparseToeplitz& other, double c)
    {
        Num_Cols = other.Num_Cols;
        Num_Rows = other.Num_Rows;
        Num_Diags = other.Num_Diags;
        Diags = other.Diags;
        Vals.resize(Num_Diags);
        for (int i = 0; i < Num_Diags; i++)
        {
            Vals[i] = static_cast<T>(c * static_cast<double>(other.Vals[i]));
        }
    }

    SparseToeplitz(int nrows, int ncols, int ndiags, const int* diags, const T* vals)
    {
        Num_Rows = nrows;
        Num_Cols = ncols;

        if (ndiags > nrows + ncols - 1)
        {
            throw std::invalid_argument("There are too many diagonals for Toeplitz matrix of this size");
        }
        Num_Diags = ndiags;
        Diags.resize(Num_Diags);
        Vals.resize(Num_Diags);

        for (int i = 0; i < ndiags; i++)
        {
            Diags[i] = diags[i];
            Vals[i] = vals[i];
        }
    }

    ~SparseToeplitz() = default;

    void print()
    {
        if (rows() <= 10 && cols() <= 10) {
            std::cout << "SparseToeplitz Matrix (" << Num_Rows << " x " << Num_Cols << "):\n";
            for (int i = 0; i < Num_Rows; i++) {
                std::cout << "[";
                for (int j = 0; j < Num_Cols; j++) {
                    double val = 0.0;
                    for (int d = 0; d < Num_Diags; d++) {
                        if (Diags[d] == j - i) {
                            val = static_cast<double>(Vals[d]);
                            break;
                        }
                    }
                    std::cout << std::setw(8) << std::fixed << std::setprecision(3) << val << " ";
                }
                std::cout << "]\n";
            } std::cout << std::endl;
        } else {
            std::cout << "Matrix too large to print. (rows=" << rows() << ", cols=" << cols() << ")\n";
        }
    }

    void matvec(const Vectord& in, Vectord& out) const override
    {
        if (in.len() != Num_Cols || out.len() != Num_Rows)
        {
            throw std::invalid_argument("Vector and Matrix size dont match");
        }
        for (int i = 0; i < Num_Rows; i++)
        {
            out[i] = 0.0;
            for (int j = 0; j < Num_Diags; j++)
            {
                int current_col = Diags[j] + i;
                if (current_col >= 0 && current_col < Num_Cols)
                {
                    out[i] += in[current_col] * static_cast<double>(Vals[j]);
                }
            }
        }
    }

    void operator*=(double c) override
    {
        for (int i = 0; i < Num_Diags; i++)
        {
            Vals[i] = static_cast<T>(static_cast<double>(Vals[i]) * c);
        }
    }

    Matrix* Kronecker(Matrix& B) override
    {
        BlockToeplitz<T>* result = new BlockToeplitz<T>(B.rows() * Num_Rows, B.cols() * Num_Cols, Num_Diags);
        for (int i = 0; i < Num_Diags; i++)
        {
            result->Diags[i] = Diags[i];
            // use double scaling for Clone
            result->Vals[i] = B.Clone(static_cast<double>(Vals[i]));
        }
        return result;
    }

    Matrix* Clone(double c) override
    {
        return new SparseToeplitz(*this, c);
    }
};

// convenience typedef for the common double version
using SparseToeplitzd = SparseToeplitz<double>;

#endif

