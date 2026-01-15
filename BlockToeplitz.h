#ifndef BLOCKTOEPLITZ_H
#define BLOCKTOEPLITZ_H

#include "Matrix.h"
#include <vector>

struct DiagEntry {
    int offset;
    const MatrixPointer blockptr;
};

template<typename T = double>
class BlockToeplitz: public Matrix
{
public:
    int Num_Diags;
    std::vector<int> Diags;
    std::vector<MatrixPointer> Vals;       //sub-Toeplitz matrices (only if recursive)

    // constructor
    BlockToeplitz(int nrows, int ncols, int ndiags);

    // constructor that initializes diagonal offsets and pointers from each diagonal
    template<int N>
    BlockToeplitz(int nrows, int ncols, const DiagEntry (&list)[N])
        : Matrix(nrows, ncols)
    {
        Num_Diags = N;
        Diags.resize(Num_Diags);
        Vals.resize(Num_Diags);
        for (int i = 0; i < Num_Diags; ++i) {
            Diags[i] = list[i].offset;
            Vals[i] = list[i].blockptr;
        }
    }

    //copy and scale constructor
    BlockToeplitz(BlockToeplitz& other, double c);
    ~BlockToeplitz() = default;

    void set_diag(int idx, int offset, Matrix* ptr);

    void matvec(const Vectord& in, Vectord& out) const override;

    void operator*=(double c) override;


    Matrix* Kronecker(Matrix& B) override;

    //Returns a cloned and scaled version of
    Matrix* Clone(double c) override;

};

using BlockToeplitzd = BlockToeplitz<double>;

#endif
