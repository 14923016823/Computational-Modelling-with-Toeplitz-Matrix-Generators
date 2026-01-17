#include "BlockToeplitz.h"

#include <stdexcept>
#include <iostream>

// default constructor
template<typename T>
BlockToeplitz<T>::BlockToeplitz(int nrows, int ncols, int ndiags)
{
    Num_Rows = nrows;
    Num_Cols = ncols;
    Num_Diags = ndiags;
    Diags.assign(Num_Diags, 0);
    Vals.assign(Num_Diags, nullptr);
}

// A helper function to set diagonal offsets to Diags[] and type Matrix* to Vals[]
template<typename T>
void BlockToeplitz<T>::set_diag(int idx, int offset, Matrix* ptr) {
    Diags[idx] = offset;
    Vals[idx] = ptr;
}

//copy and scale constructor
template<typename T>
BlockToeplitz<T>::BlockToeplitz(BlockToeplitz& other, double c)
{

    Num_Rows = other.Num_Rows;
    Num_Cols = other.Num_Rows;
    Num_Diags = other.Num_Diags;
    Diags = other.Diags;
    Vals.resize(Num_Diags);
    for (int i = 0; i < Num_Diags; i++)
    {
        Vals[i] = other.Vals[i]->Clone(c);
    }

}

template<typename T>
void BlockToeplitz<T>::matvec(const Vectord& in, Vectord& out) const
{

    if (in.len() != Num_Cols || out.len() != Num_Rows)
    {
        throw std::invalid_argument("Vector length and matrix columns don't match (block_toeplitz).");
    }

    int Num_Rows_SubMatrixes = Vals[0]->rows();
    int Num_Cols_SubMatrixes = Vals[0]->cols();

    int Num_blockrows = Num_Rows / Num_Rows_SubMatrixes;
    int Num_blockcols = Num_Cols / Num_Cols_SubMatrixes;

    for (int blockrow = 0; blockrow < Num_blockrows; blockrow++)
    {
        Vectord subres(Num_Rows_SubMatrixes); // initialized to zeros
        // accumulate contributions from each diagonal/block
        for (int j = 0; j < Num_Diags; j++)
        {
            int blockcol = Diags[j] + blockrow;
            // ensure the whole sub-block fits in input vector
            if (blockcol >= 0 && blockcol < Num_blockcols)
            {
                Vectord subinput(Num_Cols_SubMatrixes);
                int col_start = blockcol * Num_Cols_SubMatrixes;
                for (int k = 0; k < Num_Cols_SubMatrixes; ++k)
                {
                    subinput[k] = in[col_start + k];
                }
                Vectord subsubres(Num_Rows_SubMatrixes);
                Vals[j]->matvec_fft(subinput, subsubres);
                subres.axpy(1.0, subsubres);
            }
        }
        // copy accumulated block into out
        int row_start = blockrow * Num_Rows_SubMatrixes;
        for (int k = 0; k < Num_Rows_SubMatrixes; ++k)
        {
            out[row_start + k] = subres[k];
        }
    }
}

template<typename T>
void BlockToeplitz<T>::regular_matvec(const Vectord& in, Vectord& out) const
{

    if (in.len() != Num_Cols || out.len() != Num_Rows)
    {
        throw std::invalid_argument("Vector length and matrix columns don't match (block_toeplitz).");
    }

    int Num_Rows_SubMatrixes = Vals[0]->rows();
    int Num_Cols_SubMatrixes = Vals[0]->cols();

    int Num_blockrows = Num_Rows / Num_Rows_SubMatrixes;
    int Num_blockcols = Num_Cols / Num_Cols_SubMatrixes;

    for (int blockrow = 0; blockrow < Num_blockrows; blockrow++)
    {
        Vectord subres(Num_Rows_SubMatrixes); // initialized to zeros
        // accumulate contributions from each diagonal/block
        for (int j = 0; j < Num_Diags; j++)
        {
            int blockcol = Diags[j] + blockrow;
            // ensure the whole sub-block fits in input vector
            if (blockcol >= 0 && blockcol < Num_blockcols)
            {
                Vectord subinput(Num_Cols_SubMatrixes);
                int col_start = blockcol * Num_Cols_SubMatrixes;
                for (int k = 0; k < Num_Cols_SubMatrixes; ++k)
                {
                    subinput[k] = in[col_start + k];
                }
                Vectord subsubres = Vals[j]->operator*(subinput);
                subres.axpy(1.0, subsubres);
            }
        }
        // copy accumulated block into out
        int row_start = blockrow * Num_Rows_SubMatrixes;
        for (int k = 0; k < Num_Rows_SubMatrixes; ++k)
        {
            out[row_start + k] = subres[k];
        }
    }
}

template<typename T>
void BlockToeplitz<T>::operator*=(double c)
{

    for (int i = 0; i < Num_Diags; i++)
    {
        (*Vals[i]) *= c;
    }

}

template<typename T>
Matrix* BlockToeplitz<T>::Kronecker(Matrix& B)//if you add a new matrix at the bottom of the chain every Num_Rows needs to be changed
{

    BlockToeplitz* result = new BlockToeplitz(B.rows()*Num_Rows,B.cols()*Num_Cols,Num_Diags);
    std::cout << "num diags=" << Num_Diags << "\n";
    for (int i = 0; i < Num_Diags; ++i)
    {
        result->Diags[i] = Diags[i];
        result->Vals[i] = Vals[i]->Kronecker(B);
    }
    return result;
}

//Returns a cloned and scaled version of the input matrix
template<typename T>
Matrix* BlockToeplitz<T>::Clone(double c)
{
    BlockToeplitz* copy = new BlockToeplitz(*this, c);
    return copy;
}

// Explicit instantiation for double
template class BlockToeplitz<double>;

