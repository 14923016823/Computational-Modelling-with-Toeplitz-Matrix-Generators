#include "SparseToeplitz.h"
#include "BlockToeplitz.h"

#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>

SparseToeplitz::SparseToeplitz(int nrows, int ncols, int ndiags)
{
    Num_Rows=nrows;
    Num_Cols=ncols;

    if (ndiags > nrows + ncols - 1)
    {
        throw std::invalid_argument("There are too many diagonals for Toeplitz matrix of this size");
    }
    Num_Diags=ndiags;
    Diags = new int[ndiags];
    Vals = new double[ndiags];
}

SparseToeplitz::SparseToeplitz(SparseToeplitz& other, double c)
{
    Num_Cols = other.Num_Cols;
    Num_Rows=other.Num_Rows;
    Num_Diags=other.Num_Diags;
    Diags = new int[Num_Diags];
    Vals = new double[Num_Diags];
    for(int i=0;i<Num_Diags;i++)
    {
        Diags[i]=other.Diags[i];
        Vals[i]=c*other.Vals[i];
    }
}

SparseToeplitz::SparseToeplitz(int nrows, int ncols, int ndiags, int* diags, double* vals)
{
    Num_Rows=nrows;
    Num_Cols=ncols;

    if (ndiags > nrows + ncols - 1)
    {
        throw std::invalid_argument("There are too many diagonals for Toeplitz matrix of this size");
    }
    Num_Diags=ndiags;
    Diags = new int[ndiags];
    Vals = new double[ndiags];

    for (int i=0;i<ndiags;i++) 
    {
        Diags[i]=diags[i];
        Vals[i]=vals[i];
    }
}




SparseToeplitz::~SparseToeplitz()
{
    // cleanup arrays allocated with new[]
    delete[] Vals;
    delete[] Diags;
}


void SparseToeplitz::print()
{
    if (rows() <= 10 && cols() <= 10) {
        std::cout << "SparseToeplitz Matrix (" << Num_Rows << " x " << Num_Cols << "):\n";
        for (int i = 0; i < Num_Rows; i++) {
            std::cout << "[";
            for (int j = 0; j < Num_Cols; j++) {
                double val = 0.0;
                for (int d = 0; d < Num_Diags; d++) {
                    if (Diags[d] == j - i) {
                        val = Vals[d];
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
 
Vectord SparseToeplitz::operator*(Vectord& vec)
{

    if (vec.len() != Num_Cols) 
    {
        throw std::invalid_argument("Vector and Matrix size dont match");
    }
    Vectord result = Vectord(Num_Rows);
    for (int i = 0;i < Num_Rows;i++) 
    {
        for (int j = 0;j < Num_Diags;j++) 
        {
            int current_col=Diags[j] + i;
            if (current_col >= 0 && current_col < Num_Cols)
            {
                result.Vec[i] += vec.Vec[current_col] * Vals[j];
            }
        }
    }
    return result;

}

void SparseToeplitz::operator*=(double c) 
{
    for (int i=0;i<Num_Diags;i++)
    {
        Vals[i]*=c;
    }
}

Matrix* SparseToeplitz::Kronecker(Matrix& B)
{
   
    BlockToeplitz* result = new BlockToeplitz(B.rows()*Num_Rows,B.cols()*Num_Cols,Num_Diags);
    
    for(int i=0;i<Num_Diags;i++)
    {
        result->Diags[i]=Diags[i];
        result->Vals[i] = B.Clone(Vals[i]);   
    }
    // avoid printing pointer with wrong format specifier (would be undefined behavior)
    // (If needed for debugging, print an informative message using std::cout.)
    return result;
}

Matrix* SparseToeplitz::Clone(double c)
{
   return new SparseToeplitz(*this, c);
}

