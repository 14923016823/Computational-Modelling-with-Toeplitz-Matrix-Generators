#include <vector>
#include <stdexcept>
#include <iostream>
//#include <tuple>
//#include "Vector.h"
#include <iomanip>

#include "SparseToeplitz.h"
#include "BlockToeplitz.h"


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



SparseToeplitz::~SparseToeplitz()
{
    printf("delting blockToeplitz\n");
    delete Vals;
    delete Diags;
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
    
    printf("vals[i],%f\n",result->Vals);
    return result;
}

Matrix* SparseToeplitz::Clone(double c)
{
   return new SparseToeplitz(*this, c);
}

Matrix* SparseToeplitz::negativeTranspose()
{
    SparseToeplitz* negTrans = new SparseToeplitz(*this, -1.0);
    return negTrans;
}

Matrix* SparseToeplitz::printFullMatrix()
{
    std::vector<std::vector<double>> M(Num_Rows, std::vector<double>(Num_Cols, 0.0));
    for (int d = 0; d < Num_Diags; d++) {
        for (int i = 0; i < Num_Rows; i++) {
            int j = i + Diags[d];
            if (j >= 0 && j < Num_Cols)
                M[i][j] = Vals[d]; // stores diagonal value
            
            // if symmetric Toeplitz, uncomment this:
            // if (Diags[d] > 0 && i - Diags[d] >= 0)
            //    M[i][i - Diags[d]] = ValsD[d];
        }
    }

    std::cout << "\nFull Dense Expansion (" << Num_Rows << "x" << Num_Cols << ")\n";

    for (int i = 0; i < Num_Rows; ++i) {
        for (int j = 0; j < Num_Cols; ++j)
            std::cout << std::setw(4) << std::left << M[i][j];
        std::cout << "\n";
    }
    return nullptr;
}