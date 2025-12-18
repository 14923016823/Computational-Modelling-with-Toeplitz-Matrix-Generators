#include "SparseToeplitz.h"
#include "BlockToeplitz.h"
#include "VectorD.h"
    
SparseToeplitz::SparseToeplitz(int num_rows, int num_cols, int length)
{
    if (length > num_rows + num_cols - 1)
    {
        throw std::invalid_argument("There are too many diagonals for Toeplitz matrix of this size");
    }
    Num_Rows = num_rows;
    Num_Cols = num_cols;
    Num_Diags = length;
    Diags = new int[Num_Diags];
    Vals = new double[Num_Diags];
}
SparseToeplitz::SparseToeplitz(int num_rows, int num_cols, int length,int* diags,double* vals)
{
    if (length > num_rows + num_cols - 1)
    {
        throw std::invalid_argument("There are too many diagonals for Toeplitz matrix of this size");
    }
    Num_Rows = num_rows;
    Num_Cols = num_cols;
    Num_Diags = length;
    Diags = diags;
    Vals = vals;
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

Vectord SparseToeplitz::operator*(Vectord& vec)
{
    int len = vec.Length;
    if (len != Num_Cols) 
    {
        throw std::invalid_argument("Vector and Matrix size dont match");
    }
    Vectord result = Vectord(len);
    #pragma omp parallel for
    for (int i = 0;i < Num_Rows;i++) 
    {
        for (int j = 0;j < Num_Diags;j++) 
        {
            if (Diags[j] + i >= 0 && Diags[j] + i < len)
            {
                result.Vec[i] += vec.Vec[Diags[j] + i] * Vals[j];
            }
        }
    }
    return result;
}

void SparseToeplitz::print()
{
    std::cout << "Values: [";
    for(int i = 0;i<Num_Diags;i++)
    {
        std::cout << Vals[i];
        if(i<Num_Diags-1)
        std::cout << ", ";
    }
    std::cout << ']' << std::endl;

    std::cout << "Diagonals: [";
    for(int i = 0;i<Num_Diags;i++)
    {
        std::cout << Diags[i];
        if(i<Num_Diags-1)
        std::cout << ", ";
    }
    std::cout << ']' << std::endl;
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

