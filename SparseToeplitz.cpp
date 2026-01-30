#include "SparseToeplitz.h"
    
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



void SparseToeplitz::MatMulAdd(const Vectord& x, const int x_start,Vectord& result,const int result_start)
{
    #pragma omp parallel for
    for (int j = 0; j < Num_Diags; ++j) 
    {
        int d = Diags[j];
        double v = Vals[j];
        if(d>=0)
        {
            int diag_length=std::min(Num_Rows, Num_Cols-d);
            for(int i=0;i<diag_length;i++)
            {
                #pragma omp atomic update
                result[i+result_start]+=v*x[i+d+x_start];
            }
        }
        else
        {
            int diag_length=std::min(Num_Cols, Num_Rows+d);
            for(int i=0;i<diag_length;i++)
            {
                #pragma omp atomic update
                result[i-d+result_start]+=v*x[i+x_start];
            }
        }
    }
}

void SparseToeplitz::MatMul(const Vectord& x,Vectord& result)
{
    if(x.len()!=Num_Cols||result.len()!=Num_Rows)
    {
        throw std::invalid_argument("Vector and Matrix size dont match matmul (ST)");
    }
    for(int i=0;i<result.len();i++)
    {
        result[i]=0.0;
    }
    MatMulAdd(x,0,result,0);
}

Vectord SparseToeplitz::operator*(Vectord& vec)//try not to use this function
{
    Vectord result = Vectord(Num_Rows);
    MatMul(vec,result);
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
    return result;
}

double SparseToeplitz::operator()(int i, int j) const
{
    if (i < 0 || i >= Num_Rows || j < 0 || j >= Num_Cols)
        return 0.0;
    int diag = j - i;
    for (int d = 0; d < Num_Diags; ++d) {
        if (Diags[d] == diag)
            return Vals[d];
    }
    return 0.0;
}

Matrix* SparseToeplitz::Clone(double c)
{
   return new SparseToeplitz(*this, c);
}

Matrix* SparseToeplitz::negativeTranspose()
{
    // The negative transpose of a (Num_Rows x Num_Cols) Toeplitz matrix
    // is a (Num_Cols x Num_Rows) matrix where each diagonal k becomes -k
    // and the order of stored diagonals must be reversed to remain sorted ascending.
    SparseToeplitz* negTrans = new SparseToeplitz(Num_Cols, Num_Rows, Num_Diags);
    for (int d = 0; d < Num_Diags; ++d) {
        // reverse order and negate diagonal offsets and values
        negTrans->Diags[d] = -Diags[Num_Diags - 1 - d];
        negTrans->Vals[d] = -Vals[Num_Diags - 1 - d];
    }
    return negTrans;
}

void SparseToeplitz::printFullMatrix()
{
    std::cout << "\nFull Dense Expansion (" << Num_Rows << "x" << Num_Cols << ")\n";

    for (int i = 0; i < Num_Rows; i++) 
    {
        for (int j = 0; j < Num_Cols; j++) 
        {
            std::cout << std::setw(4) << this->operator()(i,j);
        }
        std::cout << "\n";
    }
}