#include "BlockToeplitz.h"

#include <stdexcept>
#include <iostream>


// constructor
BlockToeplitz::BlockToeplitz(int nrows, int ncols, int ndiags)
{
    Num_Rows = nrows;
    Num_Cols = ncols;
    Num_Diags = ndiags;
    Diags = new int[Num_Diags];
    Vals = new MatrixPointer[Num_Diags];
    //Vals.resize(Num_Diags); 
}



BlockToeplitz::~BlockToeplitz()
{
    std::cout << "deleting BlockToeplitz\n";
    for(int i = 0; i < Num_Diags; ++i)
    {
        delete Vals[i];
    }
    delete[] Vals;
    delete[] Diags;
}

// constructor

//copy and scale constructor
BlockToeplitz::BlockToeplitz(BlockToeplitz& other,double c)
{

    Num_Rows=other.Num_Rows;
    Num_Cols=other.Num_Rows;
    Num_Diags=other.Num_Diags;
    Diags=new int[Num_Diags];
    Vals = new MatrixPointer[Num_Diags];
    for(int i=0;i<Num_Diags;i++)
    {
        Diags[i]=other.Diags[i];
        Vals[i]=other.Clone(c);
    }

}

Vectord BlockToeplitz::operator*(Vectord& vec)
{

    if (vec.len() != Num_Cols)
    {
        throw std::invalid_argument("Vector length and matrix columns don't match (block_toeplitz).");
    }

    int Num_Rows_SubMatrixes=Vals[0]->rows();
    int Num_Cols_SubMatrixes=Vals[0]->cols();

    
    int Num_blockrows = Num_Rows / Num_Rows_SubMatrixes;
    int Num_blockcols = Num_Cols / Num_Cols_SubMatrixes;

    Vectord result(Num_Rows);
    for (int blockrow = 0; blockrow < Num_blockrows; blockrow++)
    {
        //int row = blockrow *;
        Vectord subres(Num_Rows_SubMatrixes); // initialized to zeros
        // accumulate contributions from each diagonal/block
        for (int j = 0; j < Num_Diags; j++)
        {
            int blockcol = Diags[j] + blockrow;
            // ensure the whole sub-block fits in input vector
            if (blockcol >= 0 && blockcol<Num_blockcols)
            {
                Vectord subinput(Num_Cols_SubMatrixes);
                int col_start=blockcol*Num_Cols_SubMatrixes;
                for (int k = 0; k < Num_Cols_SubMatrixes; ++k)
                {
                    subinput.Vec[k] = vec.Vec[col_start + k]; 
                }
                Vectord subsubres = Vals[j]->operator*(subinput);
                subres.axpy(1.0, subsubres);
            }
        }
        // copy accumulated block into result
        int row_start=blockrow*Num_Rows_SubMatrixes;
        for (int k = 0; k < Num_Rows_SubMatrixes; ++k)
        {
            result.Vec[row_start + k] = subres.Vec[k];
        }
    }
    return result;
}

/*Vectord BlockToeplitz::operator*(Vectord& vec)
{
    if (vec.len() != Num_Cols)
    {
        throw std::invalid_argument("Vector length and matrix columns don't match (block_toeplitz).");
    }
    Vectord* result=new VectordNum_Rows;

}*/

void BlockToeplitz::operator*=(double c)
{

    for(int i=0;i<Num_Diags;i++)
    {
        (*Vals[i])*=c;
    }
        
}

Matrix* BlockToeplitz::Kronecker(Matrix& B)//if you add a new matrix at the bottom of the chain every Num_Rows needs to be changed
{

    BlockToeplitz* result = new BlockToeplitz(B.rows()*Num_Rows,B.cols()*Num_Cols,Num_Diags);
    std::cout << "num diags=" << Num_Diags << "\n";
    for(int i = 0; i < Num_Diags; ++i)
    {
        result->Diags[i] = Diags[i];
        result->Vals[i] = Vals[i]->Kronecker(B);
    }
    return result;
}

//Returns a cloned and scaled version of the input matrix
Matrix* BlockToeplitz::Clone(double c)
{
    BlockToeplitz* copy = new BlockToeplitz(*this, c);
    return copy;
}


