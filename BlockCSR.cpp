#include "BlockCSR.h"


// constructor
BlockCSR::BlockCSR(int nrows, int ncols, int nvals)
{
    Num_Rows = nrows;
    Num_Cols = ncols;
    Num_Vals = nvals;
    Rows = new int[Num_Rows];
    Cols = new int[Num_Cols];
    Vals = new MatrixPointer[Num_Vals];
    //Vals.resize(Num_Diags); 
}

BlockCSR::~BlockCSR()
{
    printf("deleting blockCSR\n");
    for(int i=0;i<Num_Vals;i++)
    {
        delete Vals[i];
    }
    delete Vals;
    delete Rows;
    delete Cols;
}

// constructo

//copy and scale constructor
BlockCSR::BlockCSR(BlockCSR& other,double c)
{

    Num_Rows=other.Num_Rows;
    Num_Cols=other.Num_Cols;
    Num_Vals=other.Num_Vals;
    Cols = new int[Num_Vals];
    Rows = new int[Num_Rows+1];
    Vals = new MatrixPointer[Num_Vals];
    for(int i=0;i<Num_Vals;i++)
    {
        Cols[i]=other.Cols[i];
        Vals[i]=other.Clone(c);
    }
    for(int i=0;i<Num_Rows+1;i++)
    {
        Rows[i]=other.Rows[i];
    }

}

Vectord BlockCSR::operator*(Vectord& vec)
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
        for (int j = 0; j < Num_Vals; j++)
        {
            int blockcol = Cols[j] + blockrow;
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
                subres.Sum(subsubres);
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

void BlockCSR::operator*=(double c)
{

    for(int i=0;i<Num_Vals;i++)
    {
        (*Vals[i])*=c;
    }
        
}

Matrix* BlockCSR::Kronecker(Matrix& B)//if you add a new matrix at the bottom of the chain every Num_Rows needs to be changed
{

    BlockCSR* result = new BlockCSR(B.rows()*Num_Rows,B.cols()*Num_Cols,Num_Vals);
    printf("num diags=%d\n",Num_Vals);
    printf("vals[0],%d\n",Vals[0]);
    for(int i;i<Num_Vals;i++)
    {
        result->Cols[i]=Cols[i];
        result->Vals[i] = Vals[i]->Kronecker(B);
        
    }
    for(int i=0;i<Num_Rows+1;i++)
    {
        result->Rows[i]=Rows[i];
    }
    return result;
}

//Returns a cloned and scaled version of the input matrix
Matrix* BlockCSR::Clone(double c)
{
    BlockCSR* copy = new BlockCSR(*this, c);
    return copy;
}


