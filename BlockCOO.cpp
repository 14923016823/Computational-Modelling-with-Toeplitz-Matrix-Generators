#include "BlockCOO.h"


// constructor
BlockCOO::BlockCOO(int nrows, int ncols, int nvals)
{
    Num_Rows = nrows;
    Num_Cols = ncols;
    Num_Vals = nvals;
    Array = new mtuple[Num_Vals];
    //Vals.resize(Num_Diags); 
}

BlockCOO::~BlockCOO()
{
    printf("deleting blockCOO\n");
    delete Array;
}

// constructo

//copy and scale constructor
BlockCOO::BlockCOO(BlockCOO& other,double c)
{

    Num_Rows=other.Num_Rows;
    Num_Cols=other.Num_Cols;
    Num_Vals=other.Num_Vals;
    Array = new mtuple[Num_Vals];
    for(int i=0;i<Num_Vals;i++)
    {
        (*std::get<0>(Array[i]))*=c;
    }

}

Vectord BlockCOO::operator*(Vectord& vec)
{

    if (vec.len() != Num_Cols)
    {
        throw std::invalid_argument("Vector length and matrix columns don't match (block_toeplitz).");
    }

    int Num_Rows_SubMatrixes=std::get<0>(Array[0])->rows();
    int Num_Cols_SubMatrixes=std::get<0>(Array[0])->cols();

    
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
            int blockcol = std::get<2>(Array[j]) + blockrow;
            // ensure the whole sub-block fits in input vector
            if (blockcol >= 0 && blockcol<Num_blockcols)
            {
                Vectord subinput(Num_Cols_SubMatrixes);
                int col_start=blockcol*Num_Cols_SubMatrixes;
                for (int k = 0; k < Num_Cols_SubMatrixes; ++k)
                {
                    subinput.Vec[k] = vec.Vec[col_start + k]; 
                }
                Vectord subsubres = std::get<0>(Array[j])->operator*(subinput);
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

void BlockCOO::operator*=(double c)
{

    for(int i=0;i<Num_Vals;i++)
    {
        (*std::get<0>(Array[i]))*=c;
    }
        
}

Matrix* BlockCOO::Kronecker(Matrix& B)//if you add a new matrix at the bottom of the chain every Num_Rows needs to be changed
{

    BlockCOO* result = new BlockCOO(B.rows()*Num_Rows,B.cols()*Num_Cols,Num_Vals);
    printf("num diags=%d\n",Num_Vals);
    printf("vals[0],%d\n",std::get<0>(Array[0]));
    for(int i;i<Num_Vals;i++)
    {
        std::get<0>(result->Array[i])=std::get<0>(Array[i])->Kronecker(B);
        std::get<1>(result->Array[i])=std::get<1>(Array[i]);
        std::get<2>(result->Array[i])=std::get<2>(Array[i]);
    }
    return result;
}

//Returns a cloned and scaled version of the input matrix
Matrix* BlockCOO::Clone(double c)
{
    BlockCOO* copy = new BlockCOO(*this, c);
    return copy;
}


