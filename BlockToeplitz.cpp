#include <vector>
#include <stdexcept>
#include <iostream>
//#include <tuple>
//#include "Vector.h"
#include <iomanip>

#include "BlockToeplitz.h"


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
    printf("delting blockToeplitz\n");
    for(int i=0;i<Num_Diags;i++)
    {
        delete Vals[i];
    }
    delete Vals;
    delete Diags;
}

// constructo

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
    printf("num diags=%d\n",Num_Diags);
    printf("vals[0],%d\n",Vals[0]);
    for(int i;i<Num_Diags;i++)
    {
        result->Diags[i]=Diags[i];
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


Matrix* BlockToeplitz::negativeTranspose()
{
    BlockToeplitz* negTrans = new BlockToeplitz(*this, -1.0);
    return negTrans;
}

double BlockToeplitz::operator()(int i, int j) const
{
    //determine which block we are in
    int block_rows = Vals[0]->rows();
    int block_cols = Vals[0]->cols();

    int block_row = i / block_rows;
    int block_col = j / block_cols;
    int sub_i = i % block_rows;
    int sub_j = j % block_cols;

    //find which diagonal this is
    int diag_index = block_col - block_row;
    for (int k = 0; k < Num_Diags; k++) {
        if (Diags[k] == diag_index) {
            // Access the sub-matrix element
            return (*Vals[k])(sub_i, sub_j);
        }
    }
    return 0.0; // element is zero if not on any stored diagonal
}

Matrix* BlockToeplitz::printFullMatrix()
{
    //print function for full dense expansion
    std::vector<std::vector<double>> M(Num_Rows, std::vector<double>(Num_Cols, 0.0));

    int block_rows = Vals[0]->rows();
    int block_cols = Vals[0]->cols();

    for (int k = 0; k < Num_Diags; k++) {
        for (int i = 0; i < Num_Rows; i++) {
            for (int j = 0; j < Num_Cols; j++) {
                int block_row = i / block_rows;
                int block_col = j / block_cols;
                int sub_i = i % block_rows;
                int sub_j = j % block_cols;

                if (block_col - block_row == Diags[k]) {
                    // Access the sub-matrix element
                    // Here we assume Vals[k] is a BlockToeplitz or similar with operator() defined
                    // You may need to adjust this part based on your actual MatrixPointer implementation
                    M[i][j] += static_cast<BlockToeplitz*>(Vals[k])->operator()(sub_i, sub_j);
                }
            }
        }
    }

    std::cout << "\nFull Dense Expansion (" << Num_Rows << "x" << Num_Cols << ")\n";

    for (int i = 0; i < Num_Rows; ++i) {
        for (int j = 0; j < Num_Cols; ++j)
            std::cout << std::setw(4) << std::left << M[i][j];
        std::cout << "\n";
    }

    return nullptr; // Adjust return type as needed
}