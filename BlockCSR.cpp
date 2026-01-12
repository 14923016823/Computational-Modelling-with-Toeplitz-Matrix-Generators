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
#pragma omp parallel for
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

double BlockCSR::operator()(int i, int j) const 
{ 
    if (i < 0 || i >= Num_Rows || j < 0 || j >= Num_Cols)
        throw std::invalid_argument("Index outside of matrix boundaries");

    int block_rows = Vals[0]->rows();
    int block_cols = Vals[0]->cols();
    
    int block_row = i / block_rows;
    int block_col = j / block_cols;
    int sub_i = i % block_rows;
    int sub_j = j % block_cols;

    for (int k = Rows[block_row]; k < Rows[block_row+1]; k++) //Go through the row which the value is in
    {
        if (Cols[k] == block_col)
        { 
            // Access the sub-matrix element
            return Vals[k]->operator()(sub_i, sub_j);
        }
    }
    return 0.0;
}

Matrix* BlockCSR::negativeTranspose()
{
    BlockCSR* negTrans = new BlockCSR(Num_Cols, Num_Rows, Num_Vals);
    negTrans->Rows[0] = 0;
    int n;
    #pragma omp parallel for private(n)
    for(int c=0;c<Num_Cols;c++) //Loop over columns of original matrix
    {
        n=0;
        for(int r=0;r<Num_Rows;r++) //Loop over rows of original matrix
        {
            for(int i=Rows[r];i<Rows[r+1];i++) //Loop over values
            {
                if (Cols[i]==c) //Only add a value if the current entry of Cols is the same as the current column
                {
                    n++;
                }
            }
        }
        negTrans->Rows[c+1] = n; //Set number of transposed row entries
    }

    for(int c = 0;c<Num_Cols;c++)
    {
        negTrans->Rows[c+1] += negTrans->Rows[c];
    }
    //All of the above is just making the Rows array, below is actually inserting values with correct column indices
    int m;
    #pragma omp parallel for private(m)
    for(int c=0;c<Num_Cols;c++) //Loop over columns of original matrix
    {
        m=0;
        for(int r=0;r<Num_Rows;r++) //Loop over rows of original matrix
        {
            for(int i=Rows[r];i<Rows[r+1];i++) //Loop over values in current row
            {
                if(Cols[i]==c && m<negTrans->Rows[c+1])
                {
                    negTrans->Vals[negTrans->Rows[c]+m] = Vals[i]->negativeTranspose();
                    negTrans->Cols[negTrans->Rows[c]+m] = r;
                    m++;
                }
            }
        }
    }
    return negTrans;
}

Matrix* BlockCSR::printFullMatrix()
{
    std::cout << "\nFull Dense Expansion (" << Num_Rows << "x" << Num_Cols << ")\n";
    int br = Vals[0]->rows();
    int bc = Vals[0]->cols();

    for (int k = 0; k < Num_Vals; k++) 
    {
        if (Vals[k]->rows() != br || Vals[k]->cols() != bc) 
        {
            throw std::logic_error("Inconsistent block sizes in BlockCSR");
        }
    }

    for (int i = 0; i < Num_Rows; i++) 
    {
        for (int j = 0; j < Num_Cols; j++) 
        {
            std::cout << std::setw(4) << this->operator()(i,j);
        }
        std::cout << "\n";
    }
    
    return nullptr;
}

Matrix* BlockCSR::Add(Matrix& other)
{
    return nullptr;
}