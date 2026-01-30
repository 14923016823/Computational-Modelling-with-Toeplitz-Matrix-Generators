#include "BlockToeplitz.h"


// constructor
BlockToeplitz::BlockToeplitz(int nrows, int ncols, int ndiags)
{
    Num_Rows = nrows;
    Num_Cols = ncols;
    Num_Diags = ndiags;
    Diags = new int[Num_Diags];
    Vals = new MatrixPointer[Num_Diags];
    for(int i=0; i<Num_Diags; i++)
        Vals[i] = nullptr;
}

BlockToeplitz::~BlockToeplitz()
{
    //printf("deleting blockToeplitz\n");
    for(int i=0;i<Num_Diags;i++)
    {
        delete Vals[i];
    }
    delete[] Vals;
    delete[] Diags;
}

// constructo

//copy and scale constructor
BlockToeplitz::BlockToeplitz(BlockToeplitz& other,double c)
{

    Num_Rows=other.Num_Rows;
    Num_Cols=other.Num_Cols;
    Num_Diags=other.Num_Diags;
    Diags=new int[Num_Diags];
    Vals = new MatrixPointer[Num_Diags];
    for(int i=0;i<Num_Diags;i++)
    {
        Diags[i]=other.Diags[i];
        Vals[i]=other.Vals[i]->Clone(c);
    }

}

void BlockToeplitz::MatMulAdd(const Vectord& x,const int x_start,Vectord& result,const int result_start)
{
    fflush(stdout);
    int br = Vals[0]->rows();  // block row size
    int bc = Vals[0]->cols();  // block col size

    int numBlockRows = Num_Rows / br;
    int numBlockCols = Num_Cols / bc;

    for (int k = 0; k < Num_Diags; ++k) 
    {
        int d = Diags[k];       // block-local diagonal
        Matrix* B = Vals[k];    // the block
        //std::cout << d << std::endl;
        //std::cout << Vals[k]->rows() << " " << Vals[k]->cols() << "\n";
    
        if (d >= 0) 
        {
            int numBlocks=std::min(numBlockRows,numBlockCols-d);
            for(int i=0;i<numBlocks;i++)
            {
                B->MatMulAdd(x,x_start+(d+i)*bc,result,result_start+i*br);
            }
        } 
        else 
        {
            int numBlocks=std::min(numBlockRows+d,numBlockCols);
            for(int i=0;i<numBlocks;i++)
            {
                B->MatMulAdd(x,x_start+i*bc,result,result_start+(i-d)*br);
            }
        }
        //std::cout << d << std::endl;
        //std::cout << Vals[k]->rows() << " " << Vals[k]->cols() << "\n";
    }
}

void BlockToeplitz::MatMul(const Vectord& x,Vectord& result)
{
    if(x.len()!=Num_Cols||result.len()!=Num_Rows)
    {
        throw std::invalid_argument("Vector and Matrix size dont match matmul (BT)");
    }
    for(int i=0;i<result.len();i++)
    {
        result[i]=0.0;
    }
    MatMulAdd(x,0,result,0);
}

Vectord BlockToeplitz::operator*(Vectord& vec)//try not to use this function
{
    Vectord result = Vectord(Num_Rows);
    MatMul(vec,result);
    return result;
}
/*
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
    //int i;
    #pragma omp parallel for //private(i)
    for (int blockrow = 0; blockrow < Num_blockrows; blockrow++)
    {
        //int i = omp_get_thread_num();
    
        //printf("Hello World... from thread = %d\n", i);
        //int row = blockrow *;
        Vectord subres(Num_Rows_SubMatrixes); // initialized to zeros
        // accumulate contributions from each diagonal/block
        for (int j = 0; j < Num_Diags; j++)
        {
            //printf("i=%d, j=%d\n", i,j);
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
                Vals[j]->printFullMatrix();
                std::cout << subsubres.len() << std::endl;
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
*/
void BlockToeplitz::operator*=(double c)
{
    int i;
    //#pragma omp parallel for private(i)
    for(i=0;i<Num_Diags;i++)
    {
        Vals[i]->Clone(c);
    }
        
}

Matrix* BlockToeplitz::Kronecker(Matrix& B)//if you add a new matrix at the bottom of the chain every Num_Rows needs to be changed
{
    BlockToeplitz* result = new BlockToeplitz(B.rows()*Num_Rows,B.cols()*Num_Cols,Num_Diags);
    //printf("num diags=%d\n",Num_Diags);
    //printf("vals[0],%d\n",Vals[0]);
    int i;
//#pragma omp parallel for private(i)
    for(i=0;i<Num_Diags;i++)
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
    // Create a new BlockToeplitz that represents the negative transpose
    BlockToeplitz* negTrans = new BlockToeplitz(Num_Cols, Num_Rows, Num_Diags);
    // Reverse and negate diagonals, and negative-transpose each sub-block
    for (int d = 0; d < Num_Diags; ++d) {
        int src = Num_Diags - 1 - d;
        negTrans->Diags[d] = -Diags[src];
        // call negativeTranspose on the sub-block (returns Matrix*)
        negTrans->Vals[d] = Vals[src]->negativeTranspose();
    }
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
            return Vals[k]->operator()(sub_i, sub_j);
        }
    }
    return 0.0; // element is zero if not on any stored diagonal
}

void BlockToeplitz::printFullMatrix() 
{
    std::cout << "\nFull Dense Expansion (" << Num_Rows << "x" << Num_Cols << ")\n";
    int br = Vals[0]->rows();
    int bc = Vals[0]->cols();

    for (int k = 0; k < Num_Diags; k++) {
        if (Vals[k]->rows() != br || Vals[k]->cols() != bc) {
            throw std::logic_error("Inconsistent block sizes in BlockToeplitz");
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
}