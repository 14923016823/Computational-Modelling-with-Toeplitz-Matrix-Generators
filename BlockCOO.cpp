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

//copy and scale constructor
BlockCOO::BlockCOO(BlockCOO& other,double c)
{

    Num_Rows=other.Num_Rows;
    Num_Cols=other.Num_Cols;
    Num_Vals=other.Num_Vals;
    Array = new mtuple[Num_Vals];
    for(int i=0;i<Num_Vals;i++)
    {
        (*std::get<2>(Array[i]))*=c;
    }

}

Vectord BlockCOO::operator*(Vectord& vec)
{

    if (vec.len() != Num_Cols)
    {
        throw std::invalid_argument("Vector length and matrix columns don't match (block_toeplitz).");
    }

    int Num_Rows_SubMatrixes=std::get<2>(Array[0])->rows();
    int Num_Cols_SubMatrixes=std::get<2>(Array[0])->cols();

    
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
            int blockcol = std::get<1>(Array[j]) + blockrow;
            // ensure the whole sub-block fits in input vector
            if (blockcol >= 0 && blockcol<Num_blockcols)
            {
                Vectord subinput(Num_Cols_SubMatrixes);
                int col_start=blockcol*Num_Cols_SubMatrixes;
                for (int k = 0; k < Num_Cols_SubMatrixes; ++k)
                {
                    subinput.Vec[k] = vec.Vec[col_start + k]; 
                }
                Vectord subsubres = std::get<2>(Array[j])->operator*(subinput);
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
        (*std::get<2>(Array[i]))*=c;
    }
        
}

Matrix* BlockCOO::Kronecker(Matrix& B)//if you add a new matrix at the bottom of the chain every Num_Rows needs to be changed
{

    BlockCOO* result = new BlockCOO(B.rows()*Num_Rows,B.cols()*Num_Cols,Num_Vals);
    printf("num diags=%d\n",Num_Vals);
    printf("vals[0],%d\n",std::get<0>(Array[0]));
    for(int i;i<Num_Vals;i++)
    {
        std::get<0>(result->Array[i])=std::get<0>(Array[i]);
        std::get<1>(result->Array[i])=std::get<1>(Array[i]);
        std::get<2>(result->Array[i])=std::get<2>(Array[i])->Kronecker(B);
    }
    return result;
}

//Returns a cloned and scaled version of the input matrix
Matrix* BlockCOO::Clone(double c)
{
    BlockCOO* copy = new BlockCOO(*this, c);
    return copy;
}

double BlockCOO::operator()(int i, int j) const 
{ 
    if (i < 0 || i >= Num_Rows || j < 0 || j >= Num_Cols)
        throw std::invalid_argument("Index outside of matrix boundaries");
    
    //determine which block we are in
    int block_rows = std::get<2>(Array[0])->rows();
    int block_cols = std::get<2>(Array[0])->cols();

    int block_row = i / block_rows;
    int block_col = j / block_cols;
    int sub_i = i % block_rows;
    int sub_j = j % block_cols;

    //find which element this is
    for (int k = 0; k < Num_Vals; k++) {
        if (std::get<0>(Array[k]) == block_row && std::get<1>(Array[k]) == block_col) {
            // Access the sub-matrix element
            return (*std::get<2>(Array[k]))(sub_i, sub_j);
        }
    }
    return 0.0; // element is zero if not on any stored coordinate 
}

Matrix* BlockCOO::negativeTranspose()
{
    // Create a new BlockCOO that represents the negative transpose
    /*BlockCOO* negTrans = new BlockCOO(Num_Cols, Num_Rows, Num_Vals);
    int n = 0;
    for(int c=0;c<Num_Cols;c++)
    {
        for(int r=0;r<Num_Rows;r++)
        {
            for(int i=0;i<Num_Vals;i++)
            {
                if(r == std::get<0>(Array[i]) && c==std::get<1>(Array[i]))
                {    
                    std::get<0>(negTrans->Array[n])=std::get<1>(Array[i]);
                    std::get<1>(negTrans->Array[n])=std::get<0>(Array[i]);
                    std::get<2>(negTrans->Array[n])=std::get<2>(Array[i])->negativeTranspose();
                    n++;
                }
            }
        }
    }
    return negTrans;*/

    BlockCOO* negTrans = new BlockCOO(Num_Cols, Num_Rows, Num_Vals);
    int Rows[Num_Cols+1]; //Array to count entries per row in transposed matrix
    Rows[0] = 0;
    for(int c=0;c<Num_Cols;c++) //Go through columns of original matrix
    {
        Rows[c+1]=0;
        for(int i=0;i<Num_Vals;i++) //Go through values
        {
            if(c==std::get<1>(Array[i])) //Check whether current value is in current column
            {
                Rows[c+1]++;
            }
        }
    }
    for(int c = 0;c<Num_Cols+1;c++)
    {
        Rows[c+1] += Rows[c];
    }

    int n;
#pragma omp parallel for private(n)
    for(int c=0;c<Num_Cols;c++)
    {
        n=0;
        for(int r=0;r<Num_Rows;r++)
        {
            for(int i=Rows[r];i<Rows[r+1];i++) //Loop over values in current row
            {
                if(c==std::get<1>(Array[i]) && n<Rows[c+1])
                {
                    std::get<0>(negTrans->Array[Rows[c]+n])=std::get<1>(Array[i]);
                    std::get<1>(negTrans->Array[Rows[c]+n])=std::get<0>(Array[i]);
                    std::get<2>(negTrans->Array[Rows[c]+n]) = std::get<2>(Array[i])->negativeTranspose();
                    n++;
                }
            }
        }
    }
    return negTrans;
}

Matrix* BlockCOO::printFullMatrix()
{
    std::cout << "\nFull Dense Expansion (" << Num_Rows << "x" << Num_Cols << ")\n";
    int br = std::get<2>(Array[0])->rows();
    int bc = std::get<2>(Array[0])->cols();

    for (int k = 0; k < Num_Vals; k++) {
        if (std::get<2>(Array[k])->rows() != br || std::get<2>(Array[k])->cols() != bc) 
        {
            throw std::logic_error("Inconsistent block sizes in BlockCOO");
        }
    }

    for (int i = 0; i < Num_Rows; i++) 
    {
        for (int j = 0; j < Num_Cols; j++) 
        {
            for (int k = 0; k < Num_Vals; k++) 
            {
                std::cout << std::setw(4) << this->operator()(i,j);
            }
        }
        std::cout << "\n";
    }
    
    return nullptr;
}