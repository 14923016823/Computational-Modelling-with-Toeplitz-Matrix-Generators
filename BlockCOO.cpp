#include "BlockCOO.h"


// constructor
BlockCOO::BlockCOO(int nrows, int ncols, int nvals)
{
    Num_Rows = nrows;
    Num_Cols = ncols;
    Num_Vals = nvals;
    Array = new mtuple[Num_Vals];
    for(int i=0; i<Num_Vals; i++)
        Array[i] = std::make_tuple(0,0,nullptr);
    //Vals.resize(Num_Diags); 
}

BlockCOO::~BlockCOO()
{
    printf("deleting blockCOO\n");
    delete[] Array;
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
        std::get<2>(Array[i])=std::get<2>(other.Array[i])->Clone(c);
    }

}

Vectord BlockCOO::operator*(Vectord& vec)
{
    if (vec.len() != Num_Cols)
    {
        throw std::invalid_argument("Vector length and matrix columns don't match (block_COO).");
    }

    int Num_Rows_SubMatrixes=std::get<2>(Array[0])->rows();
    int Num_Cols_SubMatrixes=std::get<2>(Array[0])->cols();

    
    //int Num_blockrows = Num_Rows / Num_Rows_SubMatrixes;
    //int Num_blockcols = Num_Cols / Num_Cols_SubMatrixes;

    Vectord result(Num_Rows);

        //int row = blockrow *;
        int j=0;
     // initialized to zeros
        // accumulate contributions from each diagonal/block
//#pragma omp parallel for //private(subres)

    for (j = 0; j < Num_Vals; j+=1)
    {
        //printf("j=%d\n",j);
        //Vectord subres(Num_Rows_SubMatrixes);
        
        int blockcol = std::get<1>(Array[j]);
        int blockrow = std::get<0>(Array[j]);
        // ensure the whole sub-block fits in input vector
        
        Vectord subinput(Num_Cols_SubMatrixes);
        Vectord subsubres(blockcol);
        int col_start=blockcol*Num_Cols_SubMatrixes;
        for (int k = 0; k < Num_Cols_SubMatrixes; ++k)
        {
            subinput.Vec[k] = vec.Vec[col_start + k]; 
        }
        //printf("a\n");
        
        
        subsubres = std::get<2>(Array[j])->operator*(subinput);
        //printf("hi\n");
        //subres.Sum(subsubres);
        
        int row_start=blockrow*Num_Rows_SubMatrixes;
        for (int k = 0; k < Num_Rows_SubMatrixes; ++k)
        {
            //#pragma omp atomic update
            result.Vec[row_start + k] += subsubres.Vec[k];
        }
    }
        // copy accumulated block into result
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
    for(int i=0;i<Num_Vals;i++)
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
            return std::get<2>(Array[k])->operator()(sub_i, sub_j);
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
    //int block_rows = std::get<2>(Array[0])->rows();
    int block_cols = std::get<2>(Array[0])->cols();

    BlockCOO* negTrans = new BlockCOO(Num_Cols, Num_Rows, Num_Vals);
    int Rows[Num_Cols/block_cols+1]; //Array to count entries per row in transposed matrix
    for(int j=0;j<Num_Cols/block_cols+1;j++)
        Rows[j] = 0;
//#pragma omp parallel for
    for(int c=0;c<Num_Cols/block_cols;c++) //Go through columns of original matrix
    {
        for(int i=0;i<Num_Vals;i++) //Go through values
        {
            if(c==std::get<1>(Array[i])) //Check whether current value is in current column
            {
                Rows[c+1]++;
                //std::cout << c << '\n';
            }
        }
    }
    for(int c = 0;c<Num_Cols/block_cols;c++)
    {
        Rows[c+1] += Rows[c];
        //std::cout << Rows[c+1] << '\n';
    }

    int n;
    //#pragma omp parallel for private(n)
    for(int c=0;c<Num_Cols/block_cols;c++) //Cols of original matrix, so rows of transposed matrix
    {
        n=Rows[c];
        //printf("%d\n",c);
        //for(int r=0;r<Num_Rows/block_rows;r++) //Rows of original matrix, so cols of transposed matrix
        //{   
            //printf("r = %d\n",r);
            for(int i=0;i<Num_Vals;i++) //Loop over values in current row
            {
                //printf("c=%d, i=%d\n",c, i);
                //printf("%d\n",n);
                if(c==std::get<1>(Array[i]) && n<Rows[c+1])
                {
                    //printf("yes, c=%d, i=%d, n=%d\n",c, i,n);
                    //printf("Rows[%d+1] = %d, n = %d\n",c, Rows[c+1], n);
                    //printf("block (%d, %d)\n", std::get<1>(Array[i]), std::get<0>(Array[i]));
                    std::get<0>(negTrans->Array[n])=std::get<1>(Array[i]);
                    std::get<1>(negTrans->Array[n])=std::get<0>(Array[i]);
                    std::get<2>(negTrans->Array[n]) = std::get<2>(Array[i])->negativeTranspose();
                    //std::get<2>(Array[i])->printFullMatrix();
                    //printf("finished c=%d, i=%d, n=%d\n",c, i,n);
                    n++;
                }
            }
        //}
    }
    return negTrans;
}

void BlockCOO::printFullMatrix()
{
    std::cout << "\nFull Dense Expansion (" << Num_Rows << "x" << Num_Cols << ")\n";
    //std::cout << std::get<0>(Array[0]) << '\n';
    //std::cout << std::get<1>(Array[0]) << '\n';
    int br = std::get<2>(Array[0])->rows();
    int bc = std::get<2>(Array[0])->cols();

    for (int k = 0; k < Num_Vals; k++) 
    {
        //std::cout << std::get<2>(Array[k])->rows() << ", " << std::get<2>(Array[k])->cols() << " != " << br << ", " << bc << '\n';
        //std::get<2>(Array[k])->printFullMatrix();
        if (std::get<2>(Array[k])->rows() != br || std::get<2>(Array[k])->cols() != bc) 
        {
            throw std::logic_error("Inconsistent block sizes in BlockCOO");
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

void BlockCOO::MatMulAdd(const Vectord& x,const int x_start,Vectord& result,const int result_start)
{
    throw std::invalid_argument("Not implemented");
}

void BlockCOO::MatMul(const Vectord& x,Vectord& result)
{
    throw std::invalid_argument("Not implemented");
}