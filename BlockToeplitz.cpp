#include "BlockToeplitz.h"


// constructor
BlockToeplitz::BlockToeplitz(int nrows, int ncols, int ndiags)
{
    Num_Rows = nrows;
    Num_Cols = ncols;
    Num_Diags = ndiags;
    Diags = new int[Num_Diags];
    Vals = new MatrixPointer[Num_Diags];
}

BlockToeplitz::~BlockToeplitz()
{
    printf("deleting blockToeplitz\n");
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
    int i;
    #pragma omp parallel for private(i)
    for(i=0;i<Num_Diags;i++)
    {
        (*Vals[i])*=c;
    }
        
}

Matrix* BlockToeplitz::Kronecker(Matrix& B)//if you add a new matrix at the bottom of the chain every Num_Rows needs to be changed
{

    BlockToeplitz* result = new BlockToeplitz(B.rows()*Num_Rows,B.cols()*Num_Cols,Num_Diags);
    printf("num diags=%d\n",Num_Diags);
    printf("vals[0],%d\n",Vals[0]);
    int i;
#pragma omp parallel for private(i)
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
            return (*Vals[k])(sub_i, sub_j);
        }
    }
    return 0.0; // element is zero if not on any stored diagonal
}

Matrix* BlockToeplitz::printFullMatrix() {
    std::vector<std::vector<double>> M(
        Num_Rows, std::vector<double>(Num_Cols, 0.0));

    int br = Vals[0]->rows();
    int bc = Vals[0]->cols();

    for (int k = 0; k < Num_Diags; k++) {
        if (Vals[k]->rows() != br || Vals[k]->cols() != bc) {
            throw std::logic_error("Inconsistent block sizes in BlockToeplitz");
        }
    }

    int num_block_rows = Num_Rows / br;
    int num_block_cols = Num_Cols / bc;

    for (int block_row = 0; block_row < num_block_rows; block_row++) {
        for (int block_col = 0; block_col < num_block_cols; block_col++) {

            int diag = block_col - block_row;

            for (int k = 0; k < Num_Diags; k++) {
                if (Diags[k] == diag) {

                    for (int i = 0; i < br; i++) {
                        for (int j = 0; j < bc; j++) {

                            int I = block_row * br + i;
                            int J = block_col * bc + j;

                            M[I][J] += (*Vals[k])(i, j);
                        }
                    }
                }
            }
        }
    }

    std::cout << "\nFull Dense Expansion (" << Num_Rows << "x" << Num_Cols << ")\n";
    for (int i = 0; i < Num_Rows; ++i) {
        for (int j = 0; j < Num_Cols; ++j)
            std::cout << std::setw(4) << M[i][j];
        std::cout << "\n";
    }

    return nullptr;
}

Matrix* BlockToeplitz::Add(Matrix& other)
{
    // If other is a BlockToeplitz, handle directly.
    BlockToeplitz* o = dynamic_cast<BlockToeplitz*>(&other);
    if (!o) {
        // If other is a SparseToeplitz, lift it to a BlockToeplitz (identity kronecker)
        SparseToeplitz* s = dynamic_cast<SparseToeplitz*>(&other);
        if (s) {
            // create a temporary BlockToeplitz that represents identity kron with s
            BlockToeplitz* tmp = new BlockToeplitz(Num_Rows, Num_Cols, 1);
            tmp->Diags[0] = 0;
            tmp->Vals[0] = s->Clone(1.0);
            Matrix* result = this->Add(*tmp);
            // cleanup temporary
            delete tmp;
            return result;
        }
        throw std::invalid_argument("BlockToeplitz::Add expects BlockToeplitz or SparseToeplitz");
    }

    if (Num_Rows != o->Num_Rows || Num_Cols != o->Num_Cols)
        throw std::invalid_argument("BlockToeplitz::Add dimension mismatch");

    // Build union of diagonals (both Diags arrays are assumed sorted)
    std::vector<int> diagUnion;
    int ia = 0, ib = 0;
    while (ia < Num_Diags || ib < o->Num_Diags) {
        if (ia < Num_Diags && (ib == o->Num_Diags || Diags[ia] < o->Diags[ib])) {
            diagUnion.push_back(Diags[ia++]);
        } else if (ib < o->Num_Diags && (ia == Num_Diags || o->Diags[ib] < Diags[ia])) {
            diagUnion.push_back(o->Diags[ib++]);
        } else {
            diagUnion.push_back(Diags[ia]); ++ia; ++ib;
        }
    }

    int n = (int)diagUnion.size();
    BlockToeplitz* R = new BlockToeplitz(Num_Rows, Num_Cols, n);
    for (int k = 0; k < n; ++k) {
        R->Diags[k] = diagUnion[k];

        // find indices in this and other
        int idxA = -1, idxB = -1;
        for (int i = 0; i < Num_Diags; ++i) if (Diags[i] == diagUnion[k]) { idxA = i; break; }
        for (int i = 0; i < o->Num_Diags; ++i) if (o->Diags[i] == diagUnion[k]) { idxB = i; break; }

        if (idxA >= 0 && idxB >= 0) {
            // both present: add sub-blocks
            Matrix* a = Vals[idxA];
            Matrix* b = o->Vals[idxB];
            BlockToeplitz* ablock = dynamic_cast<BlockToeplitz*>(a);
            BlockToeplitz* bblock = dynamic_cast<BlockToeplitz*>(b);
            if (ablock && bblock) {
                R->Vals[k] = ablock->Add(*bblock);
                continue;
            }
            SparseToeplitz* asparse = dynamic_cast<SparseToeplitz*>(a);
            SparseToeplitz* bsparse = dynamic_cast<SparseToeplitz*>(b);
            if (asparse && bsparse) {
                R->Vals[k] = asparse->Add(*bsparse);
                continue;
            }
            throw std::logic_error("Unsupported sub-block types in BlockToeplitz::Add");
        }

        // only one present: clone the present sub-block
        if (idxA >= 0) {
            R->Vals[k] = Vals[idxA]->Clone(1.0);
        } else if (idxB >= 0) {
            R->Vals[k] = o->Vals[idxB]->Clone(1.0);
        }
    }

    return R;
}