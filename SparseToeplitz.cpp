#include "SparseToeplitz.h"
#include "BlockToeplitz.h"
#include "VectorD.h"
    
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

Vectord SparseToeplitz::operator*(Vectord& vec)
{
    int len = vec.Length;
    if (len != Num_Cols) 
    {
        throw std::invalid_argument("Vector and Matrix size dont match");
    }
    Vectord result = Vectord(Num_Rows);
    int q;
    int j;
    int i;
    #pragma omp parallel for private(q,j,i)
    for (i = 0;i < Num_Rows;i++) //This always completes the j-loop before moving on to the next i, so seems besically sequential
    {
        for (j = 0;j < Num_Diags;j++) 
        {
            if (Diags[j] + i >= 0 && Diags[j] + i < len)
            {
                result.Vec[i] += vec.Vec[Diags[j] + i] * Vals[j];
            }
        }
    }
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
    
    //printf("vals[i],%f\n",result->Vals);
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

Matrix* SparseToeplitz::printFullMatrix()
{
    std::vector<std::vector<double>> M(Num_Rows, std::vector<double>(Num_Cols, 0.0));
    for (int d = 0; d < Num_Diags; d++) {
        for (int i = 0; i < Num_Rows; i++) {
            int j = i + Diags[d];
            if (j >= 0 && j < Num_Cols)
                M[i][j] = Vals[d]; // stores diagonal value
            
            // if symmetric Toeplitz, uncomment this:
            // if (Diags[d] > 0 && i - Diags[d] >= 0)
            //    M[i][i - Diags[d]] = ValsD[d];
        }
    }

    std::cout << "\nFull Dense Expansion (" << Num_Rows << "x" << Num_Cols << ")\n";

    for (int i = 0; i < Num_Rows; ++i) {
        for (int j = 0; j < Num_Cols; ++j)
            std::cout << std::setw(4) << std::left << M[i][j];
        std::cout << "\n";
    }
    return nullptr;
}

Matrix* SparseToeplitz::Add(Matrix& other)
{
    SparseToeplitz* o = dynamic_cast<SparseToeplitz*>(&other);
    if (!o)
        throw std::invalid_argument("SparseToeplitz::Add expects SparseToeplitz");

    // Build union of diagonals (both arrays assumed sorted ascending)
    std::vector<int> diagUnion;
    int ia = 0, ib = 0;
    while (ia < Num_Diags || ib < o->Num_Diags) {
        if (ia < Num_Diags && (ib == o->Num_Diags || Diags[ia] < o->Diags[ib])) {
            diagUnion.push_back(Diags[ia++]);
        } else if (ib < o->Num_Diags && (ia == Num_Diags || o->Diags[ib] < Diags[ia])) {
            diagUnion.push_back(o->Diags[ib++]);
        } else { // equal
            diagUnion.push_back(Diags[ia]); ++ia; ++ib;
        }
    }

    int n = (int)diagUnion.size();
    SparseToeplitz* R = new SparseToeplitz(Num_Rows, Num_Cols, n);
    for (int i = 0; i < n; ++i) {
        R->Diags[i] = diagUnion[i];
        R->Vals[i] = 0.0;
    }

    // add values from this
    for (int i = 0; i < Num_Diags; ++i) {
        int d = Diags[i];
        auto it = std::lower_bound(diagUnion.begin(), diagUnion.end(), d);
        int idx = (int)std::distance(diagUnion.begin(), it);
        R->Vals[idx] += Vals[i];
    }
    // add values from other
    for (int i = 0; i < o->Num_Diags; ++i) {
        int d = o->Diags[i];
        auto it = std::lower_bound(diagUnion.begin(), diagUnion.end(), d);
        int idx = (int)std::distance(diagUnion.begin(), it);
        R->Vals[idx] += o->Vals[i];
    }

    return R;
}