#include "COO.h"

COO::COO(int num_rows, int num_cols, int num_vals)
{
    Num_Vals = num_vals;
    Array = new tuple[Num_Vals];
    Num_Rows = num_rows;
    Num_Cols = num_cols;
}

COO::COO(int num_rows, int num_cols, const std::initializer_list<tuple>& list)
: COO(num_rows, num_cols, (int)list.size())
{
    std::uninitialized_copy(list.begin(), list.end(), Array);

}

Vectord COO::operator*(Vectord& vect)
{
    Vectord result(Num_Rows);
    double value;
    int row_ind;
    int col_ind;
    int i;
    #pragma omp parallel for private(i, row_ind, col_ind, value)
    for(i=0;i<Num_Vals;i++)
    {
        row_ind = std::get<0>(Array[i]);
        col_ind = std::get<1>(Array[i]);
        value = std::get<2>(Array[i]);
        result.Vec[row_ind] += value*vect.Vec[col_ind];
    }
    return result;
}

COO::~COO()
{ 
    delete[] Array;
    Array = nullptr;
}

void COO::print()
{
    std::cout << '[';
    for(int i=0;i<Num_Vals;i++)
    {
        std::cout << '(' << std::get<0>(Array[i]) << ", "
                            << std::get<1>(Array[i]) << ", "
                            << std::get<2>(Array[i]) << ')';
    }
    std::cout << ']' << std::endl;
}

COO::COO(SparseToeplitz& ST)
{
    Num_Vals = 0;
    Num_Rows = ST.rows();
    Num_Cols = ST.cols();

    //Determine number of non-zero values
    for(int q = 0; q<ST.Num_Diags;q++)
    {
        int f = ST.Diags[q];
        if(f==0)
        {
            Num_Vals += std::min(Num_Cols,Num_Rows);
        }
        else if(f>0)
        {
            Num_Vals += Num_Cols-f;
        }
        else
        {
            Num_Vals += Num_Rows+f;
        }
    }

    Array = new tuple[Num_Vals];
    int c = 0;

    //Fill Array with correct tuples
    for(int i=0;i<Num_Rows;i++)
    {
        for(int j = 0;j<ST.Num_Diags;j++)
        {
            if(Num_Cols > ST.Diags[j]+i && ST.Diags[j]+i >= 0)
            {
                Array[c] = std::make_tuple(i,ST.Diags[j]+i,ST.Vals[j]);
                c++;
            }
        }
    }
}

void COO::operator*=(double c) 
{
    int i;
    #pragma omp parallel for private(i)
    for (i=0;i<Num_Vals;i++)
    {
        std::get<2>(Array[i])*=c;
    }
}

COO::COO(COO& other,double c)
{

    Num_Rows=other.Num_Rows;
    Num_Cols=other.Num_Cols;
    Num_Vals=other.Num_Vals;
    Array = new tuple[Num_Vals];
    int i;
    #pragma omp parallel for private(i)
    for(i=0;i<Num_Vals;i++)
    {
        std::get<0>(Array[i])=std::get<0>(other.Array[i]);
        std::get<1>(Array[i])=std::get<1>(other.Array[i]);
        std::get<2>(Array[i])=std::get<2>(other.Array[i])*c;
    }

}

double COO::operator()(int i, int j) const
{
    if (i < 0 || i >= Num_Rows || j < 0 || j >= Num_Cols)
    {
        throw std::invalid_argument("Index outside of matrix boundaries");
    }
    for (int d = 0; d < Num_Vals; ++d) 
    {
        if (std::get<0>(Array[d]) == i && std::get<1>(Array[d]) == j)
            return std::get<2>(Array[d]);
    }
    return 0.0;
}

Matrix* COO::Clone(double c)
{
   return new COO(*this, c);
}

Matrix* COO::Kronecker(Matrix& B)
{
   
    BlockCOO* result = new BlockCOO(B.rows()*Num_Rows,B.cols()*Num_Cols,Num_Vals);
    
    int i;
    #pragma omp parallel for private(i)
    for(i=0;i<Num_Vals;i++)
    {
        std::get<0>(result->Array[i])=std::get<0>(Array[i]);
        std::get<1>(result->Array[i])=std::get<1>(Array[i]);
        std::get<2>(result->Array[i])=B.Clone(std::get<2>(Array[i]));
        std::cout << std::get<2>(Array[i]) << "\n";
    }
    
    //printf("vals[i],%f\n",std::get<0>(result->Array)); //This print statement does not work
    return result;
}

Matrix* COO::negativeTranspose()
{ //This seems very inefficient, but I do not know how to get it working otherwise
    COO* negTrans = new COO(Num_Cols, Num_Rows, Num_Vals);
    int Rows[Num_Cols+1]; //Array to count entries per row in transposed matrix
    
    Rows[0] = 0;
#pragma omp parallel for
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
                    std::get<2>(negTrans->Array[Rows[c]+n]) = -std::get<2>(Array[i]);
                    n++;
                }
            }
        }
    } 
    return negTrans;
}

Matrix* COO::printFullMatrix()
{
    std::cout << "\nFull Dense Expansion (" << Num_Rows << "x" << Num_Cols << ")\n";
    for (int i = 0; i < Num_Rows; i++) {
        for (int j = 0; j < Num_Cols; j++)
        {
            std::cout << std::setw(4) << this->operator()(i,j);
        }
        std::cout << "\n";
    }
    return nullptr;
}

Matrix* COO::Add(Matrix& other)
{
    return nullptr;
}