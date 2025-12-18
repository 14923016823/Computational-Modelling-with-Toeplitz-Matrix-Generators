#include "COO.h"

COO::COO(int num_vals, int num_rows, int num_cols)
{
    Num_Vals = num_vals;
    Array = new tuple[Num_Vals];
    Num_Rows = num_rows;
    Num_Cols = num_cols;
}

COO::COO(const std::initializer_list<tuple>& list, int num_rows, int num_cols)
: COO((int)list.size(), num_rows, num_cols)
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
        value = std::get<0>(Array[i]);
        row_ind = std::get<1>(Array[i]);
        col_ind = std::get<2>(Array[i]);
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
                Array[c] = std::make_tuple(ST.Vals[j],i,ST.Diags[j]+i);
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
        std::get<0>(Array[i])*=c;
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
        std::get<0>(Array[i])*=c;
    }

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
        std::get<0>(result->Array[i])=B.Clone(std::get<0>(Array[i]));
        std::get<1>(result->Array[i])=std::get<1>(Array[i]);
        std::get<2>(result->Array[i])=std::get<2>(Array[i]);
    }
    
    //printf("vals[i],%f\n",std::get<0>(result->Array)); //This print statement does not work
    return result;
}