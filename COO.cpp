#include "COO.h"

COO::COO(int length, int num_rows, int num_cols)
{
    Length = length;
    Array = new tuple[Length];
    Num_Rows = num_rows;
    Num_Cols = num_cols;
}

COO::COO(const std::initializer_list<tuple>& list, int num_rows, int num_cols)
: COO((int)list.size(), num_rows, num_cols)
{
    std::uninitialized_copy(list.begin(), list.end(), Array);

}

virtual Vectord COO::operator*(const Vectord vect)
{
    Vectord result(Num_Rows);
    double value;
    int row_ind;
    int col_ind;
    for(int i=0;i<Length;i++)
    {
        value = std::get<0>(Array[i]);
        row_ind = std::get<1>(Array[i]);
        col_ind = std::get<2>(Array[i]);
        result.Vec[row_ind] += value*vect.Vec[col_ind];
    }
    return result;
}

~COO()
{ 
    delete[] Array;
}

void COO::print()
{
    std::cout << '[';
    for(int i=0;i<Length;i++)
    {
        std::cout << '(' << std::get<0>(Array[i]) << ", "
                            << std::get<1>(Array[i]) << ", "
                            << std::get<2>(Array[i]) << ')';
    }
    std::cout << ']' << std::endl;
}

COO::COO(SparseToeplitz ST)
{
    Length = 0;
    Num_Rows = ST.Num_Rows;
    Num_Cols = ST.Num_Cols;

    //Determine number of non-zero values
    for(int q = 0; q<ST.Length;q++)
    {
        int f = ST.Diags[q];
        if(f==0)
        {
            Length += std::min(Num_Cols,Num_Rows);
        }
        else if(f>0)
        {
            Length += Num_Cols-f;
        }
        else
        {
            Length += Num_Rows+f;
        }
    }

    Array = new tuple[Length];
    int c = 0;

    //Fill Array with correct tuples
    for(int i=0;i<Num_Rows;i++)
    {
        for(int j = 0;j<ST.Length;j++)
        {
            if(Num_Cols > ST.Diags[j]+i && ST.Diags[j]+i >= 0)
            {
                Array[c] = std::make_tuple(ST.Vals[j],i,ST.Diags[j]+i);
                c++;
            }
        }
    }
}