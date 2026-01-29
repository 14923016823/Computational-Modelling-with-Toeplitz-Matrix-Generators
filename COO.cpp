#include "COO.h"

COO::COO(int num_rows, int num_cols, int num_vals)
{
    Num_Vals = num_vals;
    Array = new tuple[Num_Vals];
    Num_Rows = num_rows;
    Num_Cols = num_cols;
    for(int q = 0;q<num_vals;q++)
        Array[q] = std::make_tuple(0,0,0);
}

COO::COO(int num_rows, int num_cols, const std::initializer_list<tuple>& list)
: COO(num_rows, num_cols, (int)list.size())
{
    std::uninitialized_copy(list.begin(), list.end(), Array);
}

Vectord COO::operator*(Vectord& vect)
{
    int len = vect.len();
    if (len != Num_Cols) 
    {
        throw std::invalid_argument("Vector and Matrix size don't match");
    }
    Vectord result(Num_Rows);
    //double value;
    //int row_ind;
    //int col_ind;
    //int i;
    //#pragma omp parallel for //private(i, row_ind, col_ind, value)
    for(int i=0;i<Num_Vals;i++)
    {
        double res_i = std::get<2>(Array[i])*vect.Vec[std::get<1>(Array[i])];
        //#pragma omp atomic update
        result.Vec[std::get<0>(Array[i])] += res_i;
        //printf("i = %d\n",i);
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
    //#pragma omp parallel for private(i)
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
    //#pragma omp parallel for private(i)
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
    //#pragma omp parallel for private(i)
    for(i=0;i<Num_Vals;i++)
    {
        std::get<0>(result->Array[i])=std::get<0>(Array[i]);
        std::get<1>(result->Array[i])=std::get<1>(Array[i]);
        std::get<2>(result->Array[i])=B.Clone(std::get<2>(Array[i]));
        //std::cout << std::get<2>(Array[i]) << "\n";
    }
    
    //printf("vals[i],%f\n",std::get<0>(result->Array)); //This print statement does not work
    return result;
}

Matrix* COO::negativeTranspose()
{ //This seems very inefficient, but I do not know how to get it working otherwise
    COO* negTrans = new COO(Num_Cols, Num_Rows, Num_Vals);
    int Rows[Num_Cols+1]; //Array to count entries per row in transposed matrix
    for(int j=0;j<Num_Cols+1;j++)
        Rows[j] = 0;
//#pragma omp parallel for
    for(int c=0;c<Num_Cols;c++) //Go through columns of original matrix
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
    for(int c = 0;c<Num_Cols;c++)
    {
        Rows[c+1] += Rows[c];
        //std::cout << Rows[c+1] << '\n';
    }

    int n;
    //#pragma omp parallel for private(n)
    for(int c=0;c<Num_Cols;c++) //Cols of original matrix, so rows of transposed matrix
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
                    std::get<2>(negTrans->Array[n]) = -std::get<2>(Array[i]);
                    //std::get<2>(Array[i])->printFullMatrix();
                    //printf("finished c=%d, i=%d, n=%d\n",c, i,n);
                    n++;
                }
            }
        //}
    }
    return negTrans;
}

void COO::printFullMatrix()
{
    std::cout << "\nFull Dense Expansion (" << Num_Rows << "x" << Num_Cols << ")\n";
    for (int i = 0; i < Num_Rows; i++) {
        for (int j = 0; j < Num_Cols; j++)
        {
            std::cout << std::setw(6) << this->operator()(i,j);
        }
        std::cout << "\n";
    }
}

COO::COO(BlockToeplitz& BT)
{ //This conversion algorithm is very inefficient, but works fine for now
    Num_Vals = 0;
    Num_Rows = BT.rows();
    Num_Cols = BT.cols();

    for (int i = 0; i < Num_Rows; i++) {
        for (int j = 0; j < Num_Cols; j++)
        {
            if(BT.operator()(i,j) != 0.0)
            {
                Num_Vals++;
            }
        }
    }

    Array = new tuple[Num_Vals]; 
    int n = 0;
    for (int i = 0; i < Num_Rows; i++) {
        for (int j = 0; j < Num_Cols; j++)
        {
            if(BT.operator()(i,j) != 0.0)
            {
                Array[n] = std::tuple(i,j,BT.operator()(i,j));
                n++;
            }
        }
    }
/*
    //Determine number of non-zero values
    for(int q = 0; q<BT.Num_Diags;q++)
    {
        int f = BT.Diags[q];
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
        for(int j = 0;j<BT.Num_Diags;j++)
        {
            if(Num_Cols > BT.Diags[j]+i && BT.Diags[j]+i >= 0)
            {
                Array[c] = std::make_tuple(i,BT.Diags[j]+i,BT.Vals[j]);
                c++;
            }
        }
    }*/
}