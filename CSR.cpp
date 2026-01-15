#include "CSR.h"

CSR::CSR(double* vals, int* cols, int* rows, int num_vals, int num_rows, int num_cols)
{
    Vals = vals;
    Cols = cols;
    Rows = rows;
    Num_Vals = num_vals;
    Num_Cols = num_cols;
    Num_Rows = num_rows;
}

CSR::CSR(int num_rows, int num_cols, int num_vals)
{
    Vals = new double[num_vals];
    Cols = new int[num_vals];
    Rows = new int[num_rows+1];
    Num_Vals = num_vals;
    Num_Cols = num_cols;
    Num_Rows = num_rows;
    for(int q = 0;q<num_rows+1;q++)
        Rows[q] = 0;
    for(int q = 0;q<num_vals;q++)
        Vals[q] = 0;
    for(int q = 0;q<num_vals;q++)
        Cols[q] = 0;
}

Vectord CSR::operator*(Vectord& vect)
{   
    int len = vect.Length;
    if (len != Num_Cols) 
    {
        throw std::invalid_argument("Vector and Matrix size don't match");
    }
    Vectord result(Num_Rows);
    //int i=0;
    //int j=0;
    #pragma omp parallel for //num_threads(12)// private(i,j) 
    for(int i=0;i<Num_Rows;i++)
    {
        //printf(" i = %d\n",i);
        //printf("%d\n",Rows[i+1]>Rows[i]);
        int res_i = 0;
        #pragma omp parallel for reduction(+ : res_i)
        for(int j=Rows[i];j<Rows[i+1];j++)
        {
            //printf(" j = %d\n",j);
            //printf("Cols[%d]=%d\n",c,Cols[c]);
            //printf("vect[%d]=%f\n",Cols[c],vect.Vec[Cols[c]]);
            res_i = vect.Vec[Cols[j]]*Vals[j];
            //printf("r[i]_j = %f\n",vect.Vec[Cols[c]]*Vals[c]);
            //printf("r[i]_tot = %f\n",result.Vec[i]);
        }
        result.Vec[i] += res_i;
    }
    return result;
}

//Destructor
CSR::~CSR()
{ 
    delete[] Vals;
    delete[] Rows;
    delete[] Cols;
    Vals = nullptr;
    Rows = nullptr;
    Cols = nullptr;
}

void CSR::print()
{   
    std::cout << "Values: [";
    for(int i = 0;i<Num_Vals;i++)
    {
        std::cout << Vals[i];
        if(i<Num_Vals-1)
        std::cout << ", ";
    }
    std::cout << ']' << std::endl;

    std::cout << "Columns: [";
    for(int i = 0;i<Num_Vals;i++)
    {
        std::cout << Cols[i];
        if(i<Num_Vals-1)
        std::cout << ", ";
    }
    std::cout << ']' << std::endl;

    std::cout << "Rows: [";
    for(int i = 0;i<Num_Rows+1;i++)
    {
        std::cout << Rows[i];
        if(i<Num_Rows)
        std::cout << ", ";
    }
    std::cout << ']' << std::endl;
}

CSR::CSR(CSR& other, double c)
{
    Num_Cols = other.Num_Cols;
    Num_Rows=other.Num_Rows;
    Num_Vals=other.Num_Vals;
    Cols = new int[Num_Vals];
    Rows = new int[Num_Rows+1];
    Vals = new double[Num_Vals];
    for(int i=0;i<Num_Vals;i++)
    {
        Cols[i]=other.Cols[i];
        Vals[i]=c*other.Vals[i];
    }
    for(int i=0;i<Num_Rows+1;i++)
    {
        Rows[i]=other.Rows[i];
    }
}

//Convert Sparse Toeplitz matrix to CSR matrix
CSR::CSR(SparseToeplitz& ST)
{
    Num_Vals = 0;
    Num_Cols = ST.cols();
    Num_Rows = ST.rows();
    int f;
    //Get the length of the Vals and Cols arrays
    for(int q = 0; q<ST.Num_Diags;q++)
    {
        f = ST.Diags[q];
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
    Vals = new double[Num_Vals];
    Cols = new int[Num_Vals];
    Rows = new int[Num_Rows+1];
    int c = 0;
    Rows[0] = 0;
    //Fill Vals and Cols arrays with correct values
    for(int i=0;i<Num_Rows;i++)
    {
        for(int j = 0;j<ST.Num_Diags;j++)
        {
            if(Num_Cols > ST.Diags[j]+i && ST.Diags[j]+i >= 0)
            {
                Vals[c] = ST.Vals[j];
                Cols[c] = ST.Diags[j]+i;
                c++;
            }
        }
        Rows[i+1] = c;
    }
}

void CSR::operator*=(double c) 
{
    for (int i=0;i<Num_Vals;i++)
    {
        Vals[i]*=c;
    }
}

Matrix* CSR::Kronecker(Matrix& B)
{
   
    BlockCSR* result = new BlockCSR(B.rows()*Num_Rows,B.cols()*Num_Cols,Num_Vals);
    int i;
    result->Rows[0]=0;
    #pragma omp parallel for private(i)
    for(i=0;i<Num_Vals;i++)
    {
        result->Cols[i] = Cols[i];
        result->Rows[i+1] = Rows[i+1];
        result->Vals[i] = B.Clone(Vals[i]);
    }
    return result;
}

double CSR::operator()(int i, int j) const
{
    if (i < 0 || i >= Num_Rows || j < 0 || j >= Num_Cols)
        throw std::invalid_argument("Index outside of matrix boundaries");
    for(int k = Rows[i]; k < Rows[i+1]; k++) 
    {
        if (Cols[k] == j)
            return Vals[k];  
    }
    return 0.0;
}

Matrix* CSR::Clone(double c)
{
   return new CSR(*this, c);
}

Matrix* CSR::negativeTranspose()
{ 
    CSR* negTrans = new CSR(Num_Cols, Num_Rows, Num_Vals);
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
                    negTrans->Vals[negTrans->Rows[c]+m] = -Vals[i];
                    negTrans->Cols[negTrans->Rows[c]+m] = r;
                    m++;
                }
            }
        }
    }
    /*
    int m = 0;
    for(int c=0;c<Num_Cols;c++) //Loop over columns of original matrix
    {
        for(int r=0;r<Num_Rows;r++) //Loop over rows of original matrix
        {
            for(int i=Rows[r];i<Rows[r+1];i++) //Loop over values in current row
            {
                if(Cols[i]==c) //Placing this in the i-loop above and replacing m by n does not work
                {
                    negTrans->Vals[m] = -Vals[i];
                    negTrans->Cols[m] = r;
                    m++;
                }
            }
        }
    } //Sequential version*/
    return negTrans;
}

void CSR::printFullMatrix()
{
    std::cout << "\nFull Dense Expansion (" << Num_Rows << "x" << Num_Cols << ")\n";
    for (int i = 0; i < Num_Rows; ++i) {
        for (int j = 0; j < Num_Cols; ++j)
        {
            std::cout << std::setw(4) << this->operator()(i,j);
        }
        std::cout << "\n";
    }
}