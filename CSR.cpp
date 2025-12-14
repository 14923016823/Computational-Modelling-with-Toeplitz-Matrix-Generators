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

Vectord CSR::operator*(Vectord& vect)
{   
    int len = vect.Length;
    if (len != Num_Cols) 
    {
        throw std::invalid_argument("Vector and Matrix size dont match");
    }
    Vectord result(Num_Rows);
    int c = 0;
    int d = 0;
    for(int i=0;i<Num_Rows;i++)
    {
        if(Rows[i+1]>Rows[i])
        {
            d = Rows[i+1]-Rows[i];
            for(int j=0;j<d;j++)
            {
                result.Vec[i] += vect.Vec[Cols[c]]*Vals[c];
                c++;
            }
        }
    }
    return result;
}

//Destructor
CSR::~CSR()
{ 
    delete[] Vals;
    delete[] Rows;
    delete[] Cols;
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
    //Get the length of the Vals and Cols arrays
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
    
    for(int i=0;i<Num_Vals;i++)
    {
        result->Cols[i]=Cols[i];
        result->Rows[i]=Rows[i];
        result->Vals[i] = B.Clone(Vals[i]);   
    }
    
    printf("vals[i],%f\n",result->Vals);
    return result;
}

Matrix* CSR::Clone(double c)
{
   return new CSR(*this, c);
}