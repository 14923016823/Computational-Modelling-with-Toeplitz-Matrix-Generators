#include "CSR.h"

class CSR: public Matrix
{
public:
    CSR(double* vals, int* cols, int* rows, int length, int num_rows, int num_cols)
    {
        Vals = vals;
        Cols = cols;
        Rows = rows;
        Length = length;
        Num_Cols = num_cols;
        Num_Rows = num_rows;
    }

    virtual Vectord operator*(const Vectord vect) const
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
    ~CSR()
    { 
        delete[] Vals;
        delete[] Rows;
        delete[] Cols;
    }

    double* Vals;
    int* Cols;
    int* Rows;
    int Length;

    void print() const
    {
        std::cout << "Values: [";
        for(int i = 0;i<Length;i++)
        {
            std::cout << Vals[i];
            if(i<Length-1)
            std::cout << ", ";
        }
        std::cout << ']' << std::endl;

        std::cout << "Columns: [";
        for(int i = 0;i<Length;i++)
        {
            std::cout << Cols[i];
            if(i<Length-1)
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

    //Convert Sparse Toeplitz matrix to CSR matrix
    CSR(SparseToeplitz ST)
    {
        Length = 0;
        Num_Cols = ST.Num_Cols;
        Num_Rows = ST.Num_Rows;
        //Get the length of the Vals and Cols arrays
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
        Vals = new double[Length];
        Cols = new int[Length];
        Rows = new int[Num_Rows+1];
        int c = 0;
        Rows[0] = 0;
        //Fill Vals and Cols arrays with correct values
        for(int i=0;i<Num_Rows;i++)
        {
            for(int j = 0;j<ST.Length;j++)
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
};