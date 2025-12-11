#include <iostream>
#include <exception>
#include <cmath>
#include <initializer_list>
#include <tuple> //khfjyfgukyfjg

typedef std::tuple<double,int,int> tuple;

class Vectord
{
public:
    int Length;
    double* Vec;
    Vectord(int length)
    {
        Length = length;
        Vec = new double[Length];
    }

    Vectord(const std::initializer_list<double>& list)
    : Vectord((int)list.size())
    {
        std::uninitialized_copy(list.begin(), list.end(), Vec);
    }

    void print()
    {
        std::cout << "[";
        for(int i=0;i<Length;i++)
        {
            std::cout << Vec[i];
            if(i<Length-1)
            std::cout << ", ";
        }
        std::cout << ']' << std::endl;
    }

    void Sum(Vectord VecIn)
    {
        if(VecIn.Length=!Length)\
        {
            throw std::invalid_argument("You can't sum vectors with different sizes");
        }
        for(int i=0;i<Length;i++)
        {
            Vec[i]+=VecIn.Vec[i];
        }   
    }
};


class Matrix
{
public:
    int Num_Rows;
    int Num_Cols;
    Matrix()
    { }
    
    virtual ~Matrix(){}

    virtual Vectord operator*(const Vectord vec) const = 0;

    virtual void print() const = 0;
};


class SparseToeplitz: public Matrix
{
public:
    int Length;// the length of the Diags and Vals vectors
    int* Diags;//an int array with all nonzero diagonals (ordered smallest to largest)
    double* Vals;//a double array with the values at those diagonals (ordered the same way as the int array)
    SparseToeplitz(int num_rows, int num_cols, int length)
    {
        if (length > num_rows + num_cols - 1)
        {
            throw std::invalid_argument("There are too many diagonals for Toeplitz matrix of this size");
        }
        Num_Rows = num_rows;
        Num_Cols = num_cols;
        Length = length;
        Diags = new int[Length];
        Vals = new double[Length];
    }
    SparseToeplitz(int num_rows, int num_cols, int length,int* diags,double* vals)
    {
        if (length > num_rows + num_cols - 1)
        {
            throw std::invalid_argument("There are too many diagonals for Toeplitz matrix of this size");
        }
        Num_Rows = num_rows;
        Num_Cols = num_cols;
        Length = length;
        Diags = diags;
        Vals = vals;
    }

    virtual Vectord operator*(const Vectord vec) const
    {
        int len = vec.Length;
        if (len != Num_Cols) 
        {
            throw std::invalid_argument("Vector and Matrix size dont match");
        }
        Vectord result = Vectord(len);
        for (int i = 0;i < Num_Rows;i++) 
        {
            for (int j = 0;j < Length;j++) 
            {
                if (Diags[j] + i >= 0 && Diags[j] + i < len)
                {
                    result.Vec[i] += vec.Vec[Diags[j] + i] * Vals[j];
                }
            }
        }
        return result;
    }

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

        std::cout << "Diagonals: [";
        for(int i = 0;i<Length;i++)
        {
            std::cout << Diags[i];
            if(i<Length-1)
            std::cout << ", ";
        }
        std::cout << ']' << std::endl;
    }
};

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
    { }

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

class COO: public Matrix
{
public:
    tuple* Array;
    int Length;

    COO(int length, int num_rows, int num_cols)
    {
        Length = length;
        Array = new tuple[Length];
        Num_Rows = num_rows;
        Num_Cols = num_cols;
    }

    COO(const std::initializer_list<tuple>& list, int num_rows, int num_cols)
    : COO((int)list.size(), num_rows, num_cols)
    {
        std::uninitialized_copy(list.begin(), list.end(), Array);

    }

    virtual Vectord operator*(const Vectord vect) const
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

    void print() const
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

    COO(SparseToeplitz ST)
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
};

int main()
{
    const int length = 3;
    double vals[length] = {1,4,6};
    int num_cols = 4;
    const int num_rows = 5;
    int cols[length] = {0,3,2};
    int rows[num_rows+1] = {0,2,2,2,2,3};

    tuple t1 = std::make_tuple(1,0,0);
    tuple t2 = std::make_tuple(4,0,3);
    tuple t3 = std::make_tuple(6,4,2);
    COO B({t1,t2,t3},num_rows,num_cols);
    B.print();

    CSR A(vals, cols, rows, length, num_rows, num_cols);
    Vectord x(num_cols);
    
    for(int i=0;i<num_cols;i++)
    {
        x.Vec[i]=i;
    }
    x.print();
        
    A.print();

    Vectord y = A*x;
    y.print();

    Vectord z(A*x);
    z.print();

    const int STlength = 3;
    int STdiags[STlength] = {-2,-1,3};
    double STvals[STlength] = {5,1,4};
    int STwidth = 4;
    int STheight = 3;
    SparseToeplitz C(STheight, STwidth, STlength, STdiags, STvals);
    C.print();

    CSR C_CSR(C);
    C_CSR.print();

    COO C_COO(C);
    C_COO.print();

    return 0;
}