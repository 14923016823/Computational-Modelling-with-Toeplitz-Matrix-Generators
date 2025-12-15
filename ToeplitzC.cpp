#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>




double sourcefunc(double x, double y,double z) 
{
    double distance = x * x + y * y + z * z;
    return exp(-4*distance);
}

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

    void PrintVector()
    {
        
        for(int i=0;i<Length;i++)
        {
            printf("vec[%d]=%f\n",i,Vec[i]);

        }
        printf("done\n");
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

class SparseToeplitz 
{
public:
    int Height;// the height of the matrix
    int Width; //the width of the matrix
    int Length;// the length of the Diags and Vals vectors
    int* Diags;//an int array with all nonzero diagonals (ordered smallest to largest)
    double* Vals;//a double array with the values at those diagonals (ordered the same way as the int array)
    SparseToeplitz(int height, int width, int length)
    {
        if (length > height + width - 1)
        {
            throw std::invalid_argument("There are too many diagonals for Toeplitz matrix of this size");
        }
        Height = height;
        Width = width;
        Length = length;
        Diags = new int[Length];
        Vals = new double[Length];
    }
    SparseToeplitz(int height, int width, int length,int* diags,double* vals)
    {
        if (length > height + width - 1)
        {
            throw std::invalid_argument("There are too many diagonals for Toeplitz matrix of this size");
        }
        Height = height;
        Width = width;
        Length = length;
        Diags = diags;
        Vals = vals;
    }

    Vectord VectorMultiplication(Vectord vec) 
    {
        int len = vec.Length;
        if (len != Width) 
        {
            throw std::invalid_argument("Vector and Matrix size dont match");
        }
        Vectord result = Vectord(len);
        for (int i = 0;i < Height;i++) 
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


};

class RecursiveToeplitz
{
    int Num_Rows;
    int Num_Cols;
    bool IsRecursive=true;

    //int* Diags=nullptr; //The placement of the main diagonal of the ToeplitzMatrixes (ordered from largest to smallest)
    //RecursiveToeplitz* Vals=nullptr;//we assume all submatrixes have the same shape
    int Num_Diags; //The number of nonzero Toeplitzdiagonals

    int Num_Rows_SubMatrixes;
    int Num_Cols_SubMatrixes;
    
    std::vector<int> Diags;                    // diagonal offsets
    std::vector<RecursiveToeplitz> Vals;       // sub-Toeplitz matrices
    std::vector<double>ValsD;

    RecursiveToeplitz(int nrows,int ncols, int ndiags,int nrsub, int ncsub)
    {
        IsRecursive=true;
        Num_Rows_SubMatrixes=nrsub;
        Num_Cols_SubMatrixes=ncsub;
        Num_Rows=nrows;
        Num_Cols=ncols;
        Num_Diags=ndiags;
        Diags.resize(ndiags);
        Vals.resize(ndiags); 
    }

    RecursiveToeplitz(int nrows,int ncols, int ndiags)
    {
        IsRecursive=false;
        Num_Rows=nrows;
        Num_Rows_SubMatrixes=Num_Rows;
        Num_Cols=ncols;
        Num_Cols_SubMatrixes=Num_Cols;
        Diags.resize(ndiags);
        ValsD.resize(ndiags);
    }
    

    Vectord VectorMultiplication(Vectord vec)
    {
        if(IsRecursive=false)
        {
            return VectorMultiplicationD(vec);
        }
        Vectord result=Vectord(Num_Rows);
        for(int i=0;i<(Num_Rows/Num_Rows_SubMatrixes);i++)
        {
            int row=i*Num_Rows_SubMatrixes;

            Vectord subres=Vectord(Num_Rows_SubMatrixes);
            Vectord subinput=Vectord(Num_Cols_SubMatrixes);
            for(int j=0;j<Num_Diags;j++)
            {
                int col=Diags[j]+i*Num_Cols_SubMatrixes;
                if(col>=0 && col<Num_Cols)
                {
                    for(int k=0;k<Num_Cols_SubMatrixes;k++)
                    {
                        subinput.Vec[k]=vec.Vec[col+k]; //Copying a relevant part of the vector to subinput
                    }
                    Vectord subsubres=Vals[j].VectorMultiplication(subinput); //The result of small part of th vector with a single matrix
                    subres.Sum(subsubres);//adding these results together
                }
                for(int k=0;k<Num_Rows_SubMatrixes;k++)
                {
                    result.Vec[row+k]=subres.Vec[k]; //Copying a relevant part of the output vector to 
                }
            }
        }
        return result;
    }

    Vectord VectorMultiplicationD(Vectord vec) 
    {
        if (vec.Length != Num_Cols) 
        {
            throw std::invalid_argument("Vector and Matrix size dont match");
        }
        Vectord result = Vectord(Num_Rows);
        for (int i = 0;i < Num_Rows;i++) 
        {
            for (int j = 0;j < Num_Diags;j++) 
            {
                if (Diags[j] + i >= 0 && Diags[j] + i < Num_Cols)
                {
                    result.Vec[i] += vec.Vec[Diags[j] + i] * ValsD[j];
                }
            }
        }
        return result;
    }
};


int main()
{
    int n=10;
    Vectord StartVec = Vectord(n);
    for(int i=0;i<n;i++)
    {
        StartVec.Vec[i]=i;
    }
    int diags[2]={-(n-1),1};
    double vals[2]={1,1};
    SparseToeplitz LeftShift=SparseToeplitz(10,10,2,diags,vals);
    Vectord leftshifted=LeftShift.VectorMultiplication(StartVec);

    StartVec.PrintVector();
    leftshifted.PrintVector();
    return 0;
}