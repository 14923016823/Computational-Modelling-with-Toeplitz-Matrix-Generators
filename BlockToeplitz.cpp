#include <iostream>
#include <stdexcept>//allows throwing exceptions
#include <cmath>
#include <vector> //for std::vec
#include <cstdio>//alows the use of printf

class Vectord
{
public:
    std::vector<double> Vec;

    Vectord(int length = 0) : Vec(length, 0.0) {}

    int Length() const { return static_cast<int>(Vec.size()); }

    void PrintVector() const
    {
        for (int i = 0; i < Length(); ++i)
        {
            printf("vec[%d]=%f\n", i, Vec[i]);
        }
        printf("done\n");
    }

    void Sum(const Vectord& VecIn)
    {
        if (VecIn.Length() != Length())
        {
            throw std::invalid_argument("You can't sum vectors with different sizes");
        }
        for (int i = 0; i < Length(); ++i)
        {
            Vec[i] += VecIn.Vec[i];
        }
    }
};


class RecursiveToeplitz
{
public:
    int Num_Rows = 0;
    int Num_Cols = 0;
    bool IsRecursive = true;

    int Num_Diags = 0;
    int Num_Rows_SubMatrixes = 0;
    int Num_Cols_SubMatrixes = 0;

    std::vector<int> Diags;                    //
    std::vector<RecursiveToeplitz> Vals;       // sub-Toeplitz matrices (only if recursive)
    std::vector<double> ValsD;                 // diagonal values (only if leaf)

    // default constructor (needed for vector::resize)
    RecursiveToeplitz() : IsRecursive(false), Num_Rows(0), Num_Cols(0),
                          Num_Diags(0), Num_Rows_SubMatrixes(0), Num_Cols_SubMatrixes(0) {}

    // recursive constructor
    RecursiveToeplitz(int nrows, int ncols, int ndiags, int nrsub, int ncsub)
    {
        IsRecursive = true;
        Num_Rows_SubMatrixes = nrsub;
        Num_Cols_SubMatrixes = ncsub;
        Num_Rows = nrows;
        Num_Cols = ncols;
        Num_Diags = ndiags;
        Diags.resize(ndiags);
        Vals.resize(ndiags); // requires configuration later on
    }

    // leaf constructor (diagonals are scalar values)
    RecursiveToeplitz(int nrows, int ncols, int ndiags)
    {
        IsRecursive = false;
        Num_Rows = nrows;
        Num_Rows_SubMatrixes = 1;
        Num_Cols = ncols;
        Num_Cols_SubMatrixes = 1;
        Num_Diags = ndiags;
        Diags.resize(ndiags);
        ValsD.resize(ndiags, 0.0);
    }

    // Corrected copy-and-scale constructor
RecursiveToeplitz(const RecursiveToeplitz& other, double c)
{
    Num_Rows = other.Num_Rows;
    Num_Cols = other.Num_Cols;
    IsRecursive = other.IsRecursive;
    Num_Diags = other.Num_Diags;
    Num_Rows_SubMatrixes = other.Num_Rows_SubMatrixes;
    Num_Cols_SubMatrixes = other.Num_Cols_SubMatrixes;

    Diags = other.Diags;  // always the same

    if (!IsRecursive)
    {
        // leaf: scale every diagonal value
        ValsD.resize(Num_Diags);
        for (int i = 0; i < Num_Diags; i++)
            ValsD[i] = other.ValsD[i] * c;
    }
    else
    {
        // recursive: COPY the submatrices, then scale them
        Vals = other.Vals;              // <--- FIX 1 (deep copy via vector copy)
        for (int i = 0; i < Num_Diags; i++)
            Vals[i].ScalarMultiplication(c);   // <--- FIX 2 (apply scaling here)
    }
}


    //This constructor can copy another Matrix
    /*
    RecursiveToeplitz(const RecursiveToeplitz& other,double c)
    {
        Num_Rows=other.Num_Rows;
        Num_Cols=other.Num_Cols;
        IsRecursive=other.IsRecursive;
        Num_Diags=other.Num_Diags;
        Num_Rows_SubMatrixes=other.Num_Rows_SubMatrixes;
        Num_Cols_SubMatrixes=other.Num_Cols_SubMatrixes;
        if(IsRecursive==false)
        {
            Diags.resize(Num_Diags);
            ValsD.resize(Num_Diags);
            for(int i=0;i<Num_Diags;i++)
            {       
                Diags[i]=other.Diags[i];
                ValsD[i]=c*other.ValsD[i];
            }
            //Vals.clear();
        }
        else
        {
            Diags.resize(Num_Diags);
            Vals.resize(Num_Diags);
            for(int i=0;i<Num_Diags;i++)
            {
                Diags[i]=other.Diags[i];
                Vals[i]=RecursiveToeplitz(other.Vals[i],c);
            }
            //ValsD.clear();
        }
    }

    */
    Vectord VectorMultiplication(const Vectord vec) const
    {
        if (!IsRecursive)
        {
            return VectorMultiplicationD(vec);
        }

        if (vec.Length() != Num_Cols)
        {
            throw std::invalid_argument("Vector length and matrix columns don't match (recursive).");
        }

        Vectord result(Num_Rows);
        int Num_blockrows = Num_Rows / Num_Rows_SubMatrixes;
        int Num_blockcols = Num_Cols / Num_Cols_SubMatrixes;

        for (int blockrow = 0; blockrow < Num_blockrows; blockrow++)
        {
            //int row = blockrow *;
            Vectord subres(Num_Rows_SubMatrixes); // initialized to zeros
            // accumulate contributions from each diagonal/block
            for (int j = 0; j < Num_Diags; j++)
            {
                int blockcol = Diags[j] + blockrow;
                // ensure the whole sub-block fits in input vector
                if (blockcol >= 0 && blockcol<Num_blockcols)
                {
                    Vectord subinput(Num_Cols_SubMatrixes);
                    int col_start=blockcol*Num_Cols_SubMatrixes;
                    for (int k = 0; k < Num_Cols_SubMatrixes; ++k)
                    {
                        subinput.Vec[k] = vec.Vec[col_start + k]; 
                    }
                    Vectord subsubres = Vals[j].VectorMultiplication(subinput);
                    subres.Sum(subsubres);
                }
            }
            // copy accumulated block into result
            int row_start=blockrow*Num_Rows_SubMatrixes;
            for (int k = 0; k < Num_Rows_SubMatrixes; ++k)
            {
                result.Vec[row_start + k] = subres.Vec[k];
            }
        }
        return result;
    }

    Vectord VectorMultiplicationD(const Vectord& vec) const
    {
        if (vec.Length() != Num_Cols)
        {
            throw std::invalid_argument("Vector and Matrix size dont match (leaf).");
        }
        Vectord result(Num_Rows);
        for (int i = 0; i < Num_Rows; i++)
        {
            double acc = 0.0;
            for (int j = 0; j < Num_Diags; j++)
            {
                int idx = Diags[j] + i;
                if (idx >= 0 && idx < Num_Cols)
                {
                    acc += vec.Vec[idx] * ValsD[j];
                }
            }
            result.Vec[i] = acc;
        }
        return result;
    }

    void ScalarMultiplication(double c)
    {
        if(IsRecursive==false)
        {
            for(int i=0;i<Num_Diags;i++)
            {
                ValsD[i]*=c;
            }
        }
        else
        {
            for(int i=0;i<Num_Diags;i++)
            {
                Vals[i].ScalarMultiplication(c);
            }
        }
    }

    void Kronecker(RecursiveToeplitz B)//if you add a new matrix at the bottom of the chain every Num_Rows needs to be changed
    {
        Num_Rows*=B.Num_Rows;
        Num_Cols*=B.Num_Cols;
        Num_Cols_SubMatrixes*=B.Num_Cols;
        Num_Rows_SubMatrixes*=B.Num_Rows;
        if(IsRecursive==true)
        {
            for(int i=0;i<Num_Diags;i++)
            {
                Vals[i].Kronecker(B);
            }
        }
        else
        {   
            IsRecursive=true; //after kronecker this will become a Blockmatrix
            Vals.resize(Num_Diags); //we will now use the Valsmatrix instead of the ValsD
            for(int i=0;i<Num_Diags;i++)
            {
                Vals[i]=RecursiveToeplitz(B,ValsD[i]);
            }
            //ValsD.clear();
        }


    }

};

void TestKronecker() //This function tests the kornicker function by generating 2d laplacian with dirichlet boundary conditions
{
    int n=3;
    RecursiveToeplitz L1d=RecursiveToeplitz(n,n,3); //Laplacian1d
    L1d.Diags[0]=-1;
    L1d.Diags[1]=0;
    L1d.Diags[2]=1;
    L1d.ValsD[0]=-1;
    L1d.ValsD[1]=2;
    L1d.ValsD[2]=-1;
    RecursiveToeplitz I=RecursiveToeplitz(n,n,1); //identity matrix
    I.Diags[0]=0;
    I.ValsD[0]=1;
    I.Kronecker(L1d);
    printf("%d,%d",I.Num_Rows,I.Num_Cols);
   

    Vectord input(n*n);
    input.Vec[0] = 1.0;
    //input.Vec[1] = 1.0;
    input.PrintVector();
    Vectord result = I.VectorMultiplication(input);
    result.PrintVector();
    fflush(stdout);
    
    for(int i=1;i<(n*n);i++)
    {
        input.Vec[i-1]=0;
        input.Vec[i]=1.0;
        Vectord result = I.VectorMultiplication(input);
        result.PrintVector();
    }
    
    

};

void TestMatmul()
{
    Vectord input(12);
    // default all zeros;

    
    //input.PrintVector();

    RecursiveToeplitz inc_1 = RecursiveToeplitz(4, 3, 2);
    RecursiveToeplitz inc0   = RecursiveToeplitz(4, 3, 2);
    RecursiveToeplitz inc1   = RecursiveToeplitz(4, 3, 2);

    inc_1.ValsD = {5.0, 5.0};
    inc_1.Diags  = {-1, 0};
    inc0.ValsD   = {-1.0, 1.0};
    inc0.Diags   = {-1, 0};
    inc1.ValsD   = {-2.0, 4.0};
    inc1.Diags   = {-1, 0};

    RecursiveToeplitz inc = RecursiveToeplitz(16, 12, 3, 4, 3);
    inc.Vals[0] = inc_1;
    inc.Vals[1] = inc0;
    inc.Vals[2] = inc1;
    inc.Diags[0] = -1;
    inc.Diags[1] = 0;
    inc.Diags[2] = 1;

    input.Vec[0] = 1.0;
    //input.Vec[1] = 1.0;
    input.PrintVector();
    Vectord result = inc.VectorMultiplication(input);
    result.PrintVector();
    
    for(int i=1;i<12;i++)
    {
        input.Vec[i-1]=0;
        input.Vec[i]=1.0;
        Vectord result = inc.VectorMultiplication(input);
        result.PrintVector();
    }

}



int main()//For testing matrix multiplication
{
    //TestMatmul();
    TestKronecker();


    

    return 0;
}
