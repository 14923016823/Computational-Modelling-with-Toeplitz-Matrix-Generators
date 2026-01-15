#include "Random.h"

int Normal_int(double mean, double var)
{
    std::random_device rd{};
    std::mt19937 gen{rd()};
 
    // Values near the mean are the most likely. Standard deviation
    // affects the dispersion of generated values from the mean.
    std::normal_distribution d{mean, var};
 
    // Draw a sample from the normal distribution and round it to an integer.
    return (int)std::lround(d(gen));
}

SparseToeplitz r_ST(int n_r, int n_c, int prob)
{
    
    int nvals = (n_r+n_c-1)/prob;
    if(nvals == 0)
        nvals = 1;
    int n=0;

    std::cout << nvals << "\n";
    SparseToeplitz D(n_r, n_c, nvals);
    int i=-n_r;
    bool twice = false;
    while(n<nvals)
    {
        if(std::experimental::randint(1,prob)==1)
        {
            D.Diags[n] = i;
            D.Vals[n] = std::experimental::randint(1,9);
            n++;
            if(D.Diags[n-1]>n_c-1 && n>1)
            {
                if(twice)
                {
                    n=0;
                    i=-n_r;
                }
                else
                {
                //std::cout << D.Diags[n-1] << " " << D.Diags[n-2] << '\n';
                n--;
                i = D.Diags[n-1];
                twice = true;
                }
            }
        }
        //if(i>n_c)
            //break;
        i++;
    }
    return D;
}

DiagonalMatrix r_Diag(int n_c, int n_r)
{
    DiagonalMatrix R(n_c,n_r);
    for(int i=0; i<n_c;i++)
    {
        R.Diag_Vals[i] = std::experimental::randint(-9,9);
    }
    return R;
}

Vectord r_Vec(int n)
{
    Vectord r(n);
    for(int i=0; i<n;i++)
    {
        r.Vec[i] = std::experimental::randint(-9,9);
    }
    return r;
}

BlockToeplitz r_BT(int n_r, int n_c, int prob)
{
    return BlockToeplitz(0,0,0);
}
