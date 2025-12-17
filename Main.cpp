#include <iostream>
#include "SparseToeplitz.h"
#include "ConjugateGradient.h"

#include <vector>

static double residual_norm(Matrix& A, const Vectord& b, Vectord& x)
{
    Vectord Ax = A * x;       // uses your Matrix::operator*(Vectord&)
    Vectord r  = b;
    r.axpy(-1.0, Ax);         // r = b - Ax
    return r.norm();
}

int main(int argc, char* argv[])
{
    // Problem size (can be set via command-line argument)
    int n = 200;
    if (argc > 1) {
        int arg_n = std::atoi(argv[1]);
        if (arg_n > 0) {
            n = arg_n;
        } else {
            std::cerr << "Invalid problem size argument. Using default n = 200.\n";
        }
    }

    // Build A: 1D Laplacian Toeplitz (Dirichlet-style interior operator)
    // diagonals: -1, 0, +1 with values: -1, 2, -1
    const int numDiags = 3;
    int diags[numDiags]   = {-1, 0, 1};
    double vals[numDiags] = {-1.0, 2.0, -1.0};

    SparseToeplitz A(n, n, numDiags, diags, vals);
    A.print();

    //Vectord b(n);
    //b.fill(5.0);
    //b.PrintVector();

    std::vector<double> b_data(n);
    for (int i = 0; i < n; ++i)
        b_data[i] = double(i*i + 3*i - 2);

    Vectord b(b_data.data(), n);

    b.PrintVector();


    // x0 = zeros
    Vectord x(n);
    x.fill(0.0);

    // CG parameters
    const int maxIters = 5000;
    const double relTol = 1e-6;
    const double absTol = 1e-6;

    std::cout << "Running CG on 1D Laplacian Toeplitz, n=" << n << "\n";

    double r0 = residual_norm(A, b, x);
    std::cout << "Initial residual norm ||r0|| = " << r0 << "\n";

    int iters = ConjugateGradient(A, b, x, maxIters, relTol, absTol);

    if (iters < 0)
    {
        std::cout << "CG failed (dimension mismatch or numerical breakdown).\n";
        return 1;
    }

    double rf = residual_norm(A, b, x);
    std::cout << "CG iterations = " << iters << "\n";
    std::cout << "Final residual norm ||r|| = " << rf << "\n";
    std::cout << "Relative residual ||r||/||r0|| = " << (rf / (r0 > 0 ? r0 : 1.0)) << "\n";

    x.PrintVector();

    return 0;
}
