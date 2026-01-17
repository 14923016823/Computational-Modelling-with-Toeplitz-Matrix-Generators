#include <vector>

#include "SparseToeplitz.h"

int main() {
    int nx = 10;
    int ny = 8;
    int Num_Diags = 5;
    std::vector<int> diags = {-9,-1, 0, 1, 7};
    std::vector<double> vals = {-1.0, -1.0, 2.0, -1.0, -1.0};

    SparseToeplitz<double> A(nx, ny, Num_Diags, diags.data(), vals.data());

    std::cout << "SparseToeplitz Matrix A (" << nx << " x " << ny << "):\n";
    A.print();

    Vectord x(ny);
    for (int i = 0; i < ny; ++i) {
        x[i] = static_cast<double>(i + 1);
    }

        for (int i = 0; i < ny; ++i) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }

    Vectord y(nx);

    A.matvec(x, y);
    std::cout << "Result of SparseToeplitz matvec:\n";
    for (int i = 0; i < nx; ++i) {
        std::cout << "y[" << i << "] = " << y[i] << std::endl;
    }   
    return 0;
}