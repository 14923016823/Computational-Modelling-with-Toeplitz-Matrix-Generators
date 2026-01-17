#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <memory>
#include <tuple>

#include "Vectord.h"
#include "SparseToeplitz.h"
#include "BlockToeplitz.h"
#include "Runge_Kutta4.h"
#include "RHS.h"

using namespace std::chrono;

struct Laplacian {
    // Start of Laplacian parameters
    int nx;
    int ny;
    double dx; 
    double dy;
    double alpha; // thermal diffusivity
    std::unique_ptr<std::tuple<BlockToeplitzd, BlockToeplitzd>> L;
    // End of Laplacian parameters

    // Constructor
    Laplacian(int nx_, int ny_, double dx_, double dy_, double alpha_)
        : nx(nx_), ny(ny_), dx(dx_), dy(dy_), alpha(alpha_) 
        {
            std::vector<int> diagsX = {-nx+1, -1, 0, 1, nx-1};
            std::vector<double> valsX = {-1.0, -1.0, 2.0, -1.0, -1.0};
            SparseToeplitzd A(nx, nx, diagsX.size(),
                              diagsX.data(), valsX.data());

            std::vector<int> diagsY = {-ny+1, -1, 0, 1, ny-1};
            std::vector<double> valsY = {-1.0, -1.0, 2.0, -1.0, -1.0};
            SparseToeplitzd B(ny, ny, diagsY.size(),
                              diagsY.data(), valsY.data());

            // Create BlockToeplitz matrices
            BlockToeplitzd Lx(nx*ny, nx*ny, diagsX.size());
            BlockToeplitzd Ly(nx*ny, nx*ny, diagsY.size());

            for (int i = 0; i < diagsX.size(); i++) {
                Lx.Diags[i] = diagsX[i];
                Ly.Diags[i] = diagsY[i];
                Lx.Vals[i] = A.Clone(static_cast<double>(valsX[i]));
                Ly.Vals[i] = B.Clone(static_cast<double>(valsY[i]));
            }

            L = std::make_unique<std::tuple<BlockToeplitzd, BlockToeplitzd>>(std::move(Lx), std::move(Ly));
        }
};

struct HeatEquationRHS: public RHS {
    private:
    
    Laplacian& laplacian;
    void operator()(double t, const Vectord& u, Vectord& rhs) override {
        const auto& [Lx, Ly] = *(laplacian.L);
        Vectord temp1(u.len());
        Vectord temp2(u.len());

        Lx.matvec(u, temp1);
        Ly.matvec(u, temp2);

        rhs = temp1;
        rhs.axpy(1.0, temp2);
        rhs.scal(laplacian.alpha / (laplacian.dx * laplacian.dy));
    }
};

int main(int argc, char* argv[])
{
   
    return 0;
}