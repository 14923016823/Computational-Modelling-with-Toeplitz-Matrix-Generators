#include "Laplacian_imp.h"

#include <iostream>

int main() {
    int n = 10;
    int blockR = 100;
    int blockC = 100;

    const int numDiagsA = 3;
    int diagsA[numDiagsA]   = {-1, 0, 1};
    double valsA[numDiagsA] = {-1.0, 2.0, -1.0};

    SparseToeplitzd A(n, n, numDiagsA, diagsA, valsA);

    const int numDiagsB = 3;
    int diagsB[numDiagsB]   = {-3, 0, 3};
    double valsB[numDiagsB] = {5.0, -3.0, 5.0};

    SparseToeplitzd B(n, n, numDiagsB, diagsB, valsB);

    constexpr int TDiagN = 3;
    DiagEntry buildToeplitzList[TDiagN] = { {-1, &B}, {0, &A}, {1, &B} };

    BlockToeplitz T1(blockR, blockC, buildToeplitzList);

    DiagEntry buildToeplitzList2[TDiagN] = { {-2, &A}, {0, &B}, {2, &A} };

    BlockToeplitz T2(blockR, blockC, buildToeplitzList2);

    Vectord x(n);
    x.fill(1.0);

    std::cout << x.len() << '\n';
    std::cout << A.cols() << '\n';
    std::cout << B.cols() << '\n'; 


    LaplacianOp<Vectord, SparseToeplitzd, SparseToeplitzd> lapOp(A, B);
    Vectord result = lapOp * x;

    result.PrintVector();

    return 0;
}
