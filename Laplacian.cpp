//NOTE: THE FOLLOWING CODE IS A WORK IN PROGRESS -- IT MAY NOT RUN CORRECTLY AS IS.
//=========== Classes for Block Toeplitz matrix representation and operations ============
#include "Laplacian.h"

Laplacian2D_ToeplitzMatrix::Laplacian2D_ToeplitzMatrix(const int rows, const int cols, Vectord& b1) {
    // Generate the Laplacian matrix with non-homogeneous spacing and variable k
    generateMatrix(rows, cols, b1);
}
    
double Laplacian2D_ToeplitzMatrix::w(double k, double a) {
    return k / (a * a);
}

double Laplacian2D_ToeplitzMatrix::k_func(double x, double y) {
    // Example variable k function; modify as needed
    return 1.0 + 0.5 * (x + y);
}

Vectord Laplacian2D_ToeplitzMatrix::generateMatrix(const int rows, const int cols, Vectord& b1) {
    // Implementation for generating the Laplacian matrix with non-homogeneous mesh spacing and variable k


    // 1.)Upper incidence matrix generation
    //     a.) Row block matrix generation
    SparseToeplitz rowBlockMatrix = SparseToeplitz(cols, cols, 3); //leaf row block matrix with 3 diagonals
    rowBlockMatrix.Diags[0] = -(cols - 1);
    rowBlockMatrix.Diags[1] = 0;
    rowBlockMatrix.Diags[2] = 1;
    rowBlockMatrix.Vals[0] = 1;
    rowBlockMatrix.Vals[1] = -1;
    rowBlockMatrix.Vals[2] = 1;

    //     b.) Identity matrix generation for upper incidence matrix
    SparseToeplitz identityMatrixRows = SparseToeplitz(rows, rows, 1); //leaf identity matrix
    identityMatrixRows.Diags[0] = 0;
    identityMatrixRows.Vals[0] = 1;

    //     c.) Identity matrix (rows) kronecker product with row block matrix
    Matrix* upperIncidenceMatrix = identityMatrixRows.Kronecker(rowBlockMatrix); // upper incidence matrix
    BlockToeplitz* rowBlockPtr = static_cast<BlockToeplitz*>(upperIncidenceMatrix);  //cast to BlockToeplitz pointer for further operations


    // 2.) Lower incidence matrix generation
    //     a.) Column block matrix generation
    SparseToeplitz columnBlockMatrix(rows, rows, 3); //leaf column block matrix with 3 diagonals
    columnBlockMatrix.Diags[0] = -(rows - 1);
    columnBlockMatrix.Diags[1] = 0;
    columnBlockMatrix.Diags[2] = 1;
    columnBlockMatrix.Vals[0] = 1;
    columnBlockMatrix.Vals[1] = -1;
    columnBlockMatrix.Vals[2] = 1;

    //     b.) Identity matrix generation for lower incidence matrix
    SparseToeplitz identityMatrixCols(cols, cols, 1); //leaf identity matrix
    identityMatrixCols.Diags[0] = 0;
    identityMatrixCols.Vals[0] = 1;

    //     c.) Column block matrix kronecker product with identity matrix (cols)
    Matrix* lowerIncidenceMatrix = columnBlockMatrix.Kronecker(identityMatrixCols); // lower incidence matrix
    BlockToeplitz* lowerIncidencePtr = static_cast<BlockToeplitz*>(lowerIncidenceMatrix);  //cast to BlockToeplitz pointer for further operations


    // 4.) Final Laplacian matrix assembly
    //     a.) Concatenate upper and lower incidence matrices
    int ndiags = 2; //number of diagonals in incidence matrix
    BlockToeplitz IncidenceMatrix(rows * cols * 2, rows * cols, ndiags); // incidence matrix with 2 block rows
    IncidenceMatrix.Diags[0] = 0; // diagonal for upper incidence matrix
    IncidenceMatrix.Diags[1] = 1; // diagonal for lower incidence matrix
    IncidenceMatrix.Vals[0] = rowBlockPtr;
    IncidenceMatrix.Vals[1] = lowerIncidencePtr;
    
    //     b.) Compute negative transpose
    Matrix* negTransPtr = IncidenceMatrix.negativeTranspose();
    BlockToeplitz* negTransBlockPtr = static_cast<BlockToeplitz*>(negTransPtr); //cast to BlockToeplitz pointer for further operations

    //     c.) W_ee matrix generation function
    int dim = 2 * rows * cols;
    DiagonalMatrix W_ee_matrix(dim);

    for (int i = 0; i < dim; ++i) {
        W_ee_matrix.Diag_Vals[i] = k_func(i % cols, i / rows);
    }

    //     d.) Matrix-vector multiplcation setup to obtain solution to Laplacian system
    Vectord b2 = IncidenceMatrix * b1;
    Vectord Wb = W_ee_matrix * b2;
    Vectord final_b = (*negTransBlockPtr) * Wb;

    return final_b;
}