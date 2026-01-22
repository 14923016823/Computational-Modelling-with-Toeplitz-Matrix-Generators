//NOTE: THE FOLLOWING CODE IS A WORK IN PROGRESS -- IT MAY NOT RUN CORRECTLY AS IS.
//=========== Classes for Block Toeplitz matrix representation and operations ============
#include "Laplacian.h"

Laplacian2D_ToeplitzMatrix::Laplacian2D_ToeplitzMatrix(const int rows, const int cols, Vectord& b1) {
    // Generate the Laplacian matrix with non-homogeneous spacing and variable k
    generateMatrix(rows, cols, b1);
}
    


double Laplacian2D_ToeplitzMatrix::k_func(double x, double y) {
    // Example variable k function; modify as needed
    return 1.0 + 0.5 * (x + y);
}

void Laplacian2D_ToeplitzMatrix::generateMatrix(const int rows, const int cols, Vectord& b1) 
{
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
    //Incidence=&IncidenceMatrix;
    
    //     b.) Compute negative transpose
    Incidence=new BlockToeplitz(IncidenceMatrix);
    Matrix* negTransPtr = IncidenceMatrix.negativeTranspose();
    //BlockToeplitz negTrans=*negTrans;
    Incidence_T=(*negTransPtr).Clone(1);
    BlockToeplitz* negTransBlockPtr = static_cast<BlockToeplitz*>(negTransPtr); //cast to BlockToeplitz pointer for further operations

    //     c.) W_ee matrix generation function
    int dim = 2 * rows * cols;

    DiagonalMatrix W_ee_matrix(dim);

double dx = 1.0 / cols;
double dy = 1.0 / rows;

for (int i = 0; i < dim; ++i) 
{
    int col = i % cols;
    int row = i / cols;  
    double x = col * dx;
    double y = row * dy;
    double k_val = k_func(x, y);
    W_ee_matrix.Diag_Vals[i] = k_val/(dx*dy);  // use dx (or dy) for spacing
}

    Diagonal=W_ee_matrix.Clone(1.0);
}

Vectord Laplacian2D_ToeplitzMatrix::Laplacian(Vectord& input)
{
    Vectord b2 = (*Incidence) * input;
    Vectord Wb = (*Diagonal) * b2;
    Vectord final_b = (*Incidence_T) * Wb;
    return final_b;
}

Laplacian3D_ToeplitzMatrix::Laplacian3D_ToeplitzMatrix(const int rows, const int cols, const int arrays, Vectord& b1) {
    // Generate the Laplacian matrix with non-homogeneous spacing and variable k
    generateMatrix(rows, cols, arrays, b1);
}



double Laplacian3D_ToeplitzMatrix::k_func(double x, double y, double z) {
    // Example variable k function; modify as needed
    return 1.0 + (x + y + z)/3.0;
}

void Laplacian3D_ToeplitzMatrix::generateMatrix(const int rows, const int cols, const int arrays, Vectord& b1) 
{
    // Implementation for generating the Laplacian matrix with non-homogeneous mesh spacing and variable k


    // 1.)Upper incidence matrix generation (2D)
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


    // 2.) Lower incidence matrix generation (2D)
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

    // 3.) Lower incidence matrix generation (3D)
    //     a.) Array block matrix generation
    SparseToeplitz arrayBlockMatrix(rows*cols, rows*cols, 3); //leaf array block matrix with 3 diagonals
    columnBlockMatrix.Diags[0] = -(rows*cols - 1);
    columnBlockMatrix.Diags[1] = 0;
    columnBlockMatrix.Diags[2] = 1;
    columnBlockMatrix.Vals[0] = 1;
    columnBlockMatrix.Vals[1] = -1;
    columnBlockMatrix.Vals[2] = 1;

    //     b.) Identity matrix generation for lower incidence matrix
    SparseToeplitz identityMatrixArrays(arrays, arrays, 1); //leaf identity matrix
    identityMatrixArrays.Diags[0] = 0;
    identityMatrixArrays.Vals[0] = 1;

    //     c.) Column block matrix kronecker product with identity matrix (arrays)
    Matrix* lowerIncidenceMatrix_3D = arrayBlockMatrix.Kronecker(identityMatrixArrays); // lower incidence matrix
    BlockToeplitz* lowerIncidencePtr_3D = static_cast<BlockToeplitz*>(lowerIncidenceMatrix_3D);  //cast to BlockToeplitz pointer for further operations


    // 4.) Final Laplacian matrix assembly
    //     a.) Concatenate upper and lower incidence matrices
    int ndiags = 2; //number of diagonals in incidence matrix
    BlockToeplitz IncidenceMatrix_2D(rows * cols * 2, rows * cols, ndiags); // incidence matrix with 2 block rows
    IncidenceMatrix_2D.Diags[0] = 0; // diagonal for upper incidence matrix
    IncidenceMatrix_2D.Diags[1] = 1; // diagonal for lower incidence matrix
    IncidenceMatrix_2D.Vals[0] = rowBlockPtr;
    IncidenceMatrix_2D.Vals[1] = lowerIncidencePtr;
    //Incidence=&IncidenceMatrix;

    Matrix* upperIncidenceMatrix_3D = identityMatrixArrays.Kronecker(IncidenceMatrix_2D); // upper incidence matrix
    BlockToeplitz* BlockPtr_2D = static_cast<BlockToeplitz*>(upperIncidenceMatrix_3D);  //cast to BlockToeplitz pointer for further operations


    BlockToeplitz IncidenceMatrix_3D(rows * cols * arrays * 3, rows * cols * arrays, ndiags); // incidence matrix with 2 block rows
    IncidenceMatrix_3D.Diags[0] = 0; // diagonal for upper incidence matrix
    IncidenceMatrix_3D.Diags[1] = 1; // diagonal for lower incidence matrix
    IncidenceMatrix_3D.Vals[0] = BlockPtr_2D;
    IncidenceMatrix_3D.Vals[1] = lowerIncidencePtr_3D;
    
    //     b.) Compute negative transpose
    Incidence=new BlockToeplitz(IncidenceMatrix_3D);
    Matrix* negTransPtr = IncidenceMatrix_3D.negativeTranspose();
    //BlockToeplitz negTrans=*negTrans;
    Incidence_T=(*negTransPtr).Clone(1);
    BlockToeplitz* negTransBlockPtr = static_cast<BlockToeplitz*>(negTransPtr); //cast to BlockToeplitz pointer for further operations

    //     c.) W_ee matrix generation function
    int dim = 3 * rows * cols * arrays;

    DiagonalMatrix W_ee_matrix(dim);

double dx = 1.0 / cols;
double dy = 1.0 / rows;
double dz = 1.0 / arrays;

for (int i = 0; i < dim; ++i) 
{
    int col = i % cols;
    int row = i / cols;  
    int arr = i;
    double x = col * dx;
    double y = row * dy;
    double z = arr * dz;
    double k_val = k_func(x, y, z);
    W_ee_matrix.Diag_Vals[i] = k_val/(dx*dy*dz);  // use dx (or dy) for spacing
}

    Diagonal=W_ee_matrix.Clone(1.0);
}

Vectord Laplacian3D_ToeplitzMatrix::Laplacian(Vectord& input)
{
    Vectord b2 = (*Incidence) * input;
    Vectord Wb = (*Diagonal) * b2;
    Vectord final_b = (*Incidence_T) * Wb;
    return final_b;
}