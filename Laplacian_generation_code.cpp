/* Laplacian matrix generation based on a looping stencil - 2D case */
#include <iostream>
#include <vector>
#include <iomanip>
#include <stdexcept>

#include "SparseToeplitz.h"
#include "BlockToeplitz.h"

//======================================== TO DOs ========================================
// - PRIORITY: Integration with new BlockToeplitz class structure -- Done
// - Implement "variable" k
// - Implment non-homogeneous mesh spacing
// - Scale up to 3D case?


//================ Classes for full matrix representation and operations =================
/*
class Matrix {
public:
    // Constructors
    Matrix() : rows(0), cols(0), data() {
        // Default constructor
    }

    Matrix(const int rows, const int cols) : rows(rows), cols(cols), data(rows,
            std::vector<double>(cols, 0.0)) {    
        // Initialize matrix with zeros
    }

    // Access operators
    double &operator()(const int i, const int j) {
        // Access element at (i, j)
        return data[i][j];
    }

    const double &operator()(const int i, const int j) const {
        // Access element at (i, j) (const version)
        return data[i][j];
    }

    int numRows() const {
        // Get number of rows from private member
        return rows;
    }

    int numCols() const {
        // Get number of columns from private member
        return cols;
    }

    // Functions
    Matrix kroneckerProduct(const Matrix& other) const {
        // Implement Kronecker product logic here
        Matrix result(this->numRows() * other.numRows(), this->numCols() * other.numCols());
    
        for (int i = 0; i < this->numRows(); ++i) {
            for (int j = 0; j < this->numCols(); ++j) {
                for (int k = 0; k < other.numRows(); ++k) {
                    for (int l = 0; l < other.numCols(); ++l) {
                        result.data[i * other.numRows() + k][j * other.numCols() + l] =
                            this->data[i][j] * other.data[k][l];
                    }
                }
            }
        }
        return result;
    }
    
    private:
    int rows;
    int cols;
    std::vector<std::vector<double>> data;
};

class Laplacian2D_FullMatrix {
public:
    Laplacian2D_FullMatrix(const int rows, const int cols) {
        // Generate the Laplacian matrix
        generateMatrix(rows, cols);
    }

private:
    void generateMatrix(const int rows, const int cols) {
        // Implementation for generating the Laplacian matrix
        
        // 1.)Upper incidence matrix generation
        //     a.) Row block matrix generation
        Matrix rowBlockMatrix(cols, cols);
        
        for (int i = 0; i < cols; ++i) {
            for (int j = 0; j < cols; ++j) {
                rowBlockMatrix(i, j) = (j == (i + 1) % cols) ? 1 : 0;
                rowBlockMatrix(i, j) = (j == i) ? -1 : rowBlockMatrix(i, j);
            }
        }

        //     b.) Identity matrices generation
        Matrix identityMatrixRows(rows, rows);
        Matrix identityMatrixCols(cols, cols);

        for (int i = 0; i < rows; ++i)
            identityMatrixRows(i, i) = 1;

        for (int i = 0; i < cols; ++i)
            identityMatrixCols(i, i) = 1;

        //     c.) Identity matrix (rows) kronecker product with row block matrix
        Matrix upperIncidenceMatrix = identityMatrixRows.kroneckerProduct(rowBlockMatrix);

        for (int i = 0; i < upperIncidenceMatrix.numRows(); ++i) {
            for (int j = 0; j < upperIncidenceMatrix.numCols(); ++j) {
            }
        }
    
        // 2.) Lower incidence matrix generation
        //     a.) Column block matrix generation
        Matrix columnBlockMatrix(rows, rows);

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < rows; ++j) {
                columnBlockMatrix(i, j) = (j == (i + 1) % rows) ? 1 : 0;
                columnBlockMatrix(i, j) = (j == i) ? -1 : columnBlockMatrix(i, j);
            }
        }

        //     b.) Column block matrix kronecker product with identity matrix (cols)
        Matrix lowerIncidenceMatrix = columnBlockMatrix.kroneckerProduct(identityMatrixCols);

        // 3.) Final Laplacian matrix assembly
        //     a.) Concatenate upper and lower incidence matrices
        Matrix concatenatedMatrix(upperIncidenceMatrix.numRows() + lowerIncidenceMatrix.numRows(),
                                   upperIncidenceMatrix.numCols());
        
        for (int i = 0; i < upperIncidenceMatrix.numRows(); ++i) {
            for (int j = 0; j < upperIncidenceMatrix.numCols(); ++j) {
                concatenatedMatrix(i, j) = upperIncidenceMatrix(i, j);
            }
        }
        for (int i = 0; i < lowerIncidenceMatrix.numRows(); ++i) {
            for (int j = 0; j < lowerIncidenceMatrix.numCols(); ++j) {
                concatenatedMatrix(i + upperIncidenceMatrix.numRows(), j) = lowerIncidenceMatrix(i, j); 
            }
        }
        
        //     b.) Compute negative transpose of the concatenated matrix
        Matrix negativeTranspose(concatenatedMatrix.numCols(), concatenatedMatrix.numRows());
        for (int i = 0; i < concatenatedMatrix.numRows(); ++i) {
            for (int j = 0; j < concatenatedMatrix.numCols(); ++j) {
                negativeTranspose(j, i) = -concatenatedMatrix(i, j);
            }
        }
        //     c.) Multiply incidence matrix with its negative transpose to get Laplacian
        Matrix laplacianMatrix(concatenatedMatrix.numCols(), concatenatedMatrix.numCols());
        for (int i = 0; i < laplacianMatrix.numRows(); ++i) {
            for (int j = 0; j < laplacianMatrix.numCols(); ++j) {
                double sum = 0.0;
                for (int k = 0; k < concatenatedMatrix.numRows(); ++k) {
                    sum += concatenatedMatrix(k, i) * negativeTranspose(j, k);
                }
                laplacianMatrix(i, j) = sum;
            }
        }
        // Print the Laplacian matrix
        std::cout << "Laplacian Matrix (" << laplacianMatrix.numRows() << "x" <<
            laplacianMatrix.numCols() << "):\n";
        for (int i = 0; i < laplacianMatrix.numRows(); ++i) {
            for (int j = 0; j < laplacianMatrix.numCols(); ++j) {
                std::cout << laplacianMatrix(i, j) << " ";
            }
            std::cout << "\n";
        }

    }
};
*/


//NOTE: THE FOLLOWING CODE IS A WORK IN PROGRESS -- IT MAY NOT RUN CORRECTLY AS IS.
//=========== Classes for Block Toeplitz matrix representation and operations ============

class Laplacian2D_ToeplitzMatrix {
    public:
    Laplacian2D_ToeplitzMatrix(const int rows, const int cols) {
        // Generate the Laplacian matrix (homogeneous mesh / default k)
        generateMatrix(rows, cols);
    }

    Laplacian2D_ToeplitzMatrix(const int rows, const int cols, const double k_ex, const double k_ey, const double dx, const double dy) {
        // Generate the Laplacian matrix with non-homogeneous spacing and variable k
        generateMatrix(rows, cols, k_ex, k_ey, dx, dy);
    }
    
    private:
    double w(double k, double a) {
        return k / (a * a);
    }

    SparseToeplitz W_e_rows(const int rows, const double k_ex, const double dx, const double dy) {
        // Generate W_ee block matrix for row connections; "multiplies" the negative transpose row block amtrix
        SparseToeplitz W_ee_rows(rows, rows, 1);
        W_ee_rows.Diags[0] = 0;
        W_ee_rows.Vals[0] = w(k_ex, dx);

        return W_ee_rows;
    }

    SparseToeplitz W_e_cols(const int rows, const int cols, const double k_ey, const double dx, const double dy) {
        // Generate W_ee matrix for column connections; "multiplies" the negative transpose lower (column) incidence matrix
        int dim = rows * cols;

        SparseToeplitz W_ee_cols(dim, dim, 1);
        W_ee_cols.Diags[0] = 0;
        W_ee_cols.Vals[0] = w(k_ey, dy);

        return W_ee_cols;
    }

    void generateMatrix(const int rows, const int cols) {
        // Implementation for generating the Laplacian matrix

        // The "matrix multiplication" between the negative transpose of the incidence matrix and the incidence matrix itself is split between the
        // contrubution of the upper part (desribing the influence of the row connections) and the lower part (describing the influence of the
        // column connections). This is done to exploit the Kronecker product structure of the incidence matrices and avoid forming the full
        // incidence matrix.


        // 1.)Upper incidence matrix generation
        //     a.) Row block matrix generation
        SparseToeplitz rowBlockMatrix = SparseToeplitz(cols, cols, 3); //leaf row block matrix with 3 diagonals
        rowBlockMatrix.Diags[0] = -(cols - 1);
        rowBlockMatrix.Diags[1] = 0;
        rowBlockMatrix.Diags[2] = 1;
        rowBlockMatrix.Vals[0] = 1;
        rowBlockMatrix.Vals[1] = -1;
        rowBlockMatrix.Vals[2] = 1;

        //     b.) Identity matrices generation for upper incidence matrix
        SparseToeplitz identityMatrixRows = SparseToeplitz(rows, rows, 1); //leaf identity matrix
        identityMatrixRows.Diags[0] = 0;
        identityMatrixRows.Vals[0] = 1;


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


        // 3.) Final Laplacian matrix assembly
        //     a.) Concatenate upper and lower incidence matrices
        //          --not required--
        
        //     b.) Compute negative transposes
        Matrix* negTransposeUpperBlock = rowBlockMatrix.negativeTranspose(); // negative transpose of the row block matrix
        SparseToeplitz* negPtr = static_cast<SparseToeplitz*>(negTransposeUpperBlock); //cast to SparseToeplitz pointer for further operations
        
        Matrix* negTransposeColumnBlockMatrix = columnBlockMatrix.negativeTranspose(); // negative transpose of the column block matrix
        Matrix* lowerIncidenceMatrixTranspose = (*negTransposeColumnBlockMatrix).Kronecker(identityMatrixCols); // negative transpose of the lower incidence matrix
        BlockToeplitz* lowerPtr = static_cast<BlockToeplitz*>(lowerIncidenceMatrixTranspose); //cast to BlockToeplitz pointer for further operations

        //     c.) Multiply incidence matrix with its negative transpose to get Laplacian
        SparseToeplitz upperIncidenceBlock = SparseToeplitz(rows, rows, 5); // upper incidence block matrix with 5 diagonals
        upperIncidenceBlock.Diags[0] = -(rows - 1);
        upperIncidenceBlock.Diags[1] = -1;
        upperIncidenceBlock.Diags[2] = 0;
        upperIncidenceBlock.Diags[3] = 1;
        upperIncidenceBlock.Diags[4] = (rows - 1);
        upperIncidenceBlock.Vals[0] = negPtr->Vals[1] * columnBlockMatrix.Vals[2];
        upperIncidenceBlock.Vals[1] = negPtr->Vals[2] * columnBlockMatrix.Vals[1];
        upperIncidenceBlock.Vals[2] = negPtr->Vals[1] * columnBlockMatrix.Vals[1] + negPtr->Vals[2] * columnBlockMatrix.Vals[0];
        upperIncidenceBlock.Vals[3] = negPtr->Vals[0] * columnBlockMatrix.Vals[1];
        upperIncidenceBlock.Vals[4] = negPtr->Vals[0] * columnBlockMatrix.Vals[1];

        Matrix* upperIncidenceLaplacian = identityMatrixRows.Kronecker(upperIncidenceBlock); // upper incidence Laplacian matrix
        BlockToeplitz* upperIncidencePtr = static_cast<BlockToeplitz*>(upperIncidenceLaplacian); //cast to BlockToeplitz pointer for further operations

        SparseToeplitz lowerIncidenceLaplacian = SparseToeplitz(rows, rows, 1); // lower incidence Laplacian block matrix with 1 diagonal
        lowerIncidenceLaplacian.Diags[0] = 0;
        { // computing product of lower incidence matrix with its negative transpose
            BlockToeplitz* L = lowerPtr;
            BlockToeplitz* LI = lowerIncidencePtr;
            SparseToeplitz* A = dynamic_cast<SparseToeplitz*>(L->Vals[1]);
            SparseToeplitz* B = dynamic_cast<SparseToeplitz*>(LI->Vals[1]);
            SparseToeplitz* C = dynamic_cast<SparseToeplitz*>(L->Vals[2]);
            SparseToeplitz* D = dynamic_cast<SparseToeplitz*>(LI->Vals[0]);
            if (!A || !B || !C || !D) {
                throw std::logic_error("expected SparseToeplitz leaf blocks in lower incidence matrices");
            }
            // Use index 0 for leaf blocks (they each have a single diagonal stored at index 0)
            lowerIncidenceLaplacian.Vals[0] = A->Vals[0] * B->Vals[0] + C->Vals[0] * D->Vals[0];
        }

        Matrix* laplacianMatrix = (*upperIncidencePtr).Add(lowerIncidenceLaplacian);
        BlockToeplitz& laplacianPtr = *static_cast<BlockToeplitz*>(laplacianMatrix);
        // resulting matrix is the Laplacian matrix with Dirichlet boundary conditions
        laplacianPtr.printFullMatrix();
            
        //     d.) Store as RecursiveToeplitz matrix

    }

    void generateMatrix(const int rows, const int cols, const double k_ex, const double k_ey, const double dx, const double dy) {
        // Implementation for generating the Laplacian matrix with non-homogeneous mesh spacing and variable k

        // The "matrix multiplication" between the negative transpose of the incidence matrix and the incidence matrix itself is split between the
        // contrubution of the upper part (desribing the influence of the row connections) and the lower part (describing the influence of the
        // column connections). This is done to exploit the Kronecker product structure of the incidence matrices and avoid forming the full
        // incidence matrix.


        // 1.)Upper incidence matrix generation
        //     a.) Row block matrix generation
        SparseToeplitz rowBlockMatrix = SparseToeplitz(cols, cols, 3); //leaf row block matrix with 3 diagonals
        rowBlockMatrix.Diags[0] = -(cols - 1);
        rowBlockMatrix.Diags[1] = 0;
        rowBlockMatrix.Diags[2] = 1;
        rowBlockMatrix.Vals[0] = 1;
        rowBlockMatrix.Vals[1] = -1;
        rowBlockMatrix.Vals[2] = 1;

        //     b.) Identity matrices generation for upper incidence matrix
        SparseToeplitz identityMatrixRows = SparseToeplitz(rows, rows, 1); //leaf identity matrix
        identityMatrixRows.Diags[0] = 0;
        identityMatrixRows.Vals[0] = 1;


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


        // 3.) Final Laplacian matrix assembly
        //     a.) Concatenate upper and lower incidence matrices
        //          --not required--
        
        //     b.) Compute negative transposes
        Matrix* negTransposeUpperBlock = rowBlockMatrix.negativeTranspose(); // negative transpose of the row block matrix
        SparseToeplitz* negPtr = static_cast<SparseToeplitz*>(negTransposeUpperBlock); //cast to SparseToeplitz pointer for further operations
        
        Matrix* negTransposeColumnBlockMatrix = columnBlockMatrix.negativeTranspose(); // negative transpose of the column block matrix
        Matrix* lowerIncidenceMatrixTranspose = (*negTransposeColumnBlockMatrix).Kronecker(identityMatrixCols); // negative transpose of the lower incidence matrix
        BlockToeplitz* lowerPtr = static_cast<BlockToeplitz*>(lowerIncidenceMatrixTranspose); //cast to BlockToeplitz pointer for further operations

        //     c.) Comput W_ee matrices
        SparseToeplitz W_ee_rows_matrix = W_e_rows(rows, k_ex, dx, dy);
        SparseToeplitz W_ee_cols_matrix = W_e_cols(rows, cols, k_ey, dx, dy);

        //     d.) Multiply incidence matrix with its negative transpose to get Laplacian
        SparseToeplitz upperIncidenceBlock = SparseToeplitz(rows, rows, 5); // upper incidence block matrix with 5 diagonals
        upperIncidenceBlock.Diags[0] = -(rows - 1);
        upperIncidenceBlock.Diags[1] = -1;
        upperIncidenceBlock.Diags[2] = 0;
        upperIncidenceBlock.Diags[3] = 1;
        upperIncidenceBlock.Diags[4] = (rows - 1);
        upperIncidenceBlock.Vals[0] = negPtr->Vals[1] * columnBlockMatrix.Vals[2];
        upperIncidenceBlock.Vals[1] = negPtr->Vals[2] * columnBlockMatrix.Vals[1];
        upperIncidenceBlock.Vals[2] = negPtr->Vals[1] * columnBlockMatrix.Vals[1] + negPtr->Vals[2] * columnBlockMatrix.Vals[0];
        upperIncidenceBlock.Vals[3] = negPtr->Vals[0] * columnBlockMatrix.Vals[1];
        upperIncidenceBlock.Vals[4] = negPtr->Vals[0] * columnBlockMatrix.Vals[1];

        for (int i = 0; i < 5; ++i)
            upperIncidenceBlock.Vals[i] *= W_ee_rows_matrix.Vals[0];

        Matrix* upperIncidenceLaplacian = identityMatrixRows.Kronecker(upperIncidenceBlock); // upper incidence Laplacian matrix
        BlockToeplitz* upperIncidencePtr = static_cast<BlockToeplitz*>(upperIncidenceLaplacian); //cast to BlockToeplitz pointer for further operations

        SparseToeplitz lowerIncidenceLaplacian = SparseToeplitz(rows, rows, 1); // lower incidence Laplacian block matrix with 1 diagonal
        lowerIncidenceLaplacian.Diags[0] = 0;
        { // computing product of lower incidence matrix with its negative transpose
            BlockToeplitz* L = lowerPtr;
            BlockToeplitz* LI = lowerIncidencePtr;
            SparseToeplitz* A = dynamic_cast<SparseToeplitz*>(L->Vals[1]);
            SparseToeplitz* B = dynamic_cast<SparseToeplitz*>(LI->Vals[1]);
            SparseToeplitz* C = dynamic_cast<SparseToeplitz*>(L->Vals[2]);
            SparseToeplitz* D = dynamic_cast<SparseToeplitz*>(LI->Vals[0]);
            if (!A || !B || !C || !D) {
                throw std::logic_error("expected SparseToeplitz leaf blocks in lower incidence matrices");
            }
            // Use index 0 for leaf blocks (they each have a single diagonal stored at index 0)
            lowerIncidenceLaplacian.Vals[0] = (A->Vals[0] * B->Vals[0] + C->Vals[0] * D->Vals[0]) * W_ee_cols_matrix.Vals[0];
        }

        Matrix* laplacianMatrix = (*upperIncidencePtr).Add(lowerIncidenceLaplacian);
        BlockToeplitz& laplacianPtr = *static_cast<BlockToeplitz*>(laplacianMatrix);
        // resulting matrix is the Laplacian matrix with Dirichlet boundary conditions
        laplacianPtr.printFullMatrix();
            
        //     d.) Store as RecursiveToeplitz matrix

    }
};


//=============================== Main function ========================================
int main() {
    int rows = 4;
    int cols = 4;
    double k_ex = 2.0;
    double k_ey = 3.0;
    double dx = 0.5;
    double dy = 0.5;
    
    //Laplacian2D_FullMatrix laplacian(rows, cols);
    //Laplacian2D_ToeplitzMatrix laplacian_toeplitz(rows, cols);
    Laplacian2D_ToeplitzMatrix laplacian_toeplitz(rows, cols, k_ex, k_ey, dx, dy);

    return 0;
}