/* Laplacian matrix generation based on a looping stencil - 2D case */
#include <iostream>
#include <vector>
#include <iomanip>
#include <stdexcept>

#include "SparseToeplitz.h"
#include "BlockToeplitz.h"
#include "DiagonalMatrix.h"
#include "Vectord.h"

//======================================== TO DOs ========================================
// - PRIORITY: Integration with new BlockToeplitz class structure -- Done
// - Implement "variable" k
// - Implment non-homogeneous mesh spacing = no need for this
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
    Laplacian2D_ToeplitzMatrix(const int rows, const int cols, Vectord& b1) {
        // Generate the Laplacian matrix with non-homogeneous spacing and variable k
        generateMatrix(rows, cols, b1);
    }
    
    private:
    double w(double k, double a) {
        return k / (a * a);
    }

    double k_func(double x, double y) {
        // Example variable k function; modify as needed
        return 1.0 + 0.5 * (x + y);
    }

    Vectord generateMatrix(const int rows, const int cols, Vectord& b1) {
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
        int dim = rows * cols;
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
};


//=============================== Main function ========================================
int main() {
    int rows = 3;
    int cols = 3;

    Vectord vec(rows * cols);
    for (int i = 0; i < rows * cols; ++i) {
        vec.Vec[i] = 1.0; // Example initialization
    }
    
    //Laplacian2D_FullMatrix laplacian(rows, cols);
    Laplacian2D_ToeplitzMatrix laplacian_toeplitz(rows, cols, vec);

    return 0;
}