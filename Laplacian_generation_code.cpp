/* Laplacian matrix generation based on a looping stencil - 2D case */
#include <iostream>
#include <vector>
#include <iomanip>
#include <stdexcept>

#include "SparseToeplitz.h"
#include "BlockToeplitz.h"

//======================================== TO DOs ========================================
// - PRIORITY: Integration with new BlockToeplitz class structure
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
        // Generate the Laplacian matrix
        generateMatrix(rows, cols);
        //RecursiveToeplitz LaplacianMatrix = RecursiveToeplitz(2, 1, 1, 1, 2); //initialize empty matrix
    }
    
    //Function to pring the laplacian matrix
    /*void PrintLaplacianMatrix(const RecursiveToeplitz& laplacian) {
    
    }*/
    
    private:
    void generateMatrix(const int rows, const int cols) {
        // Implementation for generating the Laplacian matrix

        // 1.)Upper incidence matrix generation
        //     a.) Row block matrix generation
        SparseToeplitz rowBlockMatrix = SparseToeplitz(cols, cols, 3); //leaf matrix with 3 diagonals
        rowBlockMatrix.Diags[0] = -(cols - 1);
        rowBlockMatrix.Diags[1] = 0;
        rowBlockMatrix.Diags[2] = 1;
        rowBlockMatrix.Vals[0] = 1;
        rowBlockMatrix.Vals[1] = -1;
        rowBlockMatrix.Vals[2] = 1;

        //     b.) Identity matrix generation for lower incidence matrix
        SparseToeplitz identityMatrixCols(cols, cols, 1); //leaf identity matrix
        identityMatrixCols.Diags[0] = 0;
        identityMatrixCols.Vals[0] = 1;

        // 2.) Lower incidence matrix generation
        //     a.) Column block matrix generation
        SparseToeplitz columnBlockMatrix(rows, rows, 3); //leaf matrix with 3 diagonals
        columnBlockMatrix.Diags[0] = -(rows - 1);
        columnBlockMatrix.Diags[1] = 0;
        columnBlockMatrix.Diags[2] = 1;
        columnBlockMatrix.Vals[0] = 1;
        columnBlockMatrix.Vals[1] = -1;
        columnBlockMatrix.Vals[2] = 1;

        //     b.) Column block matrix kronecker product with identity matrix (cols)
        Matrix* lowerIncidenceMatrix = columnBlockMatrix.Kronecker(identityMatrixCols);
        Matrix* columnBlockMatrixTranspose = columnBlockMatrix.negativeTranspose();
        Matrix* lowerIncidenceMatrixTranspose = (*columnBlockMatrixTranspose).Kronecker(identityMatrixCols);

        // 3.) Final Laplacian matrix assembly
        //     a.) Concatenate upper and lower incidence matrices
        
        
        //     b.) Compute negative transposes
        

        (*lowerIncidenceMatrix).printFullMatrix();

        //negTransposeUpperBlock.PrintFullToeplitz();

        //     c.) Multiply incidence matrix with its negative transpose to get Laplacian
            SparseToeplitz upperIncidenceBlock = SparseToeplitz(rows, rows, 3);

            // For "upper incidence matrix" per block, matrix product:
            // cols = number of columns in original mesh grid
            //  - Main diagonal:
            // result[i,i] = blockT[i,i]*block[i,i] + blockT[i+1,i]*block[i,i+1] for i in [0,cols-2]
            // result[cols-1,cols-1] = blockT[cols-1,cols-1]*block[cols-1,0] + block[cols-1,cols-1]*blockT[cols-1,cols-1]

            /*for (int i=0; i<cols; i++) {
                upperIncidenceBlock.Vals[0] = negTransposeUpperBlock.Vals[i]*rowBlockMatrix.Vals[i] + negTransposeUpperBlock[];*/
            

            
            //  - First upper diagonal:
            // result[i,i+1] = blockT[i,i]*block[i,i+1] for i in [0,cols-2]

            //  - Last upper diagonal:
            // result[0,cols-1] = blockT[0,cols-1]*block[cols-1,cols-1]

            //  - First lower diagonal:
            // result[i+1,i] = blockT[i+1,i]*block[i,i] for i in [0,cols-2]

            // - Last lower diagonal:
            // result[cols-1,0] = blockT[cols-1,cols-1]*block[0,cols-1]

                // The result matrix here will be of size cols x cols, and can be kronecker-multiplied with
                // identity matrix (rows) to get the full upper incidence contribution to the Laplacian.


            //For "lower incidence matrix":
            // rows = number of rows in original mesh grid
            // length = rows*rows
            // - Main diagonal:
            // result[i,i] = blockT[i,i]*block[i,i] + blockT[i,i+2*cols]*block[i+2*cols,i] for i in [0,cols-1]
            // result[i,i] = blockT[i,i-cols]*block[i-cols,i] + blockT[i,i]*block[i,i] for i in [cols,length-1]

            // - First upper diagonal:
            // result[i,i+rows] = block[i,i+rows]*blockT[i,i] for i in [0,2*rows-1]

            // - Second upper diagonal:
            // result[i,2*rows+i] = block[i+2*rows,i+2*rows]*blockT[i,i+2*rows] for i in [0,rows-1]

            // - First lower diagonal:
            // result[i+rows,i] = block[i,i]*blockT[i+rows,i] for i in [0,2*rows-1]

            // - Second lower diagonal:
            // result[2*rows+i,i] = block[i,i+2*rows]*blockT[i+2*rows,i+2*rows] for i in [0,rows-1]

                // The result matrix here wils be of size length x length, and can be directly added to the
                // upper incidence contribution to get the full Laplacian matrix.


        //     d.) Store as RecursiveToeplitz matrix

    }
};


//=============================== Main function ========================================
int main() {
    int rows = 3;
    int cols = 5;
    
    //Laplacian2D_FullMatrix laplacian(rows, cols);
    Laplacian2D_ToeplitzMatrix laplacian_toeplitz(rows, cols);

    return 0;
}
