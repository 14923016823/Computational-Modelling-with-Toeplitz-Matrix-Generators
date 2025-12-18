/* Laplacian matrix generation based on a looping stencil - 2D case */
#include <iostream>
#include <vector>
#include <iomanip>
#include <stdexcept>

//======================================== TO DOs ========================================
// - PRIORITY: Check matrix storage format (row-major vs column-major) => Toeplitz
// - Implement "variable" k
// - Implment non-homogeneous mesh spacing
// - Scale up to 3D case?


//================ Classes for full matrix representation and operations =================
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


//=========== Classes for Block Toeplitz matrix representation and operations ============
// THE FOLLOWING Vectord AND RecursiveToeplitz CLASSES HAVE BEEN COPIED FROM BlockToeplitz.cpp
class Vectord {
public:
    std::vector<double> Vec;

    Vectord(int length = 0) : Vec(length, 0.0) {}

    int Length() const { return static_cast<int>(Vec.size()); }

    void PrintVector() const {
        for (int i = 0; i < Length(); ++i) {
            printf("vec[%d]=%f\n", i, Vec[i]);
        }
        printf("done\n");
    }

    void Sum(const Vectord& VecIn) {
        if (VecIn.Length() != Length()) {
            throw std::invalid_argument("You can't sum vectors with different sizes");
        }
        for (int i = 0; i < Length(); ++i) {
            Vec[i] += VecIn.Vec[i];
        }
    }
};

class RecursiveToeplitz {
public:
    int Num_Rows = 0;
    int Num_Cols = 0;
    bool IsRecursive = true;

    int Num_Diags = 0;
    int Num_Rows_SubMatrixes = 0;
    int Num_Cols_SubMatrixes = 0;

    std::vector<int> Diags;                    //
    std::vector<RecursiveToeplitz> Vals;       // sub-Toeplitz matrices (only if recursive)
    std::vector<double> ValsD;                 // diagonal values (only if leaf)

    // default constructor (needed for vector::resize)
    RecursiveToeplitz() : IsRecursive(false), Num_Rows(0), Num_Cols(0),
                          Num_Diags(0), Num_Rows_SubMatrixes(0), Num_Cols_SubMatrixes(0) {}

    // recursive constructor
    RecursiveToeplitz(int nrows, int ncols, int ndiags, int nrsub, int ncsub) {
        IsRecursive = true;
        Num_Rows_SubMatrixes = nrsub;
        Num_Cols_SubMatrixes = ncsub;
        Num_Rows = nrows;
        Num_Cols = ncols;
        Num_Diags = ndiags;
        Diags.resize(ndiags);
        Vals.resize(ndiags); // requires configuration later on
    }

    // leaf constructor (diagonals are scalar values)
    RecursiveToeplitz(int nrows, int ncols, int ndiags) {
        IsRecursive = false;
        Num_Rows = nrows;
        Num_Rows_SubMatrixes = 1;
        Num_Cols = ncols;
        Num_Cols_SubMatrixes = 1;
        Num_Diags = ndiags;
        Diags.resize(ndiags);
        ValsD.resize(ndiags, 0.0);
    }

    // Corrected copy-and-scale constructor
    RecursiveToeplitz(const RecursiveToeplitz& other, double c){
        Num_Rows = other.Num_Rows;
        Num_Cols = other.Num_Cols;
        IsRecursive = other.IsRecursive;
        Num_Diags = other.Num_Diags;
        Num_Rows_SubMatrixes = other.Num_Rows_SubMatrixes;
        Num_Cols_SubMatrixes = other.Num_Cols_SubMatrixes;

        Diags = other.Diags;  // always the same

        if (!IsRecursive) {
            // leaf: scale every diagonal value
            ValsD.resize(Num_Diags);
            for (int i = 0; i < Num_Diags; i++)
                ValsD[i] = other.ValsD[i] * c;
        }
        else {
            // recursive: COPY the submatrices, then scale them
            Vals = other.Vals;              // <--- FIX 1 (deep copy via vector copy)
            for (int i = 0; i < Num_Diags; i++)
                Vals[i].ScalarMultiplication(c);   // <--- FIX 2 (apply scaling here)
        }
    }

    Vectord VectorMultiplication(const Vectord vec) const {
        if (!IsRecursive) {
            return VectorMultiplicationD(vec);
        }

        if (vec.Length() != Num_Cols) {
            throw std::invalid_argument("Vector length and matrix columns don't match (recursive).");
        }

        Vectord result(Num_Rows);
        int Num_blockrows = Num_Rows / Num_Rows_SubMatrixes;
        int Num_blockcols = Num_Cols / Num_Cols_SubMatrixes;

        for (int blockrow = 0; blockrow < Num_blockrows; blockrow++) {
            //int row = blockrow *;
            Vectord subres(Num_Rows_SubMatrixes); // initialized to zeros

            // accumulate contributions from each diagonal/block
            for (int j = 0; j < Num_Diags; j++) {
                int blockcol = Diags[j] + blockrow;

                // ensure the whole sub-block fits in input vector
                if (blockcol >= 0 && blockcol<Num_blockcols) {
                    Vectord subinput(Num_Cols_SubMatrixes);
                    int col_start=blockcol*Num_Cols_SubMatrixes;

                    for (int k = 0; k < Num_Cols_SubMatrixes; ++k) {
                        subinput.Vec[k] = vec.Vec[col_start + k]; 
                    }
                    Vectord subsubres = Vals[j].VectorMultiplication(subinput);
                    subres.Sum(subsubres);
                }
            }
            // copy accumulated block into result
            int row_start=blockrow*Num_Rows_SubMatrixes;
            for (int k = 0; k < Num_Rows_SubMatrixes; ++k) {
                result.Vec[row_start + k] = subres.Vec[k];
            }
        }
        return result;
    }

    Vectord VectorMultiplicationD(const Vectord& vec) const {
        if (vec.Length() != Num_Cols) {
            throw std::invalid_argument("Vector and Matrix size dont match (leaf).");
        }
        Vectord result(Num_Rows);
        for (int i = 0; i < Num_Rows; i++) {
            double acc = 0.0;
            for (int j = 0; j < Num_Diags; j++) {
                int idx = Diags[j] + i;
                if (idx >= 0 && idx < Num_Cols) {
                    acc += vec.Vec[idx] * ValsD[j];
                }
            }
            result.Vec[i] = acc;
        }
        return result;
    }

    void ScalarMultiplication(double c) {
        if(IsRecursive==false) {
            for(int i=0;i<Num_Diags;i++) {
                ValsD[i]*=c;
            }
        }
        else {
            for(int i=0;i<Num_Diags;i++) {
                Vals[i].ScalarMultiplication(c);
            }
        }
    }

    void Kronecker(RecursiveToeplitz B) {
        //if you add a new matrix at the bottom of the chain every Num_Rows needs to be changed
        Num_Rows*=B.Num_Rows;
        Num_Cols*=B.Num_Cols;
        Num_Cols_SubMatrixes*=B.Num_Cols;
        Num_Rows_SubMatrixes*=B.Num_Rows;

        if(IsRecursive==true) {
            for(int i=0;i<Num_Diags;i++) {
                Vals[i].Kronecker(B);
            }
        }
        else {   
            IsRecursive=true; //after kronecker this will become a Blockmatrix
            Vals.resize(Num_Diags); //we will now use the Valsmatrix instead of the ValsD
            
            for(int i=0;i<Num_Diags;i++) {
                Vals[i]=RecursiveToeplitz(B,ValsD[i]);
            }
            //ValsD.clear();
        }
    }

    //THE FOLLOWING FUNCTIONS AND IMPLEMENTATIONS HAVE BEEN ADDED TO THE CLASS POST-COPY - NOT INCLUDED IN
    //THE ORIGINAL BlockToeplitz.cpp file
    void negativeTranspose() {
        //negate and transpose the matrix
        for(int i=0;i<Num_Diags;i++) {
            Diags[i]=-Diags[i];
            if(IsRecursive==false) {
                ValsD[i]=-ValsD[i];
            }
            else {
                Vals[i].negativeTranspose();
            }
        }
    }
    
    void PrintRecursiveToeplitz(int level = 0) const {
        //print function for debugging
        std::string indent(level * 2, ' ');
        if (IsRecursive) {
            std::cout << indent << "Recursive Toeplitz Matrix: " << Num_Rows << "x" << Num_Cols << "\n";
            for (int i = 0; i < Num_Diags; ++i) {
                std::cout << indent << "  Diagonal " << Diags[i] << ":\n";
                Vals[i].PrintRecursiveToeplitz(level + 2);
            }
        } else {
            std::cout << indent << "Leaf Toeplitz Matrix: " << Num_Rows << "x" << Num_Cols << "\n";
            for (int i = 0; i < Num_Diags; ++i) {
                std::cout << indent << "  Diagonal " << Diags[i] << " Value: " << ValsD[i] << "\n";
            }
        }
    }

    void PrintFullToeplitz() const {
        //print function for full dense expansion
        std::vector<std::vector<double>> M = ExpandToDense();
        std::cout << "\nFull Dense Expansion (" << Num_Rows << "x" << Num_Cols << ")\n";

        for (int i = 0; i < Num_Rows; ++i) {
            for (int j = 0; j < Num_Cols; ++j)
                std::cout << std::setw(4) << std::left << M[i][j];
            std::cout << "\n";
        }
    }

    // Support function for full matrix print function
    std::vector<std::vector<double>> ExpandToDense() const {

        // Print statement for a Leaf Toeplitz matrix
        if (!IsRecursive) {
            std::vector<std::vector<double>> leafM(Num_Rows, std::vector<double>(Num_Cols, 0.0));
            for (int d = 0; d < Num_Diags; d++) {
                for (int i = 0; i < Num_Rows; i++) {
                    int j = i + Diags[d];
                    if (j >= 0 && j < Num_Cols)
                        leafM[i][j] = ValsD[d]; // stores diagonal value
                    
                    // if symmetric Toeplitz, uncomment this:
                    // if (Diags[d] > 0 && i - Diags[d] >= 0)
                    //    M[i][i - Diags[d]] = ValsD[d];
                }
            }
            return leafM;
        }

        // Print statement for Recursive Toeplitz matrix
        int block_rows = Vals[0].Num_Rows;
        int block_cols = Vals[0].Num_Cols;
        std::vector<std::vector<double>> M(Num_Rows * block_rows, std::vector<double>(Num_Cols * block_cols, 0.0));

        for (int k = 0; k < Num_Diags; k++) {
            int d = Diags[k];  // block diagonal offset
            for (int i = 0; i < Num_Rows - d; i++) {
                int j = i + d;
                auto block = Vals[k].ExpandToDense();  // recursively expand sub-block
                for (int bi = 0; bi < block_rows; bi++) {
                    for (int bj = 0; bj < block_cols; bj++) {
                        M[i * block_rows + bi][j * block_cols + bj] = block[bi][bj];
                        if (d != 0) // mirror for symmetric blocks
                            M[j * block_rows + bi][i * block_cols + bj] = block[bi][bj];
                    }
                }
            }
        }
        return M;
    }
};

class Laplacian2D_ToeplitzMatrix {
    public:
    Laplacian2D_ToeplitzMatrix(const int rows, const int cols) {
        // Generate the Laplacian matrix
        generateMatrix(rows, cols);
        //RecursiveToeplitz LaplacianMatrix = RecursiveToeplitz(2, 1, 1, 1, 2); //initialize empty matrix
    }
    
    //Function to pring the laplacian matrix
    void PrintLaplacianMatrix(const RecursiveToeplitz& laplacian) {
    
    }
    
    private:
    void generateMatrix(const int rows, const int cols) {
        // Implementation for generating the Laplacian matrix

        // 1.)Upper incidence matrix generation
        //     a.) Row block matrix generation
        RecursiveToeplitz rowBlockMatrix = RecursiveToeplitz(cols, cols, 3); //leaf matrix with 3 diagonals
        rowBlockMatrix.Diags[0] = -(cols - 1);
        rowBlockMatrix.Diags[1] = 0;
        rowBlockMatrix.Diags[2] = 1;
        rowBlockMatrix.ValsD[0] = 1;
        rowBlockMatrix.ValsD[1] = -1;
        rowBlockMatrix.ValsD[2] = 1;

        //     b.) Identity matrix generation for lower incidence matrix
        RecursiveToeplitz identityMatrixCols(cols, cols, 1); //leaf identity matrix
        identityMatrixCols.Diags[0] = 0;
        identityMatrixCols.ValsD[0] = 1;

        // 2.) Lower incidence matrix generation
        //     a.) Column block matrix generation
        RecursiveToeplitz columnBlockMatrix(rows, rows, 3); //leaf matrix with 3 diagonals

        //     b.) Column block matrix kronecker product with identity matrix (cols)
        RecursiveToeplitz lowerIncidenceMatrix = columnBlockMatrix;
        lowerIncidenceMatrix.Kronecker(identityMatrixCols);

        // 3.) Final Laplacian matrix assembly
        //     a.) Concatenate upper and lower incidence matrices
        
        
        //     b.) Compute negative transposes
        RecursiveToeplitz negTransposeUpperBlock = rowBlockMatrix;
        RecursiveToeplitz negTransposeLower = lowerIncidenceMatrix;
        negTransposeUpperBlock.negativeTranspose();
        lowerIncidenceMatrix.negativeTranspose();

        //negTransposeUpperBlock.PrintFullToeplitz();

        //     c.) Multiply incidence matrix with its negative transpose to get Laplacian
            RecursiveToeplitz upperIncidenceBlock = RecursiveToeplitz(rows, rows, 3);

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
