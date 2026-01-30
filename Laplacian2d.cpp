//NOTE: THE FOLLOWING CODE IS A WORK IN PROGRESS -- IT MAY NOT RUN CORRECTLY AS IS.
//=========== Classes for Block Toeplitz matrix representation and operations ============
#include "Laplacian2d.h"

Laplacian2D_ToeplitzMatrix::Laplacian2D_ToeplitzMatrix(const int rows, const int cols) {
    // Generate the Laplacian matrix with non-homogeneous spacing and variable k
    generateMatrix(rows, cols);
}
    
Laplacian2D_ToeplitzMatrix::~Laplacian2D_ToeplitzMatrix()
{

}

double Laplacian2D_ToeplitzMatrix::k_func(double x, double y) {
    // Example variable k function; modify as needed
    return 1.0;
}

void Laplacian2D_ToeplitzMatrix::generateMatrix(const int rows, const int cols) 
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
    IncidenceMatrix.Diags[1] = -1; // diagonal for lower incidence matrix
    IncidenceMatrix.Vals[0] = rowBlockPtr;
    IncidenceMatrix.Vals[1] = lowerIncidencePtr;
    //Incidence=&IncidenceMatrix;
    
    //     b.) Compute negative transpose
    Incidence=IncidenceMatrix.Clone(1.0); //The Clone seems necessary to make the Laplacian usable multiple times
    Matrix* negTransPtr = Incidence->negativeTranspose();
    //BlockToeplitz negTrans=*negTrans;
    Incidence_T=negTransPtr; //Maybe a  memory leak, but it works for now
    //Incidence = Incidence_T->negativeTranspose(); 
    //BlockToeplitz* negTransBlockPtr = static_cast<BlockToeplitz*>(negTransPtr); //cast to BlockToeplitz pointer for further operations

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
        W_ee_matrix.Diag_Vals[i] = k_val/(dx*dx)+k_val/(dy*dy);  // use dx (or dy) for spacing
    }

    Diagonal=W_ee_matrix.Clone(1.0);
    //Incidence->printFullMatrix();
    Tmp1= Vectord(Incidence->rows());
    Tmp2=Vectord(Incidence->rows());
    //Incidence_T->printFullMatrix();
}

void Laplacian2D_ToeplitzMatrix::Laplacian2d(const Vectord& input,Vectord& result)
{
    //std::cout << "a\n";
    //Incidence->printFullMatrix();
    //Incidence_T->printFullMatrix();
    //input.print();
    Incidence->printFullMatrix();
    Incidence_T->printFullMatrix();
    
    Incidence->MatMul(input,Tmp1);
    //b2.print();
    //std::cout << "b\n";
    //Vectord Wb(Diagonal->rows());
    Diagonal->MatMul(Tmp1,Tmp2);
    //Wb.print();
    //std::cout << "c\n";
    //Vectord final_b(Incidence_T->rows());
    Incidence_T->MatMul(Tmp2,result);

}

//Function that makes the COO Laplacian matrix
COO Laplacian2D_ToeplitzMatrix::COO_Laplacian()
{
    //COO Incidence_COO(*static_cast<BlockToeplitz*>(Incidence));
    COO Incidence_T_COO(*static_cast<BlockToeplitz*>(Incidence_T));
    DiagonalMatrix* Diag(static_cast<DiagonalMatrix*>(Diagonal));
    for(int i = 0; i<Incidence_T_COO.Num_Vals;i++)
    {
        std::get<2>(Incidence_T_COO.Array[i])*=sqrt(Diag->Diag_Vals[std::get<0>(Incidence_T_COO.Array[i])]);
    }

    int N=5*Incidence_T_COO.rows();

    COO result(Incidence_T_COO.rows(), Incidence_T_COO.rows(), N);

    int n = 0;
    double sum;
    bool changed;
    for(int r = 0; r<Incidence_T_COO.rows(); r++)
    {
        //printf("r=%d\n",r);
        for(int c = 0; c<Incidence_T_COO.rows(); c++)
        {
            //printf("c=%d\n",c);
            //Create values on and above main diagonal
            if(r<=c) 
            {
                changed=false;
                sum = 0.0;
                for(int i = 0; i<Incidence_T_COO.Num_Vals; i++) 
                //i was suupposed to start from k to avoid unnecessary steps, but something went wrong
                {
                    if(r == std::get<0>(Incidence_T_COO.Array[i]))
                    {
                        for(int j = i; j<Incidence_T_COO.Num_Vals; j++)
                        {
                            if(c == std::get<0>(Incidence_T_COO.Array[j]) 
                                && std::get<1>(Incidence_T_COO.Array[i])==std::get<1>(Incidence_T_COO.Array[j]))
                            {
                                sum += std::get<2>(Incidence_T_COO.Array[j])*std::get<2>(Incidence_T_COO.Array[i]);
                                changed= true;
                                break;
                            }
                            if(c<std::get<0>(Incidence_T_COO.Array[j]))
                            {
                                break;
                            }
                        }
                    }
                    if(r<std::get<0>(Incidence_T_COO.Array[i]))
                    {
                        if(changed==true)
                        {
                            std::get<2>(result.Array[n]) = -sum;
                            std::get<0>(result.Array[n]) = r;
                            std::get<1>(result.Array[n]) = c;
                            n++;
                        }
                        break;
                    }
                }
            }
            //Check for pre-existing values due to symmetry, insert values below main diagonal
            else
            {
                for(int i = 0; i<n; i++) 
                {   
                    if(r == std::get<1>(result.Array[i]) && c == std::get<0>(result.Array[i]))
                    {
                        std::get<1>(result.Array[n]) = std::get<0>(result.Array[i]);
                        std::get<0>(result.Array[n]) = std::get<1>(result.Array[i]);
                        std::get<2>(result.Array[n]) = std::get<2>(result.Array[i]);
                        n++;
                    }
                }
            }
        }
    }
    //Set last value
    if(changed==true)
        {
            std::get<2>(result.Array[n]) = -sum;
            std::get<0>(result.Array[n]) = Incidence_T_COO.rows()-1;
            std::get<1>(result.Array[n]) = Incidence_T_COO.rows()-1;
            n++;
        }
    return result;
}

CSR Laplacian2D_ToeplitzMatrix::CSR_Laplacian()
{
    //CSR Incidence_CSR(*static_cast<BlockToeplitz*>(Incidence));
    CSR Incidence_T_CSR(*static_cast<BlockToeplitz*>(Incidence_T));
    DiagonalMatrix* Diag(static_cast<DiagonalMatrix*>(Diagonal));
    for(int i = 0; i<Incidence_T_CSR.rows();i++)
    {
        for(int j = Incidence_T_CSR.Rows[i]; j<Incidence_T_CSR.Rows[i+1];j++)   
        { 
            Incidence_T_CSR.Vals[j]*=sqrt(Diag->Diag_Vals[i]);
        }
    }

    int N = 5*Incidence_T_CSR.rows();

    CSR result(Incidence_T_CSR.rows(), Incidence_T_CSR.rows(), N);

    double sum;
    bool changed;
    for(int r = 0; r<Incidence_T_CSR.rows(); r++)
    {
        //printf("r = %d\n", r);
        result.Rows[r+1] = result.Rows[r];
        for(int c = 0; c<Incidence_T_CSR.rows(); c++)
        {
            //Create values on and above main diagonal
            //printf("c = %d\n", c);
            if(r<=c) 
            {
                sum = 0.0;
                changed = false;
                for(int i = Incidence_T_CSR.Rows[r]; i<Incidence_T_CSR.Rows[r+1]; i++)
                {
                    for(int j = Incidence_T_CSR.Rows[c]; j<Incidence_T_CSR.Rows[c+1]; j++)
                    {
                        if(Incidence_T_CSR.Cols[i]==Incidence_T_CSR.Cols[j])
                        {
                            sum += Incidence_T_CSR.Vals[j]*Incidence_T_CSR.Vals[i];
                            changed = true;
                            break;
                        }
                    }
                }
                if(changed==true)
                {
                    result.Vals[result.Rows[r+1]] = -sum;
                    result.Cols[result.Rows[r+1]] = c;
                    result.Rows[r+1]++;
                    //break;
                }
            }
            //Check for pre-existing values due to symmetry, insert values below main diagonal
            else
            {
                for(int i = 0; i<r; i++) 
                {   
                    for(int j = result.Rows[c]; j<result.Rows[c+1];j++)
                    {
                        if(r == result.Cols[j] && c == i)
                        {
                            result.Vals[result.Rows[r+1]] = result.Vals[j];
                            result.Cols[result.Rows[r+1]] = i;
                            result.Rows[r+1]++;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    return result;
}