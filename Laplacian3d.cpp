#include "Laplacian3d.h"

Laplacian3D_ToeplitzMatrix::Laplacian3D_ToeplitzMatrix(const int rows, const int cols, const int arrays) {
    // Generate the Laplacian matrix with non-homogeneous spacing and variable k
    generateMatrix(rows, cols, arrays);
}
    
Laplacian3D_ToeplitzMatrix::~Laplacian3D_ToeplitzMatrix()
{

}


double Laplacian3D_ToeplitzMatrix::k_func(double x, double y, double z) {
    // Example variable k function; modify as needed
    return 1.0;
}

void Laplacian3D_ToeplitzMatrix::generateMatrix(const int rows, const int cols, const int arrays) 
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
    SparseToeplitz arrayBlockMatrix(arrays, arrays, 3); //leaf array block matrix with 3 diagonals
    arrayBlockMatrix.Diags[0] = -(arrays - 1);
    arrayBlockMatrix.Diags[1] = 0;
    arrayBlockMatrix.Diags[2] = 1;
    arrayBlockMatrix.Vals[0] = 1;
    arrayBlockMatrix.Vals[1] = -1;
    arrayBlockMatrix.Vals[2] = 1;

    //     b.) Identity matrix generation for lower incidence matrix
    SparseToeplitz identityMatrixArrays(arrays, arrays, 1); //leaf identity matrix
    identityMatrixArrays.Diags[0] = 0;
    identityMatrixArrays.Vals[0] = 1;
    SparseToeplitz identityMatrix_2D(rows*cols, rows*cols, 1);
    identityMatrix_2D.Diags[0] = 0;
    identityMatrix_2D.Vals[0] = 1;

    //     c.) Column block matrix kronecker product with identity matrix (arrays)
    lowerIncidenceMatrix_3D = (arrayBlockMatrix.Kronecker(identityMatrix_2D))->Clone(1.0); // lower incidence matrix
    BlockToeplitz* lowerIncidencePtr_3D = static_cast<BlockToeplitz*>(lowerIncidenceMatrix_3D);  //cast to BlockToeplitz pointer for further operations
    lowerIncidenceMatrix_3D = lowerIncidencePtr_3D->Clone(1.0);

    // 4.) Final Laplacian matrix assembly
    //     a.) Concatenate upper and lower incidence matrices
    int ndiags = 2; //number of diagonals in incidence matrix
    BlockToeplitz IncidenceMatrix_2D(rows * cols * 2, rows * cols, ndiags); // incidence matrix with 2 block rows
    IncidenceMatrix_2D.Diags[0] = 0; // diagonal for upper incidence matrix
    IncidenceMatrix_2D.Diags[1] = -1; // diagonal for lower incidence matrix
    IncidenceMatrix_2D.Vals[0] = rowBlockPtr;
    IncidenceMatrix_2D.Vals[1] = lowerIncidencePtr;
    //Incidence=&IncidenceMatrix;

    //IncidenceMatrix_2D.printFullMatrix();
    upperIncidenceMatrix_3D = (identityMatrixArrays.Kronecker(IncidenceMatrix_2D))->Clone(1.0); // upper incidence matrix
    BlockToeplitz* BlockPtr_2D = static_cast<BlockToeplitz*>(upperIncidenceMatrix_3D);  //cast to BlockToeplitz pointer for further operations
    //BlockPtr_2D->printFullMatrix();
    upperIncidenceMatrix_3D = BlockPtr_2D->Clone(1.0);
    
    BlockToeplitz IncidenceMatrix_3D(rows * cols * arrays * 3, rows * cols * arrays, ndiags); // incidence matrix with 2 block rows
    IncidenceMatrix_3D.Diags[0] = 0; // diagonal for upper incidence matrix
    IncidenceMatrix_3D.Diags[1] = -1; // diagonal for lower incidence matrix
    IncidenceMatrix_3D.Vals[0] = BlockPtr_2D;
    IncidenceMatrix_3D.Vals[1] = lowerIncidencePtr_3D;

    leftIncidenceT_Matrix_3D = (upperIncidenceMatrix_3D->negativeTranspose())->Clone(1.0);
    BlockToeplitz* left_T_3D = static_cast<BlockToeplitz*>(leftIncidenceT_Matrix_3D);
    leftIncidenceT_Matrix_3D = left_T_3D->Clone(1.0);

    rightIncidenceT_Matrix_3D = (lowerIncidenceMatrix_3D->negativeTranspose())->Clone(1.0);
    BlockToeplitz* right_T_3D = static_cast<BlockToeplitz*>(rightIncidenceT_Matrix_3D);
    rightIncidenceT_Matrix_3D = right_T_3D->Clone(1.0);
    
    //     b.) Compute negative transpose
    Incidence=IncidenceMatrix_3D.Clone(1.0);
    //Incidence->printFullMatrix(); //This print also returns an error because of inconsistent block sizes
    Matrix* negTransPtr = IncidenceMatrix_3D.negativeTranspose();
    //BlockToeplitz negTrans=*negTrans;
    Incidence_T=negTransPtr;
    //BlockToeplitz* negTransBlockPtr = static_cast<BlockToeplitz*>(negTransPtr); //cast to BlockToeplitz pointer for further operations

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
    W_ee_matrix.Diag_Vals[i] = k_val/(pow(dx,2))+k_val/(pow(dy,2))+k_val/(pow(dz,2));  // use dx (or dy) for spacing
}

    Diagonal=W_ee_matrix.Clone(1.0);
    Tmp1= Vectord(Incidence->rows());
    Tmp2=Vectord(Incidence->rows());
}

void Laplacian3D_ToeplitzMatrix::Laplacian3d(const Vectord& input,Vectord& result)
{
    upperIncidenceMatrix_3D->MatMulAdd(input,0,Tmp1,0);
    lowerIncidenceMatrix_3D->MatMulAdd(input,0,Tmp1,upperIncidenceMatrix_3D->rows());
    Diagonal->MatMul(Tmp1,Tmp2);
    leftIncidenceT_Matrix_3D->MatMulAdd(Tmp2,0,result,0);
    rightIncidenceT_Matrix_3D->MatMulAdd(Tmp2,leftIncidenceT_Matrix_3D->cols(),result,0);
}

COO Laplacian3D_ToeplitzMatrix::COO_Laplacian()
{
    //COO Incidence_COO(*static_cast<BlockToeplitz*>(Incidence));
    COO Incidence_COO(*static_cast<BlockToeplitz*>(Incidence));
    COO Incidence_T_COO(*static_cast<COO*>(Incidence_COO.negativeTranspose()));
    Incidence_T_COO.printFullMatrix();
    std::cout << "left:\n";
    leftIncidenceT_Matrix_3D->printFullMatrix();
    std::cout << "right:\n";
    rightIncidenceT_Matrix_3D->printFullMatrix();
    DiagonalMatrix* Diag(static_cast<DiagonalMatrix*>(Diagonal));
    for(int i = 0; i<Incidence_T_COO.Num_Vals;i++)
    {
        std::get<2>(Incidence_T_COO.Array[i])*=sqrt(Diag->Diag_Vals[std::get<0>(Incidence_T_COO.Array[i])]);
    }

    int N = 0;
    int val_tot = 0;
    int n_c;

    for(int c=0;c<Incidence_T_COO.rows();c++)
    {
        n_c = 0;
        for(int i = val_tot;i<Incidence_T_COO.Num_Vals;i++)
        {
            if(std::get<0>(Incidence_T_COO.Array[i])==c)
            {
                n_c++;
                val_tot++;
            }
            if(std::get<0>(Incidence_T_COO.Array[i])>c)
            {
                break;
            }
        }
        N += (2*n_c-3);
    }

    N=7*Incidence_T_COO.rows();

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

CSR Laplacian3D_ToeplitzMatrix::CSR_Laplacian()
{
    //CSR Incidence_CSR(*static_cast<BlockToeplitz*>(Incidence));
    CSR Incidence_CSR(*static_cast<BlockToeplitz*>(Incidence));
    CSR Incidence_T_CSR(*static_cast<CSR*>(Incidence_CSR.negativeTranspose()));
    DiagonalMatrix* Diag(static_cast<DiagonalMatrix*>(Diagonal));
    for(int i = 0; i<Incidence_T_CSR.rows();i++)
    {
        for(int j = Incidence_T_CSR.Rows[i]; j<Incidence_T_CSR.Rows[i+1];j++)   
        { 
            Incidence_T_CSR.Vals[j]*=sqrt(Diag->Diag_Vals[i]);
        }
    }

    int N = 0;

    for(int c=0;c<Incidence_T_CSR.rows();c++)
    {
        N += (2*(Incidence_T_CSR.Rows[c+1]-Incidence_T_CSR.Rows[c])-3);
    }

    N = 7*Incidence_T_CSR.rows();

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
                            if(r==9 && c==9)
                            {
                                printf("%d\n",Incidence_T_CSR.Cols[i]);
                            }
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