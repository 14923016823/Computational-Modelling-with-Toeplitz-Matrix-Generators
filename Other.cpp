#include "Other.h"

int main()
{
    const int length = 3;
    double vals[length] = {1,4,6};
    int num_cols = 4;
    int num_rows = 5;
    int cols[length] = {0,3,2};
    int rows[num_rows+1] = {0,2,2,2,2,3};

    tuple t1 = std::make_tuple(0,0,1);
    tuple t2 = std::make_tuple(0,3,4);
    tuple t3 = std::make_tuple(4,2,6);
    COO B(num_rows,num_cols,{t1,t2,t3});
    //B.print();

    CSR A(vals, cols, rows, length, num_rows, num_cols);
    Vectord x(num_cols);
    A.print();

    B.printFullMatrix();
    Matrix* Bt = B.negativeTranspose();
    Bt->print();
    Bt->printFullMatrix();
    
    for(int i=0;i<num_cols;i++)
    {
        x.Vec[i]=i;
    }
    x.print();
        
    A.print();

    Vectord y = A*x;
    y.print();

    Vectord z = B*x;
    z.print();

    const int STlength = 3;
    int STdiags[STlength] = {-2,-1,3};
    double STvals[STlength] = {5,1,4};
    int STwidth = 4;
    int STheight = 4;
    SparseToeplitz C(STheight, STwidth, STlength, STdiags, STvals);
    C.print();

    C.printFullMatrix();

    Vectord a = C*x;
    a.print();
    
    int Onediag[1] = {0};
    double Oneval[1] = {1};
    SparseToeplitz One(1,1,1,Onediag,Oneval);
    One.print();

    CSR C_CSR(C);
    C_CSR.print();

    COO C_COO(C);
    C_COO.print();

    COO One_COO(One);
    One_COO.print();

    CSR One_CSR(One);
    One_CSR.print();

    Matrix* Two = C.Kronecker(One);

    Matrix* Two_COO = One_COO.Kronecker(C_COO);
    *Two_COO *= 2;

    Matrix* Two_CSR = One_CSR.Kronecker(C_CSR);
    *Two_CSR *= 2;

    std::cout << "?" << "\n";

    std::cout << "?" << "\n";

    Vectord in = Vectord({0,1,2,3});

    std::cout << "?" << "\n";

    Vectord out = *Two*in;

    std::cout << "?" << "\n";

    Vectord out_COO = *Two_COO*in;

    std::cout << "?" << "\n";

    Matrix* Two_BT = One.Kronecker(C);
    *Two_BT *= 3;

    std::cout << "?" << "\n";

    Vectord out_BT = *Two_BT*in;

    Two_COO->printFullMatrix();
    Two_BT->printFullMatrix();
    Two_CSR->printFullMatrix();

    Matrix* Three_CSR = Two_CSR->negativeTranspose();
    Matrix* Three_COO = Two_COO->negativeTranspose();

    Three_CSR->printFullMatrix();
    Three_COO->printFullMatrix();
    
    out.print();
    out_COO.print();
    out_BT.print();
   
    return 0;
}