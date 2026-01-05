#include "Other.h"

int main()
{
    const int length = 3;
    double vals[length] = {1,4,6};
    int num_cols = 4;
    const int num_rows = 5;
    int cols[length] = {0,3,2};
    int rows[num_rows+1] = {0,2,2,2,2,3};

    tuple t1 = std::make_tuple(1,0,0);
    tuple t2 = std::make_tuple(4,0,3);
    tuple t3 = std::make_tuple(6,4,2);
    COO B({t1,t2,t3},num_rows,num_cols);
    B.print();

    CSR A(vals, cols, rows, length, num_rows, num_cols);
    Vectord x(num_cols);
    
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

    Vectord a = C*x;
    a.print();
    
    int Onediag[1] = {0};
    double Oneval[1] = {1};
    SparseToeplitz One(1,1,1,Onediag,Oneval);
    One.print();

    CSR C_CSR(C);
    C_CSR.print();

    CSR One_CSR(One);
    One_CSR.print();

    Matrix* Two = C.Kronecker(One);
    std::cout << (*Two).cols() << std::endl;

    Matrix* Two_CSR = One_CSR.Kronecker(C_CSR);
    *Two_CSR *= 2;

    Vectord in = Vectord({0,1,2,3});

    Vectord out = *Two*in;

    Vectord out_CSR = *Two_CSR*in;

    Matrix* Two_BT = One.Kronecker(C);
    *Two_BT *= 3;

    Vectord out_BT = *Two_BT*in;

    out.print();
    out_CSR.print();
    out_BT.print();

    COO C_COO(C);
    C_COO.print();

int i;
   
    return 0;
}