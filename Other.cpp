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

    Vectord z = A*x;
    z.print();

    const int STlength = 1;
    int STdiags[STlength] = {0};
    double STvals[STlength] = {3};
    int STwidth = 2;
    int STheight = 2;
    SparseToeplitz C(STheight, STwidth, STlength, STdiags, STvals);
    int Onediag[1] = {0};
    double Oneval[1] = {1};
    SparseToeplitz One(1,1,1,Onediag,Oneval);
    C.print();
    One.print();

    CSR C_CSR(C);
    C_CSR.print();

    CSR One_CSR(One);
    One_CSR.print();

    Matrix* Two = One.Kronecker(C);

    Matrix* Two_CSR = One_CSR.Kronecker(C_CSR);

    Vectord in = Vectord({2,5});

    Vectord out = *Two*in;

    Vectord out_CSR = *Two_CSR*in;

    out.print();
    out_CSR.print();

    COO C_COO(C);
    C_COO.print();

int i;
//omp_set_num_threads(2);
#pragma omp parallel private(i) num_threads(7)             
{
    
    i = omp_get_thread_num();
    
    printf("Hello World... from thread = %d\n", i);
} 
    
    return 0;
}