#include "Other.h"

int main()
{
    int length = 3;
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

    std::srand(10);
    std::experimental::reseed(10);
    int r_rows = 4000;
    int r_cols = 4000;
    int r_prob = 5;

    SparseToeplitz D = r_ST(r_rows,r_cols,r_prob);
    //D.printFullMatrix();
    D.print();
    Vectord x2 = r_Vec(r_cols);
    //x2.print();

    auto start = std::chrono::steady_clock::now();
    Vectord r = D*x2;
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "t = " << diff.count() << '\n';
    r.Vec[0] = r.Vec[0];
    //r.print();

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
    start = std::chrono::steady_clock::now();
    Vectord y = A*x;
    end = std::chrono::steady_clock::now();
    diff = end - start;
    std::cout << "t = " << diff.count() << '\n';
    y.print();

    start = std::chrono::steady_clock::now();
    Vectord z = B*x;
    end = std::chrono::steady_clock::now();
    diff = end - start;
    std::cout << "t = " << diff.count() << '\n';
    z.print();

    int STlength = 3;
    int STdiags[STlength] = {-2,-1,3};
    double STvals[STlength] = {5,1,4};
    int STwidth = 4;
    int STheight = 4;
    SparseToeplitz C(STheight, STwidth, STlength, STdiags, STvals);
    C.print();

    std::cout << "C: \n";
    C.printFullMatrix();

    start = std::chrono::steady_clock::now();
    Vectord a = C*x;
    end = std::chrono::steady_clock::now();
    diff = end - start;
    std::cout << "t = " << diff.count() << '\n';
    a.print();
    
    int Onediag[1] = {0};
    double Oneval[1] = {1};
    SparseToeplitz One(3,2,1,Onediag,Oneval);
    One.print();
    std::cout << "One: \n";
    One.printFullMatrix();

    CSR C_CSR(C);
    C_CSR.print();

    COO C_COO(C);
    C_COO.print();

    COO One_COO(One);
    One_COO.print();

    CSR One_CSR(One);
    One_CSR.print();

    Matrix* Two = C.Kronecker(One);
    std::cout << "Two: \n";
    Two->printFullMatrix();

    Vectord in(Two->cols());
    for(int i=0;i<Two->cols();i++)
    {
        in.Vec[i]=i;
    }

    Matrix* Two_COO = C_COO.Kronecker(One_COO);
    *Two_COO *= 2;

    Matrix* Two_CSR = C_CSR.Kronecker(One_CSR);
    *Two_CSR *= 2;

    //Vectord in = Vectord({0,1,2,3});

    start = std::chrono::steady_clock::now();
    Vectord out = *Two*in;
    end = std::chrono::steady_clock::now();
    diff = end - start;
    std::cout << "t = " << diff.count() << '\n';

    start = std::chrono::steady_clock::now();
    Vectord out_COO = *Two_COO*in;
    end = std::chrono::steady_clock::now();
    diff = end - start;
    std::cout << "t = " << diff.count() << '\n';

    Matrix* Two_BT = C.Kronecker(One);
    *Two_BT *= 3;

    start = std::chrono::steady_clock::now();
    Vectord out_BT = *Two_BT*in;
    end = std::chrono::steady_clock::now();
    diff = end - start;
    std::cout << "t = " << diff.count() << '\n';

    Two_COO->printFullMatrix();
    Two_BT->printFullMatrix();
    Two_CSR->printFullMatrix();


    start = std::chrono::steady_clock::now();
    Matrix* Three_CSR = Two_CSR->negativeTranspose();
    end = std::chrono::steady_clock::now();
    diff = end - start;
    std::cout << "t = " << diff.count() << '\n';
    
    Matrix* Three_COO = Two_COO->negativeTranspose();

    Three_CSR->printFullMatrix();

    try 
    {
    Three_COO->printFullMatrix();
    }
    catch (const char* msg)
    {
        std::cout << msg <<std::endl;
    }
    
    
    out.print();
    out_COO.print();
    out_BT.print();
   
    return 0;
}