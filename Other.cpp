#include "Other.h"

int main()
{/*
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

/*
    //B.print();

    //std::srand(10);
    //std::experimental::reseed(10);
    int r_rows = 4000;
    int r_cols = 4000;
    int r_prob = 5;

    SparseToeplitz D = r_ST(r_rows,r_cols,r_prob);
    CSR D_CSR(D);
    COO D_COO(D);
    //D.printFullMatrix();
    //D.print();
    Vectord x2 = r_Vec(r_cols);
    for(int i=0;i<r_cols;i++)
    {
        x2.Vec[i]=(i*3+8)/0.29836298;
    }
    //x2.print();

    auto start = std::chrono::steady_clock::now();
    Vectord r = D*x2;
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "t = " << diff.count() << '\n';
    r.Vec[0] = r.Vec[0];

    Vectord r_CSR = D_CSR*x2;
    Vectord r_COO = D_COO*x2;

    Vectord r_diff(r_cols);
    for(int i=0;i<r_cols;i++)
    {
        r_diff.Vec[i] = r.Vec[i]-r_CSR.Vec[i];
    }
    r_diff.print();
    //r.print();

    CSR A(vals, cols, rows, length, num_rows, num_cols);
    Vectord x(num_cols);
    A.print();

    A.printFullMatrix();
    Matrix* Bt = A.negativeTranspose();
    Bt->print();
    Bt->printFullMatrix();
    
    for(int i=0;i<num_cols;i++)
    {
        x.Vec[i]=(i*3+8);
    }
    x.print();
        
    A.print();
    //start = std::chrono::steady_clock::now();
    Vectord y = A*x;
    //end = std::chrono::steady_clock::now();
    //diff = end - start;
   // std::cout << "t = " << diff.count() << '\n';
    y.print();

    //start = std::chrono::steady_clock::now();
    Vectord z = B*x;
    //end = std::chrono::steady_clock::now();
    //diff = end - start;
    //std::cout << "t = " << diff.count() << '\n';
    z.print();
*/
    int STlength = 3;
    int STdiags[STlength] = {-1,0,1};
    double STvals[STlength] = {-1,2,-1};
    int STwidth = 4;
    int STheight = 4;
    Vectord x(STwidth);
    for(int i=0;i<STwidth;i++)
    {
        x.Vec[i]=i;
    }
    SparseToeplitz C(STheight, STwidth, STlength, STdiags, STvals);
    //C.print();

    //std::cout << "C: \n";
    //C.printFullMatrix();

    //auto start = std::chrono::steady_clock::now();
    Vectord a = C*x;
    //auto end = std::chrono::steady_clock::now();
    //std::chrono::duration<double> diff = end - start;
    //std::cout << "t = " << diff.count() << '\n';
    //a.print();
    
    int Onediag[1] = {0};
    double Oneval[1] = {1};
    SparseToeplitz One(3,2,1,Onediag,Oneval);
    //One.print();
    //std::cout << "One: \n";
    //One.printFullMatrix();

    //CSR C_CSR(C);
    //C_CSR.print();

    //COO C_COO(C);
    //C_COO.print();

    //COO One_COO(One);
    //One_COO.print();

    //CSR One_CSR(One);
    //One_CSR.print();

    //start = std::chrono::steady_clock::now();
    //Vectord a_csr = C_CSR*x;
    //end = std::chrono::steady_clock::now();
    //diff = end - start;
    //std::cout << "t = " << diff.count() << '\n';

    //start = std::chrono::steady_clock::now();
    //Vectord a_coo = C_COO*x;
    //end = std::chrono::steady_clock::now();
    //diff = end - start;
    //std::cout << "t = " << diff.count() << '\n';

    Matrix* Two = C.Kronecker(One);
    std::cout << "Two: \n";
    Two->printFullMatrix();
    Two->operator*=(5);
    Two->printFullMatrix();
    BlockToeplitz* Dos = static_cast<BlockToeplitz*>(Two);
    CSR Two_COO(*Dos);
    Two_COO.printFullMatrix();
/*
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
   */
  
    int rows = 3;
    int cols = 3;
    int arrays = 1;

    Vectord vec(cols*rows*arrays);
    for (int i = 0; i < rows * cols * arrays; ++i) {
        vec.Vec[i] = i; // Example initialization
    }
    

    //Laplacian2D_FullMatrix laplacian(rows, cols);

    Laplacian2D_ToeplitzMatrix L2d=Laplacian2D_ToeplitzMatrix(rows,cols);
    COO L2d_COO = L2d.COO_Laplacian();
    //L2d_COO.print();
    CSR L2d_CSR = L2d.CSR_Laplacian();
    //L2d_CSR.printFullMatrix();
    try
    {
        Vectord a(rows*cols);
        L2d.Laplacian2d(vec,a);
        a.print();

        L2d_COO.printFullMatrix();
        L2d_CSR.printFullMatrix();

        Vectord b = L2d_CSR*vec;
        b.print();
    }
    catch(const char* msg)
    {
        std::cout << msg <<std::endl;
    }

    /*
    std::cout << "L\n";
    Laplacian3D_ToeplitzMatrix L3d=Laplacian3D_ToeplitzMatrix(rows,cols,arrays);
    std::cout << "yess\n";
    try
    {
        Vectord b = L3d*vec;
        b=b;
    }
    catch(const char* msg)
    {
        std::cout << msg <<std::endl;
    }
    std::cout << "yes\n";
*//*
    int rows = 4;
    int cols = 3;
    Laplacian2D_ToeplitzMatrix L2d=Laplacian2D_ToeplitzMatrix(rows,cols);
    COO L2d_COO = L2d.COO_Laplacian();*/
    return 0;
}