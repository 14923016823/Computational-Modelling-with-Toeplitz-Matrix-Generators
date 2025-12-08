#pragma once

class Vectord
{
protected:
    int Length;
public:
    double* Vec;
    Vectord(int length);
    ~Vectord();
   
    int len();
    void PrintVector();
    void Sum(const Vectord& VecIn);
};
