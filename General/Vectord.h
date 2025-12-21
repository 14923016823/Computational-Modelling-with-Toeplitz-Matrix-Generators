#pragma once
#include <vector>
#include <stdexcept>
#include <iostream>
//#include <tuple>



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
