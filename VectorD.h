#pragma once

#include <initializer_list>
#include <memory>

class Vectord
{
public:
    int Length;
    double* Vec;
    Vectord(int length);

    Vectord(const std::initializer_list<double>& list);

    int len();

    void print();

    void Sum(Vectord VecIn);
};