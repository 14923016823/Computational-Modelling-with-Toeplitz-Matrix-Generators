#pragma once

#include <initializer_list>
#include <memory>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <exception>

class Vectord
{
public:
    int Length;
    double* Vec;
    Vectord(int length);

    Vectord(const std::initializer_list<double>& list);

    int len();

    void print();

    void PrintVector();

    void Sum(Vectord VecIn);
};