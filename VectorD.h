#pragma once

#include <initializer_list>
#include <memory>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <exception>
#include <vector>
#include <iomanip>
#include <chrono>
#include <numeric>

class Vectord
{
public:
    int Length;
    double* Vec;
    Vectord(int length);

    Vectord(const std::initializer_list<double>& list);

    int len() const;

    void print();

    void PrintVector();

    void Sum(Vectord VecIn);

    double dot(Vectord& VecIn);

    Vectord& scal(double c);
};