#include "Vectord.h"

#include <algorithm>
#include <cmath>
#include <iostream>

Vectord::Vectord()
    : Length(0), Vec(nullptr)
{}

Vectord::Vectord(int n)
    : Length(n), Vec(n > 0 ? new double[n] : nullptr)
{
    std::fill(Vec, Vec + Length, 0.0);
}

Vectord::Vectord(const Vectord& other)
    : Length(other.Length), Vec(other.Length > 0 ? new double[other.Length] : nullptr)
{
    std::copy(other.Vec, other.Vec + Length, Vec);
}

Vectord::Vectord(Vectord&& other) noexcept
    : Length(other.Length), Vec(other.Vec)
{
    other.Length = 0;
    other.Vec = nullptr;
}

Vectord& Vectord::operator=(const Vectord& other)
{
    if (this == &other) return *this;

    resize(other.Length);
    std::copy(other.Vec, other.Vec + Length, Vec);
    return *this;
}

Vectord& Vectord::operator=(Vectord&& other) noexcept
{
    if (this == &other) return *this;

    delete[] Vec;
    Length = other.Length;
    Vec = other.Vec;

    other.Length = 0;
    other.Vec = nullptr;
    return *this;
}

Vectord::~Vectord()
{
    delete[] Vec; 
}

int Vectord::len() const
{
    return Length;
}

void Vectord::PrintVector() const
{
    for (int i = 0; i < Length; ++i)
    {
        printf("vec[%d]=%f\n", i, Vec[i]);
    }
    printf("done\n");
}

double& Vectord::operator[](int i)
{
    return Vec[i];
}

const double& Vectord::operator[](int i) const
{

    return Vec[i];
}

double* Vectord::data()
{
    return Vec;
}

const double* Vectord::data() const
{
    return Vec;
}

void Vectord::resize(int n)
{
    if (n == Length) return;

    delete[] Vec;
    Length = n;
    Vec = (n > 0) ? new double[n] : nullptr;

    if(Vec)
        std::fill(Vec, Vec + Length, 0.0);
}

void Vectord::fill(double value)
{
    if (Vec)
        std::fill(Vec, Vec + Length, 0.0);
}

double Vectord::dot(const Vectord& other) const
{
    if (Length != other.Length)
        throw std::runtime_error("Vectord::dot size mismatch");

    double s = 0.0;
    for (int i = 0; i < Length; ++i)
        s += Vec[i] * other.Vec[i];

    return s;
}

double Vectord::norm2() const
{
    return this->dot(*this);
}

double Vectord::norm() const
{
    return std::sqrt(norm2());
}

void Vectord::scal(double a)
{
    for (int i = 0; i < Length; ++i)
        Vec[i] *= a;
}

void Vectord::axpy(double a, const Vectord& x)
{
    if (Length != x.Length)
        throw std::runtime_error("Vectord::axpy size mismatch");

    for (int i = 0; i < Length; ++i)
        Vec[i] += a * x.Vec[i];
}

void Vectord::axpby(double a, double b, const Vectord& x)
{
    if (Length != x.Length)
        throw std::runtime_error("Vectord::axpby size mismatch");

    for (int i = 0; i < Length; ++i)
        Vec[i] = a * Vec[i] + b * x.Vec[i];
}



