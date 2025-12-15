#ifndef VECTORD_H
#define VECTORD_H

#include <stdexcept>
#include <iostream>

class Vectord {
public:
    int Length;
    double* Vec;

    Vectord();
    explicit Vectord(int n);

    Vectord(const Vectord& other);
    Vectord(Vectord&& other) noexcept;

    Vectord& operator=(const Vectord& other);
    Vectord& operator=(Vectord&& other) noexcept;

    ~Vectord();

    int len() const;

    void PrintVector() const;

    double& operator[](int i);
    const double& operator[](int i) const;

    double* data();
    const double* data() const;

    void resize(int n);
    void fill(double value);

    double dot(const Vectord& other) const;

    double norm2() const;
    double norm() const;

    void scal(double a);
    void axpy(double a, const Vectord& x);
    void axpby(double a, double b, const Vectord& x);
};

#endif // VECTORD_H





