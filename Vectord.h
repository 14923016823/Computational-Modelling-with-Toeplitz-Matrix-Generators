#ifndef VECTORD_H
#define VECTORD_H

#include <Eigen/Dense>
#include <iostream>
#include <stdexcept>

template<typename T = double>
class Vector {
public:
    Vector() = default;
    explicit Vector(int n) { data_.resize(n); data_.setZero(); update_ptr(); }
    Vector(const T* data, int n) { data_ = Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>(data, n); update_ptr(); }
    Vector(const Vector& other) = default;
    Vector(Vector&& other) noexcept = default;
    Vector& operator=(const Vector& other) = default;
    Vector& operator=(Vector&& other) noexcept = default;
    ~Vector() = default;

    int len() const { return static_cast<int>(data_.size()); }

    void PrintVector() const {
        std::cout << data_.transpose() << std::endl;
    }

    T& operator[](int i) { return data_(i); }
    const T& operator[](int i) const { return data_(i); }

    T* Vec = nullptr;

    T* data() { return data_.data(); }
    const T* data() const { return data_.data(); }

    void resize(int n) { data_.resize(n); data_.setZero(); update_ptr(); }
    void fill(T value) { data_.setConstant(value); update_ptr(); }

    T dot(const Vector& other) const {
        if (len() != other.len()) throw std::runtime_error("dot: size mismatch");
        return data_.dot(other.data_);
    }

    T norm2() const { return data_.squaredNorm(); }
    T norm() const { return data_.norm(); }

    void scal(T a) {
        data_ *= a;
        update_ptr();
    }

    void axpy(T a, const Vector& x) {
        if (len() != x.len()) throw std::runtime_error("axpy: size mismatch");
        data_ += a * x.data_;
        update_ptr();
    }

    void axpby(T a, T b, const Vector& x) {
        if (len() != x.len()) throw std::runtime_error("axpby: size mismatch");
        data_ = a * data_ + b * x.data_;
        update_ptr();
    }

private:
    void update_ptr() { Vec = data_.data(); }
    Eigen::Matrix<T, Eigen::Dynamic, 1> data_;
};

// Keep old name for backwards compatibility (global alias)
using Vectord = Vector<double>;

#endif





