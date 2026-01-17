#ifndef VECTORD_H
#define VECTORD_H

#include <vector>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>

template<typename T = double>
class Vector {
public:
    Vector() = default;
    explicit Vector(int n) : data_(n, T(0)) { update_ptr(); }
    Vector(const T* data, int n) : data_(data, data + n) { update_ptr(); }
    Vector(const Vector& other) = default;
    Vector(Vector&& other) noexcept = default;
    Vector& operator=(const Vector& other) = default;
    Vector& operator=(Vector&& other) noexcept = default;
    ~Vector() = default;

    int len() const { return static_cast<int>(data_.size()); }

    void PrintVector() const {
        for (size_t i = 0; i < data_.size(); ++i) {
            std::cout << std::fixed << std::setprecision(3) << data_[i] << " ";
        }
        std::cout << std::endl;
    }

    T& operator[](int i) { return data_[i]; }
    const T& operator[](int i) const { return data_[i]; }

    T* Vec = nullptr;

    T* data() { return data_.data(); }
    const T* data() const { return data_.data(); }

    void resize(int n) { data_.resize(n, T(0)); update_ptr(); }
    void fill(T value) { std::fill(data_.begin(), data_.end(), value); update_ptr(); }

    T dot(const Vector& other) const {
        if (len() != other.len()) throw std::runtime_error("dot: size mismatch");
        return std::inner_product(data_.begin(), data_.end(), other.data_.begin(), T(0));
    }

    T norm2() const {
        return std::inner_product(data_.begin(), data_.end(), data_.begin(), T(0));
    }

    T norm() const { return std::sqrt(norm2()); }

    void scal(T a) {
        std::transform(data_.begin(), data_.end(), data_.begin(),
                      [a](T x) { return a * x; });
        update_ptr();
    }

    void axpy(T a, const Vector& x) {
        if (len() != x.len()) throw std::runtime_error("axpy: size mismatch");
        std::transform(data_.begin(), data_.end(), x.data_.begin(), data_.begin(),
                      [a](T y, T z) { return y + a * z; });
        update_ptr();
    }

    void axpby(T a, T b, const Vector& x) {
        if (len() != x.len()) throw std::runtime_error("axpby: size mismatch");
        std::transform(data_.begin(), data_.end(), x.data_.begin(), data_.begin(),
                      [a, b](T y, T z) { return a * y + b * z; });
        update_ptr();
    }

private:
    void update_ptr() { Vec = data_.data(); }
    std::vector<T> data_;
};

// Keep old name for backwards compatibility (global alias)
using Vectord = Vector<double>;

#endif





