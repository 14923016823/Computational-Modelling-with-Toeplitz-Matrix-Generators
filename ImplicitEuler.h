#ifndef IMPLICIT_EULER_H
#define IMPLICIT_EULER_H

#include "Matrix.h"
#include "Vectord.h"
#include "ConjugateGradient.h"
#include <iostream>

class ImplicitOperator : public Matrix {
private:
    Matrix* A_ptr;
    double dt;
public:
    ImplicitOperator(Matrix& a, double d) : A_ptr(&a), dt(d) {
        Num_Rows = a.rows();
        Num_Cols = a.cols();
    }

    void matvec(const Vectord& in, Vectord& out) const override {
        A_ptr->matvec(in, out);
        out.scal(-dt);
        out.axpy(1.0, in);
    }

    void operator*=(double c) override {
        dt *= c;
    }

    Matrix* Kronecker(Matrix& B) override {
        return nullptr; // Not implemented
    }

    Matrix* Clone(double c) override {
        return new ImplicitOperator(*A_ptr, dt * c);
    }
};

class ImplicitEuler {
private:
    int N;
    Vectord u_new;
    ImplicitOperator* op;
public:
    ImplicitEuler(int size) : N(size), u_new(size), op(nullptr) {}

    ~ImplicitEuler() { if (op) delete op; }

    void step(double t, double dt, Vectord& u, Matrix& A) {
        if (op == nullptr) {
            op = new ImplicitOperator(A, dt);
        } else if (op->dt != dt || op->A_ptr != &A) {
            delete op;
            op = new ImplicitOperator(A, dt);
        }
        // Solve op * u_new = u
        int iters = ConjugateGradient(*op, u, u_new);
        if (iters < 0) {
            std::cerr << "Conjugate Gradient failed to converge" << std::endl;
            // For now, copy u to u_new or something, but let's assume it works
        }
        u = u_new;
    }
};

#endif