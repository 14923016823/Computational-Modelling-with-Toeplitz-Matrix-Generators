#ifndef RUNGE_KUTTA4_H
#define RUNGE_KUTTA4_H

#include "Vectord.h"
#include "RHS.h"

class Runge_Kutta4 {
private:
    int N;
    Vectord k1, k2, k3, k4, tmp;

public:
    Runge_Kutta4(int size) : N(size), k1(size), k2(size), k3(size), k4(size), tmp(size) {}

    void step(double t, double dt, Vectord& u, RHS& rhs) {
        const double h = dt;
        const double h2 = 0.5 * dt;
        const double h6 = dt / 6.0;
        const double h3 = dt / 3.0;

        // k1 = f(t, u)
        rhs(t, u, k1);

        // tmp = u + h2*k1
        tmp = u;
        tmp.axpy(h2, k1);
        // k2 = f(t + h2, tmp)
        rhs(t + h2, tmp, k2);

        // tmp = u + h2*k2
        tmp = u;
        tmp.axpy(h2, k2);
        // k3 = f(t + h2, tmp)
        rhs(t + h2, tmp, k3);

        // tmp = u + h*k3
        tmp = u;
        tmp.axpy(h, k3);
        // k4 = f(t + h, tmp)
        rhs(t + h, tmp, k4);

        // u = u + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        u.axpy(h6, k1);
        u.axpy(h3, k2);
        u.axpy(h3, k3);
        u.axpy(h6, k4);
    }
};

#endif