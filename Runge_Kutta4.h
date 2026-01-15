#ifndef RUNGE_KUTTA4_H
#define RUNGE_KUTTA4_H

#include <tuple>
#include <utility>

template<class RHS, class State>
void rk4_step(RHS&& rhs,
              State& u, double t, double dt,
              State& k1, State& k2, State& k3, State& k4, State& tmp)
{
    const double h  = dt;
    const double h2 = 0.5 * dt;

    rhs(u, k1, t);

    tmp = u;
    tmp.axpy(h2, k1);
    rhs(tmp, k2, t + h2);

    tmp = u;
    tmp.axpy(h2, k2);
    rhs(tmp, k3, t + h2);

    tmp = u;
    tmp.axpy(h, k3);       
    rhs(tmp, k4, t + h);

    u.axpy(h / 6.0, k1);
    u.axpy(h / 3.0, k2);
    u.axpy(h / 3.0, k3);
    u.axpy(h / 6.0, k4);
}

#endif