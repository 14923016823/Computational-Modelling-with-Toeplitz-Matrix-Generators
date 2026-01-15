#ifndef RHS_H
#define RHS_H

#include <tuple>
#include "Vectord.h"

template<typename Func, typename... Args>
struct RHS {
    Func func;
    std::tuple<Args...> args;

    RHS(Func f, Args... a) : func(std::move(f)), args(std::move(a)...) {}

    template<class State>
    void operator()(const State& u, State& dudt, double t) const {
        std::apply([&](const Args&... unpacked_args) {
            func(u, dudt, t, unpacked_args...);
        }, args);
    }
};

#endif
