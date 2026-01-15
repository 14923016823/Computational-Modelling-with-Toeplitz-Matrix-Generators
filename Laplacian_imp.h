#ifndef LAPLACIAN_IMP_H
#define LAPLACIAN_IMP_H

#include <tuple>
#include <utility>

#include "SparseToeplitz.h"
#include "BlockToeplitz.h"
#include "Vectord.h"

template<typename State, typename LinOp>  
inline State& MatVecMult(const LinOp& T, State& x) 
{
    x = T * x;
    return x;
}



template<typename State, typename ... LinOp>
struct LaplacianOp {
    std::tuple<const LinOp*...> ops;

    explicit LaplacianOp(const LinOp&... ops_) : ops(&ops_...) {}

    State operator*(const State& in) const {
        State out = in;
        std::apply([&](const auto*... op){ ((out = (*op) * out), ...); }, ops);
        return out;
    }
};

#endif