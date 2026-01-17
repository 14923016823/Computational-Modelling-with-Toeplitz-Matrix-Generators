#ifndef RHS_H
#define RHS_H

#include "Vectord.h"

class RHS {
public:
    virtual ~RHS() = default;
    virtual void operator()(double t, const Vectord& u, Vectord& rhs) = 0;
};

#endif
