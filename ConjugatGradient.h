#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "Matrix.h"
#include "Vectord.h"

// Classical Conjugate Gradient for SPD systems A x = b.
// Returns number of iterations performed (>=0) on success.
// Returns -1 on dimension mismatch or numerical breakdown (e.g., p^T A p ~= 0).
int ConjugateGradient(Matrix& A,
                      const Vectord& b,
                      Vectord& x,
                      int maxIters = 1000,
                      double relTol = 1e-8,
                      double absTol = 0.0);

#endif // CONJUGATE_GRADIENT_H
