#include "ConjugateGradient.h"

#include <cmath>

static inline double safe_sqrt(double v) { return (v > 0.0) ? std::sqrt(v) : 0.0; }

int ConjugateGradient(Matrix& A,
                      const Vectord& b,
                      Vectord& x,
                      int maxIters,
                      double relTol,
                      double absTol)
{
    // Dimension checks
    if (A.rows() != b.len()) return -1;
    if (A.cols() <= 0) return -1;

    // If x is empty or wrong size, resize and set x=0
    if (x.len() != A.cols()) {
        x.resize(A.cols());
        x.fill(0.0);
    }

    // r = b - A*x
    Vectord r = b;
    Vectord Ax = A * x;
    r.axpy(-1.0, Ax);

    double b2 = b.dot(b);
    double bnorm = safe_sqrt(b2);

    double r2 = r.dot(r);
    double rnorm = safe_sqrt(r2);

    // stopping threshold: max(absTol, relTol*||b||)
    double thresh = absTol;
    double relPart = relTol * bnorm;
    if (relPart > thresh) thresh = relPart;

    if (rnorm <= thresh) return 0;

    Vectord p = r;
    double rsold = r2;

    for (int k = 0; k < maxIters; ++k)
    {
        Vectord Ap = A * p;

        double pAp = p.dot(Ap);
        if (std::fabs(pAp) < 1e-30) {
            // Numerical breakdown
            return -1;
        }

        double alpha = rsold / pAp;

        x.axpy(alpha, p);       // x = x + alpha p
        r.axpy(-alpha, Ap);     // r = r - alpha A p

        double rsnew = r.dot(r);
        double rnew_norm = safe_sqrt(rsnew);

        if (rnew_norm <= thresh) return k + 1;

        double beta = rsnew / rsold;

        // p = r + beta p
        p.scal(beta);
        p.axpy(1.0, r);

        rsold = rsnew;
    }

    return maxIters;
}
