#pragma once

#include "COO.h"
#include "CSR.h"
#include "DiagonalMatrix.h"
#include <experimental/random>
#include <random>

SparseToeplitz r_ST(int n_r, int n_c, int prob);

BlockCSR r_BCSR(int n_r, int n_c);

BlockCOO r_BCOO(int n_r, int n_c);

BlockToeplitz r_BT(int n_r, int n_c, int prob);

DiagonalMatrix r_Diag(int n_r, int n_c);

Vectord r_Vec(int n);