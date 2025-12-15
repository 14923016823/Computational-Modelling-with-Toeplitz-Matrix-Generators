#include <vector>
#include <stdexcept>
#include <iostream>
//#include <tuple>

#include "Matrix.h"

int Matrix::rows() const { return Num_Rows; }
int Matrix::cols() const { return Num_Cols; }

/*Matrix::Matrix(int cols, int rows)
{
    Num_Cols=cols;
    Num_Rows=rows;
}*/

Matrix::Matrix()
{
    Num_Cols=0;
    Num_Rows=0;
}
