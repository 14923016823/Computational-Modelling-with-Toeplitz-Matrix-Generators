#include <vector>
#include <stdexcept>
#include <iostream>
//#include <tuple>

#include "Vectord.h"


Vectord::Vectord(int length = 0) 
{
    Length=length;
    Vec=new double[Length];
}

int Vectord::len()
{
    return Length;
}

Vectord::~Vectord()
{
    delete Vec;
}

/*~Vectord()
{
    Delete
}*/

void Vectord::PrintVector()
{
    for (int i = 0; i < Length; ++i)
    {
        printf("vec[%d]=%f\n", i, Vec[i]);
    }
    printf("done\n");
}

void Vectord::Sum(const Vectord& VecIn)
{
    if (VecIn.Length != Length)
    {
        throw std::invalid_argument("You can't sum vectors with different sizes");
    }
    for (int i = 0; i < Length; i++)
    {
        Vec[i] += VecIn.Vec[i];
    }
}


