#include <iostream>
#include <exception>
#include <cmath>
#include <initializer_list>

#include "VectorD.h"

Vectord::Vectord(int length)
{
    Length = length;
    Vec = new double[Length];
}

Vectord::Vectord(const std::initializer_list<double>& list)
: Vectord((int)list.size())
{
    std::uninitialized_copy(list.begin(), list.end(), Vec);
}

int Vectord::len()
{
    return Length;
}

void Vectord::print()
{
    std::cout << "[";
    for(int i=0;i<Length;i++)
    {
        std::cout << Vec[i];
        if(i<Length-1)
        std::cout << ", ";
    }
    std::cout << ']' << std::endl;
}

void Vectord::Sum(Vectord VecIn)
{
    if(VecIn.Length=!Length)\
    {
        throw std::invalid_argument("You can't sum vectors with different sizes");
    }
    for(int i=0;i<Length;i++)
    {
        Vec[i]+=VecIn.Vec[i];
    }   
}