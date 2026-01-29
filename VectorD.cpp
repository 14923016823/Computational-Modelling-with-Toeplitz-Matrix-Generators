#include "VectorD.h"

Vectord::Vectord(int length)
{
    Length = length;
    Vec = new double[Length];
    for(int i=0;i<Length;i++)
    {
        Vec[i]=0;
    }
}

double& Vectord::operator[](const int i)
{
    if(i<0 || i>=Length)
    {
        throw std::invalid_argument("trying to acces outside of vector");
    }
    return Vec[i];
}

double Vectord::operator[](const int i) const {
    if(i<0 || i>=Length)
    {
        throw std::invalid_argument("trying to acces outside of vector");
    }
    return Vec[i]; // returns by value, read-only context
}


Vectord::Vectord(const std::initializer_list<double>& list)
: Vectord((int)list.size())
{
    std::uninitialized_copy(list.begin(), list.end(), Vec);
}

int Vectord::len() const
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

void Vectord::PrintVector()
{
    for (int i = 0; i < Length; ++i)
    {
        printf("vec[%d]=%f\n", i, Vec[i]);
    }
    printf("done\n");
}

void Vectord::Sum(Vectord VecIn)
{
    if(VecIn.Length!=Length)
    {
        throw std::invalid_argument("You can't sum vectors with different sizes");
    }
    //int i;
    //#pragma omp for private(i)
    for(int i=0;i<Length;i++)
    {
        Vec[i]+=VecIn.Vec[i];
    }   
}

double Vectord::dot(const Vectord& VecIn) const
{
    if(VecIn.Length!=Length)
    {
        throw std::invalid_argument("You can't dot vectors with different sizes");
    }
    double s=0;
    for(int i=0;i<Length;i++)
    {
        s+=Vec[i]*VecIn.Vec[i];
    }
    return s;   
}

Vectord& Vectord::scal(double c)
{
    for(int i=0;i<Length;i++)
    {
        Vec[i]*=c;
    }
    return (*this);
}

void Vectord::axpy(const double  c,const Vectord& VecIn)
{
    if(VecIn.Length!=Length)
    {
        throw std::invalid_argument("You can't sum vectors with different sizes");
    }

    for(int i=0;i<Length;i++)
    {
        Vec[i]+=c*VecIn.Vec[i];
    }
}