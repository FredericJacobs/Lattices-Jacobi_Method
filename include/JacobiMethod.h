#ifndef JACOBIMETHOD_H
#define JACOBIMETHOD_H

#include <newNTL/mat_ZZ.h>
#include <newNTL/RR.h>

class JacobiMethod
{
    public:
        static void reduceLatticeDouble (newNTL::Mat<double> &matrix,double omega = 0.9);
        static long reduceLattice (newNTL::mat_ZZ &matrix, newNTL::RR omega = 0.9);
    protected:
    private:
};

#endif // JACOBIMETHOD_H
