#ifndef JACOBIMETHOD_H
#define JACOBIMETHOD_H

#include <newNTL/mat_ZZ.h>

class JacobiMethod
{
    public:
        static newNTL::mat_ZZ reduceLattice (newNTL::mat_ZZ &matrix);
        JacobiMethod();
        virtual ~JacobiMethod();
    protected:
    private:
};

#endif // JACOBIMETHOD_H
