#ifndef JACOBIMETHOD_H
#define JACOBIMETHOD_H

#include <newNTL/mat_ZZ.h>

class JacobiMethod
{
    public:
        static void reduceLattice (newNTL::mat_ZZ &matrix);
    protected:
    private:
};

#endif // JACOBIMETHOD_H
