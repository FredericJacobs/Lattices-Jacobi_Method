#ifndef MATRIXFACTORY_H
#define MATRIXFACTORY_H

#include <newNTL/LLL.h>
#include <newNTL/mat_ZZ.h>

class MatrixFactory
{
public:
    static newNTL::mat_ZZ makeHNFMatrix(long n, long bit);
    static newNTL::mat_ZZ makeRandomSquareMatrix(long n, long bit);
protected:
private:
};

#endif // MATRIXFACTORY_H
