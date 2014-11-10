#ifndef MATRIXFACTORY_H
#define MATRIXFACTORY_H

#include <newNTL/mat_ZZ.h>
#include <newNTL/matrix.h>

class MatrixFactory
{
public:
    static newNTL::mat_ZZ      makeHNFMatrix(long n, long bit);
    static newNTL::Mat<double> makeRandomSquareMatrixDouble(long n, long bit);
    static newNTL::mat_ZZ      makeRandomSquareMatrixZZ(long n, long bit);
protected:
private:
};

#endif // MATRIXFACTORY_H
