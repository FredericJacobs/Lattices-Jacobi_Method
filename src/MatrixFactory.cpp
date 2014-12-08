
#include <cstddef>
#include "MatrixFactory.h"
#include <newNTL/LLL.h>
#include <newNTL/mat_ZZ.h>
#include <newNTL/RR.h>
#include <iostream>
#include <sstream>
#include "JacobiMethod.h"
#include <time.h>

using namespace std;
using namespace newNTL;

#include "tools.h"

ZZ getSeed(){
    srandom(time(NULL));
    ZZ seed = rand();
    return seed;
}

RR computeDeterminant(Mat<double> &matrix){
    mat_RR matrixRR;
    long cols = matrixRR.NumCols();
    long rows = matrixRR.NumRows();
    matrixRR.SetDims(cols, rows);

    for (long i = 1; i <= cols; i++){
        for (long j = 1; j <= rows; j++){
            matrixRR(i,j)=matrix(i,j);
        }
    }

    return determinant(matrixRR);
}

mat_ZZ MatrixFactory::makeRandomHNFMatrix(long n, long bit){
    SetSeed(getSeed());

    mat_ZZ randomMatrix;
    randomMatrix.SetDims(n, n);

    for (int i = 1; i <= n; i++){
        for (int j = i; j <= n; j++ ){
            randomMatrix(j,i) = RandomBits_ZZ(bit);
        }
    }

    return randomMatrix;
}

mat_ZZ MatrixFactory::makePrimeHNFMatrix(long n, long bit){
    ZZ seed = getSeed();
    vec_ZZ v;
    generate_random_HNF(v,n,bit,seed);
    mat_ZZ B;
    B.SetDims(n,n);
    clear(B);

    B(1,1) = v(1);
    for (int i=2; i<=n; i++) {
        B(i,1)=v(i);
        B(i,i)=1;
    }
    return B;
}

vec_ZZ getRandomVectorZZ(long size, long bit) {
    vec_ZZ randomVec;
    randomVec.SetLength(size);
    SetSeed(getSeed());

    for (int i = 1; i <= size; i++) {
        randomVec(i) = RandomBits_ZZ(bit);
    }

    return randomVec;
}

Mat<double> MatrixFactory::makeRandomSquareMatrixDouble(long n, long bit) {
    ZZ base = 2;
    ZZ maxDouble;
    conv(maxDouble, DBL_MAX);
    assert(power(base, bit) <  maxDouble);

    Mat<double> randomMatrix;
    randomMatrix.SetDims(n,n);
    do {
        SetSeed(getSeed());

        for (int i = 1; i <= n; i++)
            conv(randomMatrix(i),getRandomVectorZZ(n, bit));

    } while (computeDeterminant(randomMatrix) == 0);

    return randomMatrix;
}

mat_ZZ MatrixFactory::makeRandomSquareMatrixZZ(long n, long bit) {
    mat_ZZ randomMatrix;
    randomMatrix.SetDims(n,n);

    do {
        for (int i = 1; i <= n; i++) {
            randomMatrix(i) = getRandomVectorZZ(n, bit);
        }

    } while (determinant(randomMatrix) == 0);

    return randomMatrix;
}
