#include "MatrixFactory.h"
#include <newNTL/LLL.h>
#include <newNTL/mat_ZZ.h>
#include <newNTL/RR.h>
#include <iostream>
#include <sstream>
#include "JacobiMethod.h"

using namespace std;
using namespace newNTL;

#include "tools.h"

ZZ getSeed(){
    ZZ seed = rand();
    return seed;
}

mat_ZZ MatrixFactory::makeHNFMatrix(long n, long bit) {

    vec_ZZ v;
    generate_random_HNF(v,n,bit, getSeed());
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

vec_ZZ getRandomVector(long size, long bit) {
    vec_ZZ randomVec;
    randomVec.SetLength(size);
    SetSeed(getSeed());

    for (int i = 1; i <= size; i++) {
        randomVec(i) = RandomBits_ZZ(bit);
    }

    return randomVec;
}

mat_ZZ MatrixFactory::makeRandomSquareMatrix(long n, long bit) {
    mat_ZZ randomMatrix;
    randomMatrix.SetDims(n,n);

    do {
        SetSeed(getSeed());

        for (int i = 1; i <= n; i++) {
            randomMatrix(i) = getRandomVector(n, bit);
        }
    } while (determinant(randomMatrix) == 0);


    return randomMatrix;
}
