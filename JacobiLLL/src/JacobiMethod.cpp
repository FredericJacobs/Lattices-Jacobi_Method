#include "JacobiMethod.h"

using namespace std;
using namespace newNTL;
#include <newNTL/mat_RR.h>
#include <newNTL/mat_ZZ.h>
#include <newNTL/vec_ZZ.h>
#include <newNTL/vector.h>

// Questions for Anja :
// 1) No easy accessor for columns?

JacobiMethod::JacobiMethod()
{
    //ctor
}

/// LAGRANGE ALGORITHM

ZZ multiplyVectorsAsMatrices(vec_ZZ &a1, vec_ZZ &a2){ // Does ((a_i)^T*a_j) where a_i and a_j are vectors
    mat_ZZ mat_a1; mat_a1.SetDims(1, a1.length());
    mat_ZZ mat_a2; mat_a2.SetDims(a2.length(), 1);
    return (mat_a1*mat_a2)(1,1);
}

ZZ computeQ(vec_ZZ &a1, vec_ZZ &a2) {
    RR qNumerator   = multiplyVectorsAsMatrices(a1, a2); // <== Hardcoded transpose(a1)*a2
    RR qDenominator = normsq(a2);
    RR division = qNumerator/qDenominator;

    return RoundToZZ(division);
}

mat_ZZ matrixWithA1A2(vec_ZZ &a1, vec_ZZ &a2) {
    mat_ZZ a1a2; a1a2.SetDims(2, a1.length()); a1a2(1) = a1; a1a2(2) = a2;
    return transpose(a1a2);
}

mat_ZZ lagrangeAlgorithm (vec_ZZ &a1, vec_ZZ &a2)
{
    // In the paper, they compare the square root of the norms, here we will just compare the norms, which should be equivalent since they are both positive integers.
    if ( normsq(a1) < normsq(a2) ){
        swap(a1, a2);
    }

    mat_ZZ sample_Zmatrix; sample_Zmatrix.SetDims(2,2); sample_Zmatrix(1,1) = 0; sample_Zmatrix(2,1) = 1; sample_Zmatrix(1,2)=1;

    do {
        ZZ q = computeQ(a1, a2);
        mat_ZZ M = sample_Zmatrix;
        M(2,2)= q;
        M = matrixWithA1A2(a1, a2)* M;
        a1 = transpose(M)(1); a2 = transpose(M)(2);
    } while (!(normsq(a1) <= normsq(a2)));

    return matrixWithA1A2(a1, a2);
}

/// Generic Jacobi Method

bool genericJacobiMethodLoopShouldRun(mat_ZZ &matrix){

    // Rows accessors are easier than column accessors.

    mat_ZZ transposedM = transpose(matrix);

    for (int i = 1; i < matrix.NumRows()-1; i++){


        if (normsq(transposedM(i)) <= normsq(transposedM(i+1))){
            return true;
        }

        for (int j = i; j < matrix.NumRows()-1; i++){
            if (abs(multiplyVectorsAsMatrices(transposedM(j), transposedM(i))) <= (normsq(transposedM(i)))){
                return true;
            }
        }
    }

    return false;
}

mat_ZZ genericJacobiMethod(mat_ZZ &matrix) {

    int n = matrix.NumCols();

    while (genericJacobiMethodLoopShouldRun(matrix)){
        for (int i = 1; i < n-1; i++){
            for (int j = i + 1; j < n; j++){
                lagrangeAlgorithm(transpose(matrix)(i), transpose(matrix)(j));
            }
        }
    }

    return matrix;
}

mat_ZZ JacobiMethod::reduceLattice (mat_ZZ &matrix){
    vec_ZZ a1; vec_ZZ a2;
    a1.SetLength(2); a2.SetLength(2);
    a1(1) = 2; a1(2) = 2; a2(1) = 0; a2(2) = 0;

    mat_ZZ matrix2 = lagrangeAlgorithm(a2, a1);

    return matrix;
}

