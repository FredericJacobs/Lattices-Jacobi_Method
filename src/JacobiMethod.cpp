#include "JacobiMethod.h"

using namespace std;
using namespace newNTL;

#include <newNTL/LLL.h>
#include <newNTL/RR.h>
#include <newNTL/ZZ.h>
#include <newNTL/mat_RR.h>
#include <newNTL/mat_ZZ.h>
#include <newNTL/vec_ZZ.h>
#include <newNTL/vector.h>

#pragma mark Lagrange Algorithm

ZZ computeQ(vec_ZZ &a1, vec_ZZ &a2) {
    RR qNumerator   = a1 * a2;
    RR qDenominator = normsq(a2);
    RR division = qNumerator/qDenominator;

    return RoundToZZ(division);
}

void lagrange (vec_ZZ &a1, vec_ZZ &a2) {
    // In the paper, they compare the square root of the norms, here we will just compare the norms, which should be equivalent since they are both positive integers.
    if ( normsq(a1) < normsq(a2) ){
        swap(a1, a2);
    }

    do {
        ZZ q = computeQ(a1, a2);
        vec_ZZ a2_prev = a2;
        a2 = a1 - (q*a2);
        a1 = a2_prev;
    } while (!(normsq(a1) <= normsq(a2)));
}

#pragma mark Generic Jacobi Method

bool genericJacobiMethodLoopShouldRun(mat_ZZ &matrix) { // Checking it this way is sub-optimal, could do a lazier approach of only evaluating when required.
    long nRows = matrix.NumRows();

    for (long i = 1; i < nRows; i++) {
        for (long j = i + 1; j <= nRows; j++) {

            if (!(normsq(matrix(i)) <= normsq(matrix(j)))) {
                return true;
            }

            if (!(abs(matrix(i) * matrix(j)) <= (1/2)*(normsq(matrix(i))))) {
                return true;
            }
        }
    }

    return false;
}

mat_ZZ genericJacobiMethod(mat_ZZ &matrix) {
    long n = matrix.NumRows();

    while (genericJacobiMethodLoopShouldRun(matrix)) {
        for (long i = 1; i < n; i++) {
            for (long j = i + 1; j <= n; j++){
                lagrange(matrix(i), matrix(j));
            }
        }
    }

    return matrix;
}

#pragma mark LagrangeIT Algorithm

void lagrangeIT (mat_ZZ &g, mat_ZZ &z, int i, int j) {
    int s = i;
    int l = j;

    if (g(i, i) > g(j,j)){
        s = j;
        l = i;
    }

    RR qNum   = g(j,i);
    RR qDenom = g(s, s);
    ZZ q = (RoundToZZ(qNum/qDenom));

    mat_ZZ zij = ident_mat_ZZ(g.NumCols());
    zij(s,l)= -q;

    g = zij * transpose(g) * transpose(zij);
    z = z*zij;
}

#pragma mark Fast Jacobi Method

bool fastJacobiMethodLoopShouldRun(mat_ZZ &g, RR &omega) {
    int nRows = g.NumRows();

    for (int i = 1; i < nRows; i++) {
        for (int j = i + 1; j <= nRows; j++) {

            if (!(abs(RoundToZZ( (g(i) * g(j))/(g(i,i)) )) <= 1)) { // ss == ii?
                return true;
            }

            if (!(((omega*omega)*g(j,j)) <= ( g(i,i) + g(j,j) - 2*(abs(g(i,j)))))) { // ll == jj?
                return true;
            }
        }
    }

    return false;
}

// Returns Z, the unimodular reduction matrix

mat_ZZ fastJacobiMethod(mat_ZZ &matrix, mat_ZZ &a, mat_ZZ &z, RR omega) {
    int n = matrix.NumRows();
    mat_ZZ g = a * transpose(a);
    z = ident_mat_ZZ(matrix.NumCols());

    while (fastJacobiMethodLoopShouldRun(matrix, omega)) {
        for (int i = 1; i < n; i++) {
            for (int j = i + 1; j <= n; j++) {
                lagrangeIT(g, z, i, j);
            }
        }
    }

    return z;
}

#pragma mark Reduce Lattice

void JacobiMethod::reduceLattice (mat_ZZ &matrix) {
    mat_ZZ a;
    mat_ZZ z;

    fastJacobiMethod(matrix, a, z, 0.6);

    matrix = a*z;
}
