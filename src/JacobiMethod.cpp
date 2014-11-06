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
bool verbose =false;

ZZ computeQ(vec_ZZ &a1, vec_ZZ &a2) {
    RR qNumerator   = a1 * a2;
    RR qDenominator = normsq(a2);
    RR division = qNumerator/qDenominator;

    return RoundToZZ(division);
}

void lagrange (vec_ZZ &a1, vec_ZZ &a2) {
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

bool lagrangeIT (mat_ZZ &g, mat_ZZ &z, int i, int j, RR &omega) {
    
    
    int s = i;
    int l = j;
    if(verbose) {
        cout << "? i,j " <<s<< " " << l<< endl;
        cout <<z <<endl;
        //cout <<g <<endl;
    }
    if (g(i, i) > g(j,j)){
        s = j;
        l = i;
    }

    RR gij = g(s,l);
    RR gss = g(s,s);
    RR gll = g(l,l);
    ZZ q = RoundToZZ(gij/gss);
    if(verbose) cout << "q " << q <<endl;
    if (abs(q) <= 1 && (((omega*omega)*gll) <= ( gss + gll - 2*(abs(gij))))) {
        if(verbose) cout << "no red" << endl;
        return false;
    }
    //cout << "old z[l] " << normsq(z(l))  << "old g[l] " <<normsq(g(l)) <<endl;
    z(l) -= q * z(s);
    g(l) -= q * g(s);
    //cout << "new z[l] " <<normsq(z(l) )<< "new g[l] " <<normsq(g(l)) <<endl;
    for (int k = 1; k <= g.NumRows(); k++ ) {
        g(k,l) = g(l,k);
    }

    g(l,l) -= q * g(l,s);
    assert(g == z * transpose(z));
    if(verbose) {
        cout << "new i,j " <<s<< " " << l<<endl;
        cout <<z <<endl;
        //cout <<g <<endl;
    }

    return true;
}

#pragma mark Fast Jacobi Method

bool fastJacobiMethodLoopShouldRun(mat_ZZ &g, RR omega) {
    int nRows = g.NumRows();

    for (int i = 1; i < nRows; i++) {
        for (int j = i + 1; j <= nRows; j++) {

            int s = i;
            int l = j;

            if (g(i, i) > g(j,j)){
                s = j;
                l = i;
            }

            RR gij = g(s,l);
            RR gss = g(s,s);
            RR gll = g(l,l);

            ZZ condition1 = abs(RoundToZZ(gij/gss));
            if (condition1 > 1) {
                //cout << "indexes : " << s << ", " << l << endl;
                //cout << "gij, gss, gll : " << gij << ", " << gss << ", " << gll << endl;
                return true;
            }

            if (((omega * omega)*gll) > ( gss + gll - 2*(abs(gij)))) {
                //cout << "indexes : " << s << ", " << l << endl;
                //cout << "gij, gss, gll : " << gij << ", " << gss << ", " << gll << endl;
                return true;
            }
        }
    }

    return false;
}

// Returns Z, the unimodular reduction matrix

mat_ZZ fastJacobiMethod(mat_ZZ &basis, RR omega) {
    int n = basis.NumRows();
    mat_ZZ g = basis * transpose(basis);
    bool didReplace = true;
    long count = 0;
    while(didReplace){
        count++;
        if(count%1000==0)cout << count<<" ";
        didReplace = false;
        for (int i = 1; i < n; i++) {
            for (int j = i + 1; j <= n; j++) {
                didReplace =  didReplace || lagrangeIT(g, basis, i, j, omega);
            }
        }

    }
    cout <<count << " while loop runs"<<endl;
    return basis;
}

#pragma mark Reduce Lattice

void JacobiMethod::reduceLattice (mat_ZZ &matrix, RR omega) {
    mat_ZZ reducedLattice = fastJacobiMethod(matrix, omega);
    matrix = reducedLattice;
}
