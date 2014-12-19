#include "JacobiMethod.h"
#include <cstddef>
using namespace std;
using namespace newNTL;

#include <newNTL/LLL.h>
#include <newNTL/RR.h>
#include <newNTL/ZZ.h>
#include <newNTL/mat_RR.h>
#include <newNTL/mat_ZZ.h>
#include <newNTL/vec_ZZ.h>
#include <newNTL/vector.h>
#include "ReductionQualityChecker.h"

#pragma mark Lagrange Algorithm
bool verbose = true;

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
    
    if(verbose) {
        assert(g == z * transpose(z));
        cout << "new i,j " <<s<< " " << l<<endl;
        cout <<z <<endl;
        //cout <<g <<endl;
    }
    //cout <<" . " ;
    return true;
}


// Returns Z, the unimodular reduction matrix

long fastJacobiMethod(mat_ZZ &basis, RR omega) {
        int n = basis.NumRows();
    //cout << "dim "<<n  <<endl;
    mat_ZZ g = basis * transpose(basis);
    bool didReplace = true;
    long count = 0;
    RR prevortho=0;
    while(didReplace){
        count++;
        //if(count%1000==0)cout << count<<" ";
        didReplace = false;
        for (int i = 1; i < n; i++) {
            for (int j = i + 1; j <= n; j++) {
                didReplace =  didReplace || lagrangeIT(g, basis, i, j, omega);
            }
        }
        /*if(!(count%10)){
            RR currentortho = ReductionQualityChecker::computeOrthogonalityDefect(basis);
            cout <<  "  " <<count << " "  << prevortho - currentortho <<endl;
            prevortho = currentortho;
        }*/
    }
    cout <<"end " << count << " loops "<<  ReductionQualityChecker::computeOrthogonalityDefect(basis) << " defect " <<endl;
    return count;
}


#pragma mark Double Implementation

bool doublelagrangeIT (Mat<double> &g, Mat<double> &z, int i, int j, double &omega) {
    int s = i;
    int l = j;

    if (g(i, i) > g(j,j)){
        s = j;
        l = i;
    }

    double gij = g(s,l);
    double gss = g(s,s);
    double gll = g(l,l);
    double q      = rint(gij/gss);

    if (abs(q) <= 1 && (((omega*omega)*gll) <= ( gss + gll - 2*(abs(gij))))) {
        return false;
    }
    z(l) -= q * z(s);
    g(l) -= q * g(s);
    for (int k = 1; k <= g.NumRows(); k++ ) {
        g(k,l) = g(l,k);
    }

    g(l,l) -= q * g(l,s);
    if(g != z * transpose(z)){
        cerr<< "g " <<  g << endl;
        cerr<<  "z z^t " << z *transpose(z) << endl;
        exit(0);
    }

    return true;
}

Mat<double> fastJacobiMethod(Mat<double> &basis, double omega) {
    int n = basis.NumRows();
    Mat<double> g = basis * transpose(basis);
    bool didReplace = true;
    long count = 0;
    while(didReplace){
        count++;
        didReplace = false;
        for (int i = 1; i < n; i++) {
            for (int j = i + 1; j <= n; j++) {
                didReplace =  didReplace || doublelagrangeIT(g, basis, i, j, omega);
            }
        }
       
    }
    cout <<count << " while loop rounds"<<endl;
    return basis;
}

#pragma mark Reduce Lattice

long JacobiMethod::reduceLattice (mat_ZZ &matrix, RR omega) {
    return fastJacobiMethod(matrix, omega);
}

void JacobiMethod::reduceLatticeDouble(newNTL::Mat<double> &matrix, double omega) {
    Mat<double> reducedLattice = fastJacobiMethod(matrix, omega);
    matrix = reducedLattice;
}
