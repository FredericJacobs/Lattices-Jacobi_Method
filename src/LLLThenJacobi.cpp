#include "macros.h"
#include "MatrixFactory.h"
#include "JacobiMethod.h"
#include <newNTL/LLL.h>
#include <newNTL/vector.h>
#include <newNTL/vec_double.h>

using namespace newNTL;
using namespace std;

#define DIMENSIONS   100
#define BITS_SIZE    30

// - take an HNF basis of a random lattice
// - reduce the basis with LLL
// - compute HF and orthogonality defect (OD)
// - reduce the LLL-reduced basis and check if HF and OD are better
// for the Hermite factor please take not necessarily the first basis
// vector but the shortest
//
// you should search for it within all n vectors in the reduced basis (or
// sort the basis vectors)


double testLLLReductionRandom (){

    // - take an HNF basis of a random lattice

    mat_ZZ mat = MatrixFactory::makePrimeHNFMatrix(DIMENSIONS, BITS_SIZE);

    mat_ZZ reducedMat = mat;

    // - reduce the basis with LLL

    LLL_fplll(reducedMat);

    // - compute HF and orthogonality defect (OD)

    RR orthogonalityDefect_LLL = orthogonalityDefect(reducedMat);
    RR hermiteFactor_LLL       = hermiteFactor(reducedMat);

    // - reduce the LLL-reduced basis and check if HF and OD are better <== I assume you mean with Jacobi?

    JacobiMethod::reduceLattice(reducedMat, 0.8);

    RR orthogonalityDefect_Jacobi = orthogonalityDefect(reducedMat);
    RR hermiteFactor_Jacobi       = hermiteFactor(reducedMat);


    cout << " LLL Defect:  " << orthogonalityDefect_LLL << " Jacobi Defect:  " << orthogonalityDefect_Jacobi << endl;
    cout << " LLL Hermite: " << hermiteFactor_LLL       << " Jacobi Hermite: " << hermiteFactor_Jacobi << endl;
}

double computeAverage (Vec<double> &result){

    RR average = 0;

    for (int i = 1; i <= result.length() ; i++){
        average += result(i);
    }

    return to_double(average/result.length());
}

int main() {

    int testsToRun = 200;

    Vec<double> results;
    results.SetLength(200);

    for (int i = 1; i <= testsToRun; i++){
        results (i) = testLLLReductionRandom();
    }

    cout << "Average improvement: " << computeAverage(results) << endl;

    return 0;
}