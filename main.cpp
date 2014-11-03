#include <iostream>
#include <sstream>
#include <newNTL/LLL.h>
#include <newNTL/mat_ZZ.h>
#include "JacobiMethod.h"

using namespace std;
using namespace newNTL;

#include "tools.h"

const long BENCHMARK_MATRIX_SIZE  = 30;
const long BENCHMARK_MATRIX_BIT   = 10;
const ZZ   BENCHMARK_MATRIX_SEED  = 0;


# pragma mark Utility Methods

mat_ZZ generateRandomLatticeBase(long n, long bit, ZZ seed) {
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


RR computeHermiteFactor(mat_ZZ &mat) {
    RR numerator   = sqrt(normsq(mat(1)));
    RR denominator = sqrt(determinant(mat * transpose(mat)));
    RR result =  numerator/denominator;
    return result;
}

RR computeOrthogonalityDefect(mat_ZZ &mat) {
    RR multResult = sqrt(normsq(mat(1)));

    for (int i = 2; i <= mat.NumRows(); i++) {
        multResult *= sqrt(normsq(mat(i)));
    }

    RR denominator = sqrt(determinant(mat * transpose(mat)));

    return pow(multResult/denominator, 1./mat.NumRows());
}


void compareHermite(RR oldHermite, RR newHermite) {



}

void compareDefects(RR oldDefect, RR newDefect) {
    if (oldDefect > newDefect) {
        cout << "The orthogonality defect was improved. We now have : " << newDefect << " used to be : " << oldDefect << endl;
    } else if (oldDefect < newDefect) {
        cout << "The orthogonality defect was worsened. This is probably a bug" << endl;
    } else{
        cout << "The orthogonality defect didn't improve." << endl;
    }
}

# pragma mark Reduction Tests

void lllThenJacobi() {
    mat_ZZ matrix = generateRandomLatticeBase(BENCHMARK_MATRIX_SIZE, BENCHMARK_MATRIX_BIT, BENCHMARK_MATRIX_SEED);
    LLL_fplll(matrix);

    RR defectLLL  = computeOrthogonalityDefect(matrix);
    RR hermiteLLL = computeHermiteFactor(matrix);

    JacobiMethod::reduceLattice(matrix);

    RR defectLLLJacobi  = computeOrthogonalityDefect(matrix);
    RR hermiteLLLJacobi = computeHermiteFactor(matrix);

    cout << "Reduced LLL then Jacobi" << endl;
    compareDefects(defectLLL, defectLLLJacobi);
}

void jacobiThenLLL () {
    mat_ZZ matrix = generateRandomLatticeBase(BENCHMARK_MATRIX_SIZE, BENCHMARK_MATRIX_BIT, BENCHMARK_MATRIX_SEED);

    JacobiMethod::reduceLattice(matrix);

    RR defectJacobi  = computeOrthogonalityDefect(matrix);
    RR hermiteJacobi = computeHermiteFactor(matrix);

    LLL_fplll(matrix);

    RR defectLLLJacobi  = computeOrthogonalityDefect(matrix);
    RR hermiteLLLJacobi = computeHermiteFactor(matrix);

    cout << "Reduced Jacobi then LLL" << endl;
    compareDefects(defectJacobi, defectLLLJacobi);
}

void jacobiVSLLL () {
    mat_ZZ matrix = generateRandomLatticeBase(BENCHMARK_MATRIX_SIZE, BENCHMARK_MATRIX_BIT, BENCHMARK_MATRIX_SEED);
    
    mat_ZZ matrixLLL = matrix;
    LLL_fplll(matrixLLL);
    RR defectLLL  = computeOrthogonalityDefect(matrixLLL);
    RR hermiteLLL = computeHermiteFactor(matrixLLL);
    cout << "defectLLL " << defectLLL <<" hermiteLLL " <<hermiteLLL <<endl;

    JacobiMethod::reduceLattice(matrix);
    RR defectJacobi  = computeOrthogonalityDefect(matrix);
    RR hermiteJacobi = computeHermiteFactor(matrix);
    cout << "defectJacobi " << defectJacobi <<" hermiteJacobi " <<hermiteJacobi  <<" norm b1 " <<sqrt(normsq(matrix(1)))<<endl;
    
    if (defectJacobi > defectLLL) {
        cout << "LLL has a better defect" << endl;
    } else if (defectLLL > defectJacobi){
        cout << "Jacobi has a better defect" << endl;
    } else{
        cout << "Jacobi and LLL have same defect" << endl;
    }
}

void findParameterForSameDefect () {

    RR defectJacobi;
    RR defectLLL;
    RR omega = 0.6;

    mat_ZZ matrix = generateRandomLatticeBase(BENCHMARK_MATRIX_SIZE, BENCHMARK_MATRIX_BIT, BENCHMARK_MATRIX_SEED);
    mat_ZZ matrixLLL = matrix;

    LLL_fplll(matrixLLL);
    defectLLL  = computeOrthogonalityDefect(matrixLLL);


    do {
        omega += 0.05;
        cout << "Trying with omega: " << omega << endl;

        mat_ZZ matrixToReduce = matrix;

        JacobiMethod::reduceLattice(matrixToReduce, omega);

    } while (defectJacobi > defectLLL);

    cout << "Jacobi defect is " << defectJacobi << " LLL defect is " << defectLLL << endl;
}


int main()
{
    mat_ZZ matrix =   generateRandomLatticeBase(40, 20, BENCHMARK_MATRIX_SEED);
    //matrix(40,1)=1;//add very short vector 1,0,...,0,0,1
    //cout << " matrix contains short vector 1,0, ...,0,1" << endl;

    mat_ZZ matrixJacobi = matrix;

    cout << "LLL, dim 40, seed 0" <<endl;
    LLL_fplll(matrix);
    
    RR defectLLL  = computeOrthogonalityDefect(matrix);
    RR hermiteLLL = computeHermiteFactor(matrix);
    cout << "defectLLL " << defectLLL <<" hermiteLLL " <<hermiteLLL <<" norm b1 " <<sqrt(normsq(matrix(1)))<<endl;
    
   matrix =   generateRandomLatticeBase(100, 20, BENCHMARK_MATRIX_SEED);
    cout << "LLL, dim 100, seed 0" <<endl;
    LLL_fplll(matrix);
    
    defectLLL  = computeOrthogonalityDefect(matrix);
    hermiteLLL = computeHermiteFactor(matrix);
    cout << "defectLLL " << defectLLL <<" hermiteLLL " <<hermiteLLL <<" norm b1 " <<sqrt(normsq(matrix(1)))<<endl;

    
    matrix =   generateRandomLatticeBase(500, 20, BENCHMARK_MATRIX_SEED);
    LLL_fplll(matrix);
    
    defectLLL  = computeOrthogonalityDefect(matrix);
     hermiteLLL = computeHermiteFactor(matrix);
    cout << "LLL, dim 500, seed 0" <<endl;
    cout << "defectLLL " << defectLLL <<" hermiteLLL " <<hermiteLLL <<" norm b1 " <<sqrt(normsq(matrix(1)))<<endl;
    
    
    cout << "Jacobi, dim 40, seed 0, matrix contains short vector 1,0, ...,0,1" << endl;
    JacobiMethod::reduceLattice(matrixJacobi);
    RR defectJacobi  = computeOrthogonalityDefect(matrixJacobi);
    RR hermiteJacobi = computeHermiteFactor(matrixJacobi);
    cout << "defectJacobi " << defectJacobi <<" hermiteJacobi " <<hermiteJacobi  <<" norm b1 " <<sqrt(normsq(matrixJacobi(1)))<<endl;

    return 0;
}
