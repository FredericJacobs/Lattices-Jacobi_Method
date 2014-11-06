#include <iostream>
#include <sstream>
#include <chrono>
#include <fstream>
#include <newNTL/LLL.h>
#include <newNTL/mat_ZZ.h>
#include "JacobiMethod.h"
#include "MatrixFactory.h"


using namespace std;
using namespace newNTL;
using namespace std::chrono;


# pragma mark Utility Methods
const long BENCHMARK_MATRIX_BIT    = 10;
const long BENCHMARK_MATRIX_SIZE    = 30;
const int  BENCHMARK_DIMENSION_MIN = 20;
const int  BENCHMARK_DIMENSION_MAX = 400;


RR computeHermiteFactor(mat_ZZ &mat) {
    RR numerator   = sqrt(normsq(mat(1)));
    RR denominator = pow(abs(determinant(mat)), (1./mat.NumRows()));
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
    mat_ZZ matrix = MatrixFactory::makeHNFMatrix(BENCHMARK_MATRIX_SIZE, BENCHMARK_MATRIX_BIT);
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
    mat_ZZ matrix = MatrixFactory::makeHNFMatrix(BENCHMARK_MATRIX_SIZE, BENCHMARK_MATRIX_BIT);

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
    mat_ZZ matrix = MatrixFactory::makeHNFMatrix(BENCHMARK_MATRIX_SIZE, BENCHMARK_MATRIX_BIT);

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
    RR omega = 0.9;

    mat_ZZ matrix = MatrixFactory::makeRandomSquareMatrix(100, 50);
    mat_ZZ matrixLLL = matrix;

    LLL_fplll(matrixLLL);
    defectLLL  = computeOrthogonalityDefect(matrixLLL);

    mat_ZZ matrixToReduce = matrix;

    do {
        omega += 0.05;
        cout << "Trying with omega: " << omega << endl;

        JacobiMethod::reduceLattice(matrixToReduce, omega);

        defectJacobi = computeOrthogonalityDefect(matrixToReduce);

    } while (defectJacobi > defectLLL);

    cout << "Jacobi defect is " << defectJacobi << " LLL defect is " << defectLLL << endl;
}

// This is the first graph of the Fast Jacobi paper, it benchmarks the hermite factors of FJ reduced basis with LLL's.
// QUESTION: HOW MANY BITS? NOT MENTIONED IN THE PAPER?





void generateHermiteDataRandomMatrix () {

    ofstream jacobiDataFile, lllDataFile;

    jacobiDataFile.open ("dataHermiteRandomMatrix_Jacobi.txt");
    lllDataFile.open("dataHermiteRandomMatrix_LLL.txt");

    for (int i = BENCHMARK_DIMENSION_MIN ; i < BENCHMARK_DIMENSION_MAX; i++) {
        mat_ZZ randomMatrix  = MatrixFactory::makeRandomSquareMatrix(i, BENCHMARK_MATRIX_BIT);


        /// Jacobi-Reduction
        mat_ZZ jacobiReduced = randomMatrix;

        high_resolution_clock::time_point jacobiT1 = high_resolution_clock::now();
        JacobiMethod::reduceLattice(jacobiReduced, 0.99);
        high_resolution_clock::time_point jacobiT2 = high_resolution_clock::now();
        auto jacobiDuration = std::chrono::duration_cast<std::chrono::microseconds>(jacobiT2 - jacobiT1).count();
        jacobiDataFile << i << " " << computeHermiteFactor(jacobiReduced) << " " << computeOrthogonalityDefect(jacobiReduced) << " " << jacobiDuration << endl;

        // LLL-reduction
        mat_ZZ lllReduced = randomMatrix;
        high_resolution_clock::time_point lllT1 = high_resolution_clock::now();
        LLL_fplll(lllReduced, 0.99);
        high_resolution_clock::time_point lllT2 = high_resolution_clock::now();
        auto lllDuration = std::chrono::duration_cast<std::chrono::microseconds>(lllT2-lllT1).count();
        lllDataFile << i << " " << computeHermiteFactor(lllReduced) << " " << computeOrthogonalityDefect(lllReduced) << " " << lllDuration << endl;
    }

    jacobiDataFile.close();
    lllDataFile.close();
}

int main()
{
    //generateHermiteDataRandomMatrix();

    //LLL_fplll(randomMatrix);



    mat_ZZ matrix = MatrixFactory::makeHNFMatrix(40, 20);
    //matrix(40,1)=1;//add very short vector 1,0,...,0,0,1
    //cout << " matrix contains short vector 1,0, ...,0,1" << endl;

    mat_ZZ matrixJacobi = matrix;

    cout << "LLL, dim 40, seed 0" <<endl;
    LLL_fplll(matrix);

    RR defectLLL  = computeOrthogonalityDefect(matrix);
    RR hermiteLLL = computeHermiteFactor(matrix);
    cout << "defectLLL " << defectLLL <<" hermiteLLL " <<hermiteLLL <<" norm b1 " <<sqrt(normsq(matrix(1)))<<endl;
//
//   matrix = MatrixFactory::makeHNFMatrix(100, 20);
//    cout << "LLL, dim 100, seed 0" <<endl;
//    LLL_fplll(matrix);
//
//    defectLLL  = computeOrthogonalityDefect(matrix);
//    hermiteLLL = computeHermiteFactor(matrix);
//    cout << "defectLLL " << defectLLL <<" hermiteLLL " <<hermiteLLL <<" norm b1 " <<sqrt(normsq(matrix(1)))<<endl;
//
//
//    matrix = MatrixFactory::makeHNFMatrix(500, 20);
//    LLL_fplll(matrix);
//
//    defectLLL  = computeOrthogonalityDefect(matrix);
//     hermiteLLL = computeHermiteFactor(matrix);
//    cout << "LLL, dim 500, seed 0" <<endl;
//    cout << "defectLLL " << defectLLL <<" hermiteLLL " <<hermiteLLL <<" norm b1 " <<sqrt(normsq(matrix(1)))<<endl;
//
//
//    cout << "Jacobi, dim 40, seed 0, matrix contains short vector 1,0, ...,0,1" << endl;
//    JacobiMethod::reduceLattice(matrixJacobi);
//    RR defectJacobi  = computeOrthogonalityDefect(matrixJacobi);
//    RR hermiteJacobi = computeHermiteFactor(matrixJacobi);
//    cout << "defectJacobi " << defectJacobi <<" hermiteJacobi " <<hermiteJacobi  <<" norm b1 " <<sqrt(normsq(matrixJacobi(1)))<<endl;


    return 0;
}
