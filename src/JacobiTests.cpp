#include <iostream>
#include <sstream>
#include <chrono>
#include <fstream>
#include <newNTL/LLL.h>
#include <ReductionQualityChecker.h>
#include "JacobiMethod.h"
#include "MatrixFactory.h"

#include "macros.h"

using namespace std;
using namespace newNTL;
using namespace std::chrono;

# pragma mark Utility Methods

const long BENCHMARK_MATRIX_SIZE   = 30;
const long BENCHMARK_MATRIX_BIT    = 10;

const int  BENCHMARK_DIMENSION_MIN = 20;
const int  BENCHMARK_DIMENSION_MAX = 400;

mat_ZZ toMatZZ(Mat<double> matDouble){
    mat_ZZ matrixZZ;
    long cols = matrixZZ.NumCols();
    long rows = matrixZZ.NumRows();
    matrixZZ.SetDims(rows, cols);

    for (long i = 1; i <= cols; i++){
        for (long j = 1; j <= rows; j++){
            matrixZZ(i,j)= matDouble(i,j);
        }
    }
    return matrixZZ;
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
    mat_ZZ matrix = MatrixFactory::makePrimeHNFMatrix(BENCHMARK_MATRIX_SIZE, BENCHMARK_MATRIX_BIT);
    LLL_fplll(matrix);

    RR defectLLL  = orthogonalityDefect(matrix);
    RR hermiteLLL = hermiteFactor(matrix);

    JacobiMethod::reduceLattice(matrix);

    RR defectLLLJacobi  = orthogonalityDefect(matrix);
    RR hermiteLLLJacobi = hermiteFactor(matrix);

    cout << "Reduced LLL then Jacobi" << endl;
    compareDefects(defectLLL, defectLLLJacobi);
}

void jacobiThenLLL () {
    mat_ZZ matrix = MatrixFactory::makePrimeHNFMatrix(BENCHMARK_MATRIX_SIZE, BENCHMARK_MATRIX_BIT);

    JacobiMethod::reduceLattice(matrix);

    RR defectJacobi  = orthogonalityDefect(matrix);
    RR hermiteJacobi = hermiteFactor(matrix);

    LLL_fplll(matrix);

    RR defectLLLJacobi  = orthogonalityDefect(matrix);
    RR hermiteLLLJacobi = hermiteFactor(matrix);

    cout << "Reduced Jacobi then LLL" << endl;
    compareDefects(defectJacobi, defectLLLJacobi);
}

void jacobiVSLLL () {
    mat_ZZ matrix = MatrixFactory::makePrimeHNFMatrix(BENCHMARK_MATRIX_SIZE, BENCHMARK_MATRIX_BIT);

    mat_ZZ matrixLLL = matrix;
    LLL_fplll(matrixLLL);
    RR defectLLL  = orthogonalityDefect(matrixLLL);
    RR hermiteLLL = hermiteFactor(matrixLLL);
    cout << "defectLLL " << defectLLL <<" hermiteLLL " <<hermiteLLL <<endl;

    JacobiMethod::reduceLattice(matrix);
    RR defectJacobi  = orthogonalityDefect(matrix);
    RR hermiteJacobi = hermiteFactor(matrix);
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

    mat_ZZ matrix = MatrixFactory::makeRandomSquareMatrixZZ(100, 50);
    mat_ZZ matrixLLL = matrix;

    LLL_fplll(matrixLLL);
    defectLLL  = orthogonalityDefect(matrixLLL);

    mat_ZZ matrixToReduce = matrix;

    do {
        omega += 0.05;
        cout << "Trying with omega: " << omega << endl;

        JacobiMethod::reduceLattice(matrixToReduce, omega);

        defectJacobi = orthogonalityDefect(matrixToReduce);

    } while (defectJacobi > defectLLL);

    cout << "Jacobi defect is " << defectJacobi << " LLL defect is " << defectLLL << endl;
}

// This is the first graph of the Fast Jacobi paper, it benchmarks the hermite factors of FJ reduced basis with LLL's.
// QUESTION: HOW MANY BITS? NOT MENTIONED IN THE PAPER?


void generateHermiteDataRandomMatrix () {

    ofstream jacobiDataFile, lllDataFile;

    jacobiDataFile.open ("dataHermiteRandomMatrix_Jacobi.txt");
    lllDataFile.open("dataHermiteRandomMatrix_LLL.txt");


    for (int i = BENCHMARK_DIMENSION_MIN ; i < BENCHMARK_DIMENSION_MAX; i=i+10) {
        mat_ZZ randomMatrix  = MatrixFactory::makeRandomSquareMatrixZZ(i, BENCHMARK_MATRIX_BIT);
        cout << "matrix gener" <<endl;

        /// Jacobi-Reduction
        mat_ZZ jacobiReduced = randomMatrix;

        high_resolution_clock::time_point jacobiT1 = high_resolution_clock::now();
        JacobiMethod::reduceLattice(jacobiReduced, 0.99);

        high_resolution_clock::time_point jacobiT2 = high_resolution_clock::now();
        auto jacobiDuration = std::chrono::duration_cast<std::chrono::microseconds>(jacobiT2 - jacobiT1).count();

        jacobiDataFile << i << " " << hermiteFactor(jacobiReduced) << " " << orthogonalityDefect(jacobiReduced) << " " << jacobiDuration << endl;
    }

    jacobiDataFile.close();
    lllDataFile.close();
}

void generateHermiteDoubleDataRandomMatrix () {

    ofstream jacobiDataFile, lllDataFile;

    jacobiDataFile.open ("dataHermiteRandomMatrix_Jacobi.txt");
    lllDataFile.open("dataHermiteRandomMatrix_LLL.txt");

    for (int i = BENCHMARK_DIMENSION_MIN ; i < BENCHMARK_DIMENSION_MAX; i++) {
        Mat<double> randomMatrix  = MatrixFactory::makeRandomSquareMatrixDouble(i, BENCHMARK_MATRIX_BIT);

        /// Jacobi-Reduction
        Mat<double> jacobiReduced = randomMatrix;

        high_resolution_clock::time_point jacobiT1 = high_resolution_clock::now();
        JacobiMethod::reduceLatticeDouble(jacobiReduced, 0.99);
        high_resolution_clock::time_point jacobiT2 = high_resolution_clock::now();

        
        auto jacobiDuration = std::chrono::duration_cast<std::chrono::microseconds>(jacobiT2 - jacobiT1).count();
        mat_ZZ jacobiReducedZZ;
        conv(jacobiReducedZZ,jacobiReduced);

        jacobiDataFile << i << " " << hermiteFactor(jacobiReducedZZ) << " " << orthogonalityDefect(jacobiReducedZZ) << " " << jacobiDuration << endl;

        // LLL-reduction
        mat_ZZ lllReduced;
        conv(lllReduced, randomMatrix);
        
        high_resolution_clock::time_point lllT1 = high_resolution_clock::now();
        LLL_fplll(lllReduced, 0.99);
        high_resolution_clock::time_point lllT2 = high_resolution_clock::now();
        auto lllDuration = std::chrono::duration_cast<std::chrono::microseconds>(lllT2-lllT1).count();
        lllDataFile << i << " " << hermiteFactor(lllReduced) << " " << orthogonalityDefect(lllReduced) << " " << lllDuration << endl;
        cout << "i = "  << i <<" ended" <<endl;
    }

    jacobiDataFile.close();
    lllDataFile.close();
}

int main()
{
    generateHermiteDoubleDataRandomMatrix();
    return 0;
}
