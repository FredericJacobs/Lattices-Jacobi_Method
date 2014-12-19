#include <iostream>
#include <sstream>
#include <fstream>
#include <cstddef>
#include <newNTL/LLL.h>
#include "ReductionQualityChecker.h"
#include "JacobiMethod.h"
#include "MatrixFactory.h"
#include "macros.h"
#include <chrono>
using namespace std;
using namespace newNTL;
using namespace std::chrono;

#include "vec_double.h"

# pragma mark Utility Methods

const long BENCHMARK_MATRIX_SIZE   = 30;
const long BENCHMARK_MATRIX_BIT    = 10;

const int  BENCHMARK_DIMENSION_MIN = 50;
const int  BENCHMARK_DIMENSION_MAX = 200;

RR computeAverageRR (vec_RR &result){

    RR average = 0;

    for (int i = 1; i <= result.length() ; i++){
        average += result(i);
    }

    return average/result.length();
}

double computeAverageDouble(Vec<double> &result){

    RR average = 0;

    for (int i = 1; i <= result.length() ; i++){
        average += result(i);
    }

    return to_double(average/result.length());
}

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





double lllvsJacobiRandomLattice(int dimension){

    cout << "Starting reduction of dimension " << dimension << endl;

    ofstream jacobiDataFile, lllDataFile;

    jacobiDataFile.open ("dataRandomLattice_Jacobi.txt", ios::app);
    lllDataFile.open("dataRandomLattice_LLL.txt",        ios::app);

    int testsToRun = 50;

    vec_RR orthogonalityLLL;
    vec_RR orthogonalityJacobi;
    vec_RR hermiteLLL;
    vec_RR hermiteJacobi;

    vec_double timeJacobi;
    vec_double timeLLL;

    vec_double countsJacobi;

    orthogonalityJacobi.SetLength(testsToRun);
    orthogonalityLLL.SetLength(testsToRun);
    hermiteJacobi.SetLength(testsToRun);
    hermiteLLL.SetLength(testsToRun);
    timeJacobi.SetLength(testsToRun);
    timeLLL.SetLength(testsToRun);
    countsJacobi.SetLength(testsToRun);

    for (int i = 1; i <= testsToRun; i++){

        cout << "Running loop " << i << " for dimension " << dimension << endl;

        // - take an HNF basis of a random lattice

        mat_ZZ mat = MatrixFactory::makePrimeHNFMatrix(dimension, BENCHMARK_MATRIX_BIT);

        mat_ZZ lllReducedMat = mat;
        mat_ZZ jacobiReducedMat = mat;

        // - reduce the basis with LLL

        high_resolution_clock::time_point lllT1= high_resolution_clock::now();
        LLL_fplll(lllReducedMat);
        high_resolution_clock::time_point lllT2 = high_resolution_clock::now();

        cout << "LLL Reduced" << endl;

        double lllDuration = std::chrono::duration_cast<std::chrono::microseconds>(lllT2 - lllT1).count();
        timeLLL(i) = lllDuration;
        orthogonalityLLL(i) = orthogonalityDefect(lllReducedMat);
        hermiteLLL(i) = hermiteFactor(lllReducedMat);

        // - reduce the LLL-reduced basis and check if HF and OD are better <== I assume you mean with Jacobi?

        double omega  = 0.8;

        high_resolution_clock::time_point jacobiT1= high_resolution_clock::now();
        double count = JacobiMethod::reduceLattice(jacobiReducedMat, omega);
        high_resolution_clock::time_point jacobiT2 = high_resolution_clock::now();
        cout << "Jacobi reduced" << endl;

        double jacobiDuration = std::chrono::duration_cast<std::chrono::microseconds>(jacobiT2 - jacobiT1).count();

        orthogonalityJacobi(i)              = jacobiDuration;
        countsJacobi(i)                     = count;
        orthogonalityJacobi(i)              = orthogonalityDefect(jacobiReducedMat);
        hermiteJacobi(i)                    = hermiteFactor(jacobiReducedMat);
    }

    jacobiDataFile << dimension << " " << computeAverageRR(orthogonalityJacobi) << " " << computeAverageRR(hermiteJacobi) << " " << computeAverageDouble(timeJacobi) << " " << computeAverageDouble(countsJacobi) << endl;
    lllDataFile    << dimension << " " << computeAverageRR(orthogonalityLLL)    << " " << computeAverageRR(hermiteLLL)    << " " << computeAverageDouble(timeLLL)    << endl;
}

void benchmarkRandomLattice () {

    ofstream jacobiDataFile, lllDataFile;

    jacobiDataFile.open ("dataRandomLattice_Jacobi.txt");
    lllDataFile.open("dataRandomLattice_LLL.txt");

    for (int i = BENCHMARK_DIMENSION_MIN ; i < BENCHMARK_DIMENSION_MAX; i++) {
        lllvsJacobiRandomLattice(i);
    }

    jacobiDataFile.close();
    lllDataFile.close();
}

int main()
{
    benchmarkRandomLattice();
    return 0;
}
