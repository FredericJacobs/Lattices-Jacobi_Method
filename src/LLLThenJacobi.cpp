
#include <cstddef>
#include "macros.h"
#include "MatrixFactory.h"
#include "JacobiMethod.h"
#include <fstream>

#include <iostream>
#include <sstream>
#include <chrono>
#include <newNTL/LLL.h>
#include <newNTL/vector.h>
#include <newNTL/vec_double.h>

using namespace newNTL;
using namespace std;
using namespace chrono;

#define DIMENSIONS   100
#define BITS_SIZE    10


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


double testLLLReductionRandomDimension (int dimension){

    cout << "Starting reduction of dimension " << dimension << endl;

    ofstream jacobiDataFile, lllDataFile;

    jacobiDataFile.open ("jacobiReductionData.txt", ios::app);
    lllDataFile.open("lllReductionData.txt",    ios::app);

    int testsToRun = 75;

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

        mat_ZZ mat = MatrixFactory::makePrimeHNFMatrix(dimension, BITS_SIZE);

        mat_ZZ reducedMat = mat;

        // - reduce the basis with LLL

        high_resolution_clock::time_point lllT1= high_resolution_clock::now();
        LLL_fplll(reducedMat);
        high_resolution_clock::time_point lllT2 = high_resolution_clock::now();

        double lllDuration = std::chrono::duration_cast<std::chrono::microseconds>(lllT2 - lllT1).count();
        timeLLL(i) = lllDuration;

        // - compute HF and orthogonality defect (OD)

        orthogonalityLLL(i) = orthogonalityDefect(reducedMat);
        hermiteLLL(i) = hermiteFactor(reducedMat);

        // - reduce the LLL-reduced basis and check if HF and OD are better <== I assume you mean with Jacobi?


        double omega_min  = 0.75;
        double omega_max  = 0.95;
        double iterations = (omega_max*100) - (omega_min*100);

        vec_RR JorthogonalityDefectsJacobi;
        vec_RR JhermiteFactorsJacobi;
        vec_double JreductionTime;
        vec_double JjacobiIterations;

        JorthogonalityDefectsJacobi.SetLength(iterations);
        JhermiteFactorsJacobi.SetLength(iterations);
        JreductionTime.SetLength(iterations);
        JjacobiIterations.SetLength(iterations);

        mat_ZZ lllReducedMat = reducedMat;

        for (double j = 1; j <= iterations; j++){
            double omega = omega_min+(j/100);

            cout << "Running reduction on jacobi with omega " << omega << endl;

            mat_ZZ reductionMat = lllReducedMat;

            high_resolution_clock::time_point jacobiT1= high_resolution_clock::now();
            double count = JacobiMethod::reduceLattice(reductionMat, omega);
            high_resolution_clock::time_point jacobiT2 = high_resolution_clock::now();

            double jacobiDuration = std::chrono::duration_cast<std::chrono::microseconds>(jacobiT2 - jacobiT1).count();

            JreductionTime(j)              = jacobiDuration;
            JjacobiIterations(j)           = count;
            JorthogonalityDefectsJacobi(j) = orthogonalityDefect(reductionMat);
            JhermiteFactorsJacobi(j)       = hermiteFactor(reductionMat);
        }

        orthogonalityJacobi(i) = computeAverageRR(JorthogonalityDefectsJacobi);
        hermiteJacobi(i)       = computeAverageRR(JhermiteFactorsJacobi);
        timeJacobi(i)          = computeAverageDouble(JreductionTime);
        countsJacobi(i)        = computeAverageDouble(JjacobiIterations);
    }

    jacobiDataFile << dimension << " " << computeAverageRR(orthogonalityJacobi) << " " << computeAverageRR(hermiteJacobi) << " " << computeAverageDouble(timeJacobi) << " " << computeAverageDouble(countsJacobi) << endl;
    lllDataFile    << dimension << " " << computeAverageRR(orthogonalityLLL)    << " " << computeAverageRR(hermiteLLL)    << " " << computeAverageDouble(timeLLL)    << endl;
}


// - take an HNF basis of a random lattice
// - reduce the basis with LLL
// - compute HF and orthogonality defect (OD)
// - reduce the LLL-reduced basis and check if HF and OD are better
// for the Hermite factor please take not necessarily the first basis
// vector but the shortest
//
// you should search for it within all n vectors in the reduced basis (or
// sort the basis vectors)


double testLLLReductionRandom (double omega){

    ofstream jacobiDataFile, lllDataFile;

    jacobiDataFile.open ("jacobiReductionData.txt", ios::app);
    lllDataFile.open("lllReductionData.txt",    ios::app);

    int testsToRun = 100;

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

        // - take an HNF basis of a random lattice

        mat_ZZ mat = MatrixFactory::makePrimeHNFMatrix(DIMENSIONS, BITS_SIZE);

        mat_ZZ reducedMat = mat;

        // - reduce the basis with LLL

        high_resolution_clock::time_point lllT1= high_resolution_clock::now();
        LLL_fplll(reducedMat);
        high_resolution_clock::time_point lllT2 = high_resolution_clock::now();

        double lllDuration = std::chrono::duration_cast<std::chrono::microseconds>(lllT2 - lllT1).count();
        timeLLL(i) = lllDuration;

        // - compute HF and orthogonality defect (OD)

        orthogonalityLLL(i) = orthogonalityDefect(reducedMat);
        hermiteLLL(i) = hermiteFactor(reducedMat);

        // - reduce the LLL-reduced basis and check if HF and OD are better <== I assume you mean with Jacobi?


        high_resolution_clock::time_point jacobiT1= high_resolution_clock::now();
        double count = JacobiMethod::reduceLattice(reducedMat, omega);
        high_resolution_clock::time_point jacobiT2 = high_resolution_clock::now();

        double jacobiDuration = std::chrono::duration_cast<std::chrono::microseconds>(jacobiT2 - jacobiT1).count();
        timeJacobi(i) = jacobiDuration;

        countsJacobi(i) = count;

        orthogonalityJacobi(i) = orthogonalityDefect(reducedMat);
        hermiteJacobi(i)       = hermiteFactor(reducedMat);

//        cout << " LLL    Defect: "  << orthogonalityDefect_LLL          << " LLL    Hermite: "  << hermiteFactor_LLL    << " Duration: " << lllDuration    << endl;
//        cout << " Jacobi Defect: "  << orthogonalityDefect_Jacobi       << " Jacobi Hermite: "  << hermiteFactor_Jacobi << " Duration: " << jacobiDuration << endl;

    }

    jacobiDataFile << omega << " " << computeAverageRR(orthogonalityJacobi) << " " << computeAverageRR(hermiteJacobi) << " " << computeAverageDouble(timeJacobi) << " " << computeAverageDouble(countsJacobi) << endl;
    lllDataFile    << omega << " " << computeAverageRR(orthogonalityLLL)    << " " << computeAverageRR(hermiteLLL)    << " " << computeAverageDouble(timeLLL)    << endl;
}


double LLLJacobiOrtho(int dimension){
    
    cout << "Starting reduction of dimension " << dimension << endl;
    
    ofstream jacobiDataFile, lllDataFile;
    
    jacobiDataFile.open ("jacobiReductionData.txt", ios::app);
    lllDataFile.open("lllReductionData.txt",    ios::app);
    
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
        
        mat_ZZ mat = MatrixFactory::makePrimeHNFMatrix(dimension, BITS_SIZE);
        
        mat_ZZ reducedMat = mat;
        
        // - reduce the basis with LLL
        
        high_resolution_clock::time_point lllT1= high_resolution_clock::now();
        LLL_fplll(reducedMat);
        high_resolution_clock::time_point lllT2 = high_resolution_clock::now();
        
        double lllDuration = std::chrono::duration_cast<std::chrono::microseconds>(lllT2 - lllT1).count();
        timeLLL(i) = lllDuration;
        
        // - compute HF and orthogonality defect (OD)
        
        orthogonalityLLL(i) = orthogonalityDefect(reducedMat);
        hermiteLLL(i) = hermiteFactor(reducedMat);
        
        // - reduce the LLL-reduced basis and check if HF and OD are better <== I assume you mean with Jacobi?
        
        
        double omega_min  = 0.98;
        double omega_max  = 0.999;
        double iterations = (omega_max*100) - (omega_min*100);
        
        vec_RR JorthogonalityDefectsJacobi;
        vec_RR JhermiteFactorsJacobi;
        vec_double JreductionTime;
        vec_double JjacobiIterations;
        
        JorthogonalityDefectsJacobi.SetLength(iterations);
        JhermiteFactorsJacobi.SetLength(iterations);
        JreductionTime.SetLength(iterations);
        JjacobiIterations.SetLength(iterations);
        
        mat_ZZ lllReducedMat = reducedMat;
        
        for (double j = 1; j <= iterations; j++){
            double omega = omega_min+(j/100);
            
            cout << "Running reduction on jacobi with omega " << omega << endl;
            
            mat_ZZ reductionMat = lllReducedMat;
            
            high_resolution_clock::time_point jacobiT1= high_resolution_clock::now();
            double count = JacobiMethod::reduceLattice(reductionMat, omega);
            high_resolution_clock::time_point jacobiT2 = high_resolution_clock::now();
            
            double jacobiDuration = std::chrono::duration_cast<std::chrono::microseconds>(jacobiT2 - jacobiT1).count();
            
            JreductionTime(j)              = jacobiDuration;
            JjacobiIterations(j)           = count;
            JorthogonalityDefectsJacobi(j) = orthogonalityDefect(reductionMat);
            JhermiteFactorsJacobi(j)       = hermiteFactor(reductionMat);
        }
        
        orthogonalityJacobi(i) = computeAverageRR(JorthogonalityDefectsJacobi);
        hermiteJacobi(i)       = computeAverageRR(JhermiteFactorsJacobi);
        timeJacobi(i)          = computeAverageDouble(JreductionTime);
        countsJacobi(i)        = computeAverageDouble(JjacobiIterations);
    }
    
    jacobiDataFile << dimension << " " << computeAverageRR(orthogonalityJacobi) << " " << computeAverageRR(hermiteJacobi) << " " << computeAverageDouble(timeJacobi) << " " << computeAverageDouble(countsJacobi) << endl;
    lllDataFile    << dimension << " " << computeAverageRR(orthogonalityLLL)    << " " << computeAverageRR(hermiteLLL)    << " " << computeAverageDouble(timeLLL)    << endl;
}



int main() {

    ofstream jacobiDataFile, lllDataFile;

    jacobiDataFile.open ("jacobiReductionData.txt");
    lllDataFile.open("lllReductionData.txt");

    jacobiDataFile.close();
    lllDataFile.close();

    for (int dimension = 100; dimension <= 100; dimension += 5){
        LLLJacobiOrtho(dimension);
    }



    return 0;
}