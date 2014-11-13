#include "ReductionQualityChecker.h"

using namespace std;
using namespace newNTL;

RR findShortestNormVector (mat_ZZ &mat) {
    ZZ  shortestVectNormSq = -1;

    for (int i = 1; i <= mat.NumRows(); i++){

        ZZ normSq = normsq(mat(i));

        if((shortestVectNormSq == -1) || shortestVectNormSq > normSq){
            shortestVectNormSq  = normSq;
        }
    }

    return sqrt(shortestVectNormSq);
}

RR ReductionQualityChecker::computeHermiteFactor(mat_ZZ &mat) {
    RR numerator   = findShortestNormVector(mat);
    RR denominator = pow(abs(determinant(mat)), (1./mat.NumRows()));
    RR result =  numerator/denominator;
    return result;
}

RR ReductionQualityChecker::computeOrthogonalityDefect(mat_ZZ &mat) {
    RR multResult = sqrt(normsq(mat(1)));

    for (int i = 2; i <= mat.NumRows(); i++) {
        multResult *= sqrt(normsq(mat(i)));
    }

    RR denominator = sqrt(determinant(mat * transpose(mat)));

    return pow(multResult/denominator, 1./mat.NumRows());
}
