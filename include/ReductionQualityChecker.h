#ifndef REDUCTIONQUALITYCHECKER_H
#define REDUCTIONQUALITYCHECKER_H

#include <newNTL/mat_ZZ.h>
#include <newNTL/matrix.h>
#include <newNTL/RR.h>

class ReductionQualityChecker
{
public:
    static newNTL::RR computeHermiteFactor(newNTL::mat_ZZ &mat);
    static newNTL::RR computeOrthogonalityDefect(newNTL::mat_ZZ &mat);
protected:
private:
};

#endif // MATRIXFACTORY_H


