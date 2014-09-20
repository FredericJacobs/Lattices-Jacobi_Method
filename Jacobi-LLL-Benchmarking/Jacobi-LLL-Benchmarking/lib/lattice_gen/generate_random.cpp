#include <sstream>
#include <iostream>
#include <NTL/LLL.h>

NTL_CLIENT;

#include "tools.h"

mat_ZZ genSVPChallengeLattice(){
    long n = 80;
    long bit = 10;
    ZZ seed; seed = 0;

    vec_ZZ v; generate_random_HNF(v,n,bit,seed);
    mat_ZZ B; B.SetDims(n,n); clear(B);
    B(1,1) = v(1);
    for (int i=2; i<=n; i++)
    {
	B(i,1)=v(i);
	B(i,i)=1;
    }
    return B;
}
