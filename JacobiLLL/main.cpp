#include <iostream>
#include <sstream>
#include <newNTL/LLL.h>
#include <newNTL/mat_ZZ.h>

using namespace std;
using namespace newNTL;

#include "../lib/tools.h"

mat_ZZ generateRandomLatticeBase (long n, long bit, ZZ seed) {
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


int main()
{
    cout << generateRandomLatticeBase(80, 10, 0);
    return 0;
}
