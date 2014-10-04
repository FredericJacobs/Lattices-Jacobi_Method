//compile as g++ test.cpp -o test -lnewntl -lfplll -lmpfr -lgmp -lm
#include <newNTL/LLL.h>
#include <newNTL/mat_ZZ.h>
#include "newskel.h"

using namespace std;
using namespace newNTL;

int main() {
    long n = 4;
    mat_ZZ M; M.SetDims(n, n);
    
    cout << " provide matrix of dimension " << n << endl;
    mat_ZZ T; T.SetDims(2,2);
    T[0][0] = 1;T[0][1] = 2;T[1][0] = 3;T[1][1] = 4;
    cout << " Example: [ [1 2] [ 3 4] ] becomes matrix" << endl; cout << T << endl;
    cin >> M ;
    //use BKZ from fplll lib
    LLLBKZ_fplll(M);
    //use LLL from fplll lib
    LLL_fplll(M);
    cout <<"reduced matrix:" << endl;
    cout<< M << endl;
    return 0;
}
