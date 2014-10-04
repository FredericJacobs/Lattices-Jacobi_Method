#include <newNTL/LLL.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include <fplll.h>
#include <stdint.h>
#include <cassert>

#define PLANTE(errormsg) { \
cerr << errormsg << "at line " << __LINE__ << endl; \
cout << errormsg << "at line " << __LINE__ << endl; \
abort(); }


namespace newNTL {
    void LLLBKZ_fplll(newNTL::mat_ZZ& B, double beta=0.99, int blocksize=20, double eta=0.51) {
        ZZ_mat<mpz_t> M;
        std::ostringstream os; os << B;
        std::istringstream is(os.str()); is>>M;
        fplll::lllReduction(M,beta,eta,fplll::LM_WRAPPER);
        blocksize = B.NumRows() - 10;
        int n = B.NumRows();
        if (blocksize>60) blocksize=60;
        BKZParam bkzparam;
        bkzparam.b=&M;
        bkzparam.flags=BKZ_NO_LLL|BKZ_MAX_TIME;
        bkzparam.blockSize=20;
        bkzparam.maxTime=30*blocksize;
        if (blocksize < 20) blocksize = 20;
        for (int block= 20; block<=blocksize; block++)  
        {
            cerr << "BKZ-" << block << endl;
            bkzparam.blockSize=block;
            bkzparam.maxTime=3*block;
            if (block == 45) {
                cerr << "using a semi-linear pruning" << endl;
                bkzparam.pruning.resize(n);
                for (int i=0; i<n; i++) {
                    bkzparam.pruning[i]=1-double(i)/(2.*n);
                }
            }
            if (block == 55) {
                cerr << "using a quasi-linear pruning" << endl;
                bkzparam.pruning.resize(n);
                for (int i=0; i<n; i++) {
                    bkzparam.pruning[i]=1-double(i)/double(1.1*n);
                }
            }
            fplll::bkzReduction(bkzparam);
        }
        std::ostringstream os2; os2<<M;
        std::istringstream is2(os2.str()); is2>>B;
    }
}

/**
 * Int24 definition
 */
const int INT24_MAX = 8388607;

class Int24
{
    private:
        unsigned char m_Internal[3];
    public:
        inline void setval(int input) {
            m_Internal[0]   = ((unsigned char*)&input)[0];
            m_Internal[1]   = ((unsigned char*)&input)[1];
            m_Internal[2]   = ((unsigned char*)&input)[2];
        }
        inline int getval() const {
            if ( m_Internal[2] & 0x80 ) // Is this a negative?  Then we need to siingn extend.
            {
                return (0xff << 24) | (m_Internal[2] << 16) | (m_Internal[1] << 8) | (m_Internal[0] << 0);
            }
            else
            {
                return (m_Internal[2] << 16) | (m_Internal[1] << 8) | (m_Internal[0] << 0);
            }
        }
        inline void copy(const Int24& input) {
            m_Internal[0]   = input.m_Internal[0];
            m_Internal[1]   = input.m_Internal[1];
            m_Internal[2]   = input.m_Internal[2];
        }

    public:
        Int24() {}
        Int24( const int val ) { setval(val); }
        Int24( const Int24& val ) { copy(val); }
        operator int() const { return getval(); }
        operator long() const { return getval(); }
        operator double() const { return getval(); }
        Int24& operator =( const Int24& input ) { copy(input); return *this; }
        Int24& operator =( const int input ) { setval(input); return *this; }

        /***********************************************/

        operator bool() const { return getval() != 0; }
        bool operator !() const { return !(getval()); }
        Int24 operator -() { return Int24(-getval()); }

        bool operator ==( const Int24& val ) const { return getval() == val.getval(); }
        bool operator !=( const Int24& val ) const { return getval() != val.getval(); }
        bool operator >=( const Int24& val ) const { return getval() >= val.getval(); }
        bool operator <=( const Int24& val ) const { return getval() <= val.getval(); }
        bool operator > ( const Int24& val ) const { return getval() >  val.getval(); }
        bool operator < ( const Int24& val ) const { return getval() <  val.getval(); }

        Int24 operator +( const Int24& val ) const {return Int24( getval() + val.getval() );}
        Int24 operator -( const Int24& val ) const {return Int24( getval() - val.getval() );}
        Int24 operator *( const Int24& val ) const {return Int24( getval() * val.getval() );}
        Int24 operator /( const Int24& val ) const {return Int24( getval() / val.getval() );}
        Int24 operator %( const Int24& val ) const {return Int24( getval() % val.getval() );}
        Int24 operator >>( const Int24& val ) const {return Int24( getval() >> val.getval() );}
        Int24 operator <<( const Int24& val ) const {return Int24( getval() << val.getval() );}

        Int24& operator +=( const Int24& val ) { setval(getval() + val.getval()); return *this; }
        Int24& operator -=( const Int24& val ) { setval(getval() - val.getval()); return *this; }
        Int24& operator *=( const Int24& val ) { setval(getval() * val.getval()); return *this; }
        Int24& operator /=( const Int24& val ) { setval(getval() / val.getval()); return *this; }
        Int24& operator %=( const Int24& val ) { setval(getval() % val.getval()); return *this; }
        Int24& operator <<=( const Int24& val ) { setval(getval() << val.getval()); return *this; }
        Int24& operator >>=( const Int24& val ) { setval(getval() >> val.getval()); return *this; }

        friend std::ostream& operator<<(std::ostream& out, const Int24& v) {out << v.getval(); return out;}
        friend std::istream& operator>>(std::istream& in, Int24& v) { int x; in>>x; v.setval(x); return in;}

};

#define decl_external_Int24_cmp_operator(NUM) \
bool operator ==( const Int24& a, NUM b) { return a.getval() == b; } \
bool operator !=( const Int24& a, NUM b) { return a.getval() != b; } \
bool operator >=( const Int24& a, NUM b) { return a.getval() >= b; } \
bool operator <=( const Int24& a, NUM b) { return a.getval() <= b; } \
bool operator > ( const Int24& a, NUM b) { return a.getval() >  b; } \
bool operator < ( const Int24& a, NUM b) { return a.getval() <  b; } \
bool operator ==( NUM b, const Int24& a) { return b == a.getval(); } \
bool operator !=( NUM b, const Int24& a) { return b != a.getval(); } \
bool operator >=( NUM b, const Int24& a) { return b >= a.getval(); } \
bool operator <=( NUM b, const Int24& a) { return b <= a.getval(); } \
bool operator > ( NUM b, const Int24& a) { return b >  a.getval(); } \
bool operator < ( NUM b, const Int24& a) { return b <  a.getval(); }

#define decl_external_Int24_arith_operator(NUM) \
NUM operator +( const Int24& val, NUM b ) {return val.getval() + b;} \
NUM operator -( const Int24& val, NUM b ) {return val.getval() - b;} \
NUM operator *( const Int24& val, NUM b ) {return val.getval() * b;} \
NUM operator /( const Int24& val, NUM b ) {return val.getval() / b;} \
NUM operator +( NUM b, const Int24& val ) {return b + val.getval();} \
NUM operator -( NUM b, const Int24& val ) {return b - val.getval();} \
NUM operator *( NUM b, const Int24& val ) {return b * val.getval();} \
NUM operator /( NUM b, const Int24& val ) {return b / val.getval();}

#define decl_external_Int24_intarith_operator(NUM) \
NUM operator %( const Int24& val, NUM b ) {return val.getval() % b;} \
Int24 operator >>( const Int24& val, NUM b ) {return Int24( val.getval() >> b );} \
Int24 operator <<( const Int24& val, NUM b ) {return Int24( val.getval() << b );} \
NUM operator %( NUM b, const Int24& val ) {return b % val.getval();} \
NUM operator >>( NUM b, const Int24& val ) {return b >> val.getval();} \
NUM operator <<( NUM b, const Int24& val ) {return b << val.getval();}

decl_external_Int24_cmp_operator(int);
decl_external_Int24_cmp_operator(long);
decl_external_Int24_cmp_operator(double);
decl_external_Int24_arith_operator(int);
decl_external_Int24_arith_operator(long);
decl_external_Int24_arith_operator(double);
decl_external_Int24_intarith_operator(int);
decl_external_Int24_intarith_operator(long);

newNTL_CLIENT;

typedef vector<double> vec_double;
//conversion string vers n'importe quoi
static istringstream* global_iss = 0;  
    template <typename T>
istringstream& operator>>(const std::string& s,T& obj)
{
    if (global_iss) delete global_iss;
    global_iss = new istringstream (s);
    *(global_iss) >> obj;
    return *(global_iss);
}

//parsing des arguments de main
#define PARSE_MAIN_ARGS ostringstream thesyntax; \
    for (int i=1;i<argc;i++)

#define MATCH_MAIN_ARGID(parsewhat,towhere) \
    if (i==1) thesyntax << "[ " << parsewhat << " " << towhere << " ] "; \
if (0==strcmp(argv[i],(parsewhat))) {string(argv[++i]) >> (towhere); continue;}

#define DETECT_MAIN_ARGID(parsewhat,towhere,value) \
    if (i==1) thesyntax << "[ " << parsewhat  << " ] "; \
if (0==strcmp(argv[i],(parsewhat))) {towhere = value; continue;}


#define SYNTAX(whatcorrect) \
    cerr << "Syntax error at argument" << i << endl; \
cerr << "Syntax: " << argv[0] << " " << thesyntax.str() << endl << whatcorrect << endl; \
exit(1);

//-------------------------------------------------------

//DEBUG functions
#define DECLARE_DEBUG_PRINT(T) \
    void PrintElement(T elt) {  cout << elt << endl; }

DECLARE_DEBUG_PRINT(ZZ);
DECLARE_DEBUG_PRINT(RR);
DECLARE_DEBUG_PRINT(vec_ZZ);
DECLARE_DEBUG_PRINT(vec_RR);
DECLARE_DEBUG_PRINT(mat_ZZ);
DECLARE_DEBUG_PRINT(mat_RR);

// functions which should be in newNTL

//vec % m
vec_ZZ operator%(const vec_ZZ& x, const ZZ& m) {
    int n = x.length();
    vec_ZZ r; r.SetLength(n);
    for (int i=1; i<=n; i++) {r(i)=x(i)%m;}
    return r;
}

//positive modulo (answer in [0,y-1])
//int posmod(RR x, ZZ y) {
//  return to_RR((x%y + y) % y);
//}
long posmod(ZZ x, long y) {
    return to_long((x%y + y) % y);
}

//square
template<typename T> T sqr(T x) {return x*x;}

//square norm
template<typename T>
inline T normsq(const vector<T>& v){
    return v*v;
}
inline RR norm(const vec_RR& v){
    return sqrt(v*v);
}
inline RR norm(const vec_ZZ& v){
    return sqrt(to_RR(v*v));
}

