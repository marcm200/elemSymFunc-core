
/*

    expressing any complete symmetric polynomial in terms
    of elementary symmetric functions in the variables

    to construct polynomials in two variables to describe
    hyperbolic components of exact period p for the
    quadratic Mandelbrot case z^12+c

    Marc Meidlinger
    June 2021

*/

#include "stdio.h"
#include "stdlib.h"
#include "stdint.h"
#include "math.h"
#include "time.h"
#include "string.h"


// const
const int32_t MINPERIOD=2;
const int32_t MAXSYLDIMENSION=1024;
const int32_t INT32MAX=0b01111111111111111111111111111111;
const int32_t SYMTERMSPERBLOCK=256;
const int32_t AZTERMSPERBLOCK=512;
const int32_t ANZSYMVARS=12;
const int32_t MAXEQ=1024;

const char elemsymname[]="ABCDEFGHIJKLMNOPQRSTUVWXYZ";


// globals
FILE *flog=NULL,*fdescr=NULL;
int64_t time0;
int64_t gmbiginthighestusedidx=0;
int32_t PERIOD=8; // period, number of roots in a cycle
int32_t BAREISS_STEP=2;
int32_t OUTPUTLEN=4;
int32_t REPEATITER=-1; // <= 0 => infinity, > 0 => only that number
int8_t AUTOMATICSIMPLIFY=1;
char fnbase[256];


// #defines
// no BRACKETS
#define NEWLOKAL \
    Lokal *var=new Lokal;\
    if (!var) return -1;

#define RETL(WW) \
{\
    delete var;\
    return (WW);\
}

#define WRITE32(FF,WW) \
{\
    int32_t tmp32=(WW);\
    if (fwrite(&tmp32,1,sizeof(tmp32),FF) != sizeof(tmp32) ) return 0;\
}

#define READ32(FF,ERG) \
{\
    int32_t tmp32;\
    if (fread(&tmp32,1,sizeof(tmp32),FF) != sizeof(tmp32) ) return 0;\
    (ERG)=tmp32;\
}

#define LOGMSG(TT) \
{\
    printf(TT);\
    if (flog) { fprintf(flog,TT); fflush(flog); } \
}

#define LOGMSG2(TT,AA) \
{\
    printf(TT,AA);\
    if (flog) { fprintf(flog,TT,AA); fflush(flog); } \
}

#define LOGMSG3(TT,AA,BB) \
{\
    printf(TT,AA,BB);\
    if (flog) { fprintf(flog,TT,AA,BB); fflush(flog); } \
}

#define LOGMSG4(TT,AA,BB,CC) \
{\
    printf(TT,AA,BB,CC);\
    if (flog) { fprintf(flog,TT,AA,BB,CC); fflush(flog); } \
}


// includes cpp
#include "dynslowstring.cpp"
#include "bigint.cpp"


// structs
struct Consts {
    BigInt bigm1,big1;
};

struct SymTerm {
    BigInt factor;
    int32_t cexponent;
    int32_t bexponent[ANZSYMVARS]; // b0..bMAX

    void setToZero(void);
    void setFactor(BigInt&);
    int32_t setFactor(const int64_t);
    int32_t setFactor(const char*);
    int32_t get_nzb(void);
    void set_bExponent_id_w(const int32_t,const int32_t);
    void set_cExponent(const int32_t);
    void sort_bExponents_dec(void);
    void clear_bExponents(void);
    void ausgabe(FILE*);
    void ausgabeB(FILE*);
    void copyFrom(SymTerm&);
    void shift_bExponents_left(void);

    int8_t save(FILE*);
    int8_t load(FILE*);

};

// stores expressions in elementary symmetric functions
struct AZterm {
    // small c and big C are DIFFERENT variables
    // small c is th seed, big C is the third elementryy
    BigInt factor;
    int32_t cexponent;
    int32_t azexponent[ANZSYMVARS]; // variables A,B,C,...Z

    void setToZero(void);
    void setFactor(BigInt&);
    int32_t setFactor(const int64_t);
    int32_t setFactor(const char*);
    void set_azExponent_id_w(const int32_t,const int32_t);
    void set_azExponent_ch_w(const char,const int32_t);
    void set_cExponent(const int32_t);
    void ausgabe(FILE*);
    void copyFrom(AZterm&);
    void invertSign(void);

    int8_t save(FILE*);
    int8_t load(FILE*);

};

struct SymPoly {
    int32_t anzterms;
    int32_t anzmem;
    SymTerm* terms;

    // using lexicographic ordering
    int32_t* sortidx;
    int8_t useLex;

    SymPoly();
    virtual ~SymPoly();

    void setToZero(void); // clears terms
    void setToOne(void);
    void setMemory(const int32_t);
    void ausgabe(FILE*,const int32_t);
    void copyFrom(SymPoly&);

    int32_t addTerm(SymTerm&);
    int32_t get_maxBexponent(void);
    int32_t subTerm(SymTerm&);
    int32_t addPoly(SymPoly&);
    int32_t subPoly(SymPoly&);
    int8_t provideMemory1more(void);
    int8_t searchExactTerm_A(SymTerm&);
    int8_t predZero(void);

    int8_t save(FILE*);
    int8_t load(FILE*);

    int32_t getLexicographicMax_TI(int32_t&,const int32_t);
    int32_t lexRemoveTermAtSort(const int32_t);

    int8_t findTermLexPosition_TA(int32_t&,SymTerm&);

};

struct AZPoly {
    int32_t anzterms;
    int32_t anzmem;
    AZterm* terms;

    // using lexicographic ordering
    int32_t* sortidx;
    int8_t useLex;

    AZPoly();
    virtual ~AZPoly();

    void setToZero(void); // clears terms
    void setToOne(void);
    void setMemory(const int32_t);
    void ausgabe(FILE*,const int32_t);
    void copyFrom(AZPoly&);

    char* getVariables(char*);

    int32_t addTerm(AZterm&);
    int32_t addTerm_FCVE(const int64_t,const int32_t,const int32_t,const int32_t);
    int32_t subTerm(AZterm&);
    int32_t addPoly(AZPoly&);
    int32_t subPoly(AZPoly&);
    int32_t get_hipow(const int32_t);
    int32_t coeff_TID(AZPoly&,const int32_t,const int32_t);
    int32_t getLexicographicMax_TI(int32_t&,const int32_t);
    int8_t provideMemory1more(void);
    int8_t predZero(void);
    int8_t nohigherthan(const int32_t);
    int8_t predSimple(void);

    int8_t save(FILE*);
    int8_t load(FILE*);
    int8_t findTermLexPosition_TA(int32_t&,AZterm&);
    int32_t lexRemoveTermAtSort(const int32_t);

    AZPoly* get_coeffList_TV(int32_t&,const int32_t);

    void invertSign(void);

};

struct MatrixPolynom {
	int32_t dim;

	AZPoly *entryYX; // accessed via Y-index * dim + X-index

	MatrixPolynom();
	virtual ~MatrixPolynom();

	void setDimension(const int32_t);
	void ausgabe(FILE*);
	void setConstant0polynom(void);
	int8_t entryZero_YX(const int32_t,const int32_t);
	AZPoly* getPointer1based_YX(const int32_t,const int32_t);
	int8_t save(const char*);
	int8_t load(const char*,const int32_t,const int32_t);

};

typedef MatrixPolynom *PMatrixPolynom;

struct DivHlp {
    BigInt irem;
    AZPoly hlprest,hlpt2,B2;
};

struct DivHlpSym {
    BigInt irem;
    SymPoly hlprest,hlpt2;
};

struct AZRational {
    // simple rational type
    AZPoly zaehler,nenner;

    void ausgabe(FILE*,const int32_t);
    void copyFrom(AZRational&);
    void setToZero(void);

    int8_t save(FILE*);
    int8_t load(FILE*);

};

struct Solution {
    int8_t solved;
    int32_t varidx; // elemsymname
    AZRational eq;

    Solution();
    void setToUnsolved(void);
    void ausgabe(FILE*,const int32_t);
    void setSolution_id(const int32_t,AZRational&);

    int8_t save(FILE*);
    int8_t load(FILE*);

};

struct Equation0 {
    int32_t id; // globally defining id, must be distinct to all other
    // every computed equations
    AZPoly azleft;
    SymPoly nonright;
    // equation of the azleft = nonright
    // the goal of the manipulations is the get a non-part of 0

    Equation0();
    void setAll_ABC(const int32_t,AZPoly&,SymPoly&);
    void ausgabe(FILE*,const int32_t);
    void copyFrom(Equation0&);

    char* getVariablesLeft(char*);

    int32_t shiftPureCtermsToLeft(void);
    int32_t reduce_gcd_c(void);
    int32_t expressElemRight(void);

    int32_t fullSimplify(void);

    int32_t solveForLinear_TI(AZRational&,const int32_t);

    int8_t save(FILE*);
    int8_t load(FILE*);

};

struct EquationList {
    Equation0 aseq[MAXEQ];
    int32_t anzeq;

    EquationList();
    void fastClear(void);

    void ausgabe(FILE*,const int32_t);
    void ausgabeVars(FILE*,const int32_t);
    int32_t addEquation(const int32_t,AZPoly&,SymPoly&);
    int32_t addEquation(Equation0&);
    int32_t removeEquationById(const int32_t);
    int32_t get_VariablesLeft(char*);
    int32_t copyFrom(EquationList&);
    int32_t makeAddEquation_fullSymmetric_IFCBs(const int32_t,const int64_t,const int32_t,const char*);
    Equation0* getEquationPtr_id(const int32_t);

    int8_t save(FILE*);
    int8_t load(FILE*);

};


// globals 2
Consts *consts=NULL;
SymPoly sigma[ANZSYMVARS+1]; // explicit elementary symmetric functions


// forward
int8_t getFile_TF(DynSlowString&,const char*);
void initSigma(const int32_t);

int32_t eqmultiply_TIAFC(Equation0&,const int32_t,Equation0&,const int32_t,const int32_t);
int32_t eqadd_TIAB(Equation0&,const int32_t,Equation0&,Equation0&);

int8_t symterm_bcExponentVgl_AB(SymTerm&,SymTerm&);
int8_t azterm_azcExponentVgl_AB(AZterm&,AZterm&);
int8_t symtermAllVgl_AB(SymTerm&,SymTerm&);
int8_t termVgl_only_bExponent_AB(SymTerm&,SymTerm&);
int32_t substitute_TAB(AZRational&,AZPoly&,Solution&);
int32_t eliminate_TTTIEAB(AZPoly&,AZPoly&,AZPoly&,const int32_t,const int32_t,AZPoly&,AZPoly&);

int32_t polyLtbexponent_TA(SymTerm&,SymPoly&);
int32_t termMul_TAB(SymTerm&,SymTerm&,SymTerm&);
int32_t polyAdd_TAB(SymPoly&,SymPoly&,SymPoly&);
int32_t polySub_TAB(SymPoly&,SymPoly&,SymPoly&);
int32_t polyMul_TAB(SymPoly&,SymPoly&,SymPoly&);
int32_t polyPow_TAE(SymPoly&,SymPoly&,const int32_t);
int32_t polyMul_TTermA(SymPoly&,SymTerm&,SymPoly&);
int32_t polyDiv_TRABLH(SymPoly&,SymPoly&,SymPoly&,SymPoly&,const int32_t,DivHlpSym&);

int32_t fractionFreeGaussMultistep2_TABI(AZPoly&,AZPoly&,AZPoly&,const int32_t);

int32_t aztermMul_TAB(AZterm&,AZterm&,AZterm&);
int32_t azpolyAdd_TAB(AZPoly&,AZPoly&,AZPoly&);
int32_t azpolySub_TAB(AZPoly&,AZPoly&,AZPoly&);
int32_t azpolyMul_TAB(AZPoly&,AZPoly&,AZPoly&);
int32_t azpolyMul_TTermA(AZPoly&,AZterm&,AZPoly&);
int32_t convertSymToAZ_TA(AZPoly&,SymPoly&);
int32_t azpolyPow_TNA(AZPoly&,const int32_t,AZPoly&);
int32_t azpolyDiv_TRABLH(AZPoly&,AZPoly&,AZPoly&,AZPoly&,const int32_t,DivHlp&);
int32_t azpolySqr_TA(AZPoly&,AZPoly&);

int32_t createFullySymmetric_TA(SymPoly&,SymTerm&);
int32_t replace_iter_TA(SymPoly&,SymPoly&);
int32_t replace_repeated_iter_TA(SymPoly&,SymPoly&);
int32_t replace_repeated_iter(SymPoly&);
int32_t replace_iter_TTerm(SymPoly&,SymTerm&);
int32_t replace_iter_T_bid_expo(SymPoly&,const int32_t,const int32_t);
int32_t cyclicSum_TA(SymPoly&,SymTerm&);
int64_t fakultaet(const int64_t);
int64_t getbPermutationCount(SymTerm&);
int32_t split_TsymTnonA(SymPoly&,SymPoly&,SymPoly&);
int32_t expressByElem_Texp_Tnon_A(AZPoly&,SymPoly&,SymPoly&);
int32_t makeEquationSymPoly_TAI(Equation0&,SymPoly&,const int32_t);
int32_t convertStrToBTerm_TA(SymTerm&,const char*);
int32_t constructBetaBin_TI_prodA0A1B_FsumT(Equation0&,const int32_t,const int32_t,const int32_t,const int32_t,const int32_t,const char*);

int8_t predPermutation_AB(int32_t*,int32_t*);
inline int32_t maximumI(const int32_t,const int32_t);

int8_t setSylvesterMatrix_TABCDVS(MatrixPolynom&,AZPoly&,AZPoly&,AZPoly*,AZPoly*,const int32_t,const int32_t);


// routines
int32_t maximumI(const int32_t a,const int32_t b) {
    if (a > b) return a;

    return b;
}

int32_t getElemSymNameIdx(const char ch) {
    char up=ch;
    if ( (up >= 'a') && (up <= 'z') ) up=ch-'a'+'A';

    for(int32_t k=strlen(elemsymname);k>=0;k--) {
        if (elemsymname[k] == up) return k;
    }

    return -1; // error

}

char* chomp(char* s) {
	if (!s) return 0;
	for(int i=strlen(s);i>=0;i--) if (s[i]<32) s[i]=0; else break;
	return s;
}

char* upper(char* s) {
	if (!s) return 0;
	for(unsigned int i=0;i<strlen(s);i++) {
		if ((s[i]>='a')&&(s[i]<='z')) s[i]=s[i]-'a'+'A';
		if (s[i]=='ö') s[i]='Ö';
		if (s[i]=='ü') s[i]='Ü';
		if (s[i]=='ä') s[i]='Ä';
	}

	return s;
}

void initConstants(void) {
    consts=new Consts;
    if (!consts) {
        fprintf(stderr,"\nerror. memory. consts.\N");
        exit(99);
    }

    consts->bigm1.setstr("-1");
    consts->big1.setstr("1");

}

// struct SymTerm
void SymTerm::clear_bExponents(void) {
    for(int a=0;a<PERIOD;a++) bexponent[a]=0;
}

void SymTerm::shift_bExponents_left(void) {
    int32_t b0=bexponent[0];
    for(int32_t k=1;k<PERIOD;k++) {
        bexponent[k-1]=bexponent[k];
    }
    bexponent[PERIOD-1]=b0;
}

int32_t SymTerm::get_nzb(void) {
    // counts the number of non-zero b-exponents
    int32_t erg=0;
    for(int32_t k=0;k<PERIOD;k++) {
        if (bexponent[k] != 0) erg++;
    }

    return erg;
}

void SymTerm::copyFrom(SymTerm& A) {
    factor.copyFrom(A.factor);
    cexponent=A.cexponent;
    for(int32_t k=0;k<PERIOD;k++) {
        bexponent[k]=A.bexponent[k];
    } // k

}

int32_t SymTerm::setFactor(const char* A) {
    // ret 0 => success, !=0 = >error
    if (factor.setstr(A) != 0) return -1;

    return 0;
}

int32_t SymTerm::setFactor(const int64_t A) {
    // ret 0 => success, !=0 = >error
    char tt[32];
    sprintf(tt,"%I64d",A);

    if (factor.setstr(tt) != 0) return -1;

    return 0;
}

void SymTerm::setFactor(BigInt& A) {
    factor.copyFrom(A);
}

void SymTerm::setToZero(void) {
    factor.setToZero();
    for(int32_t i=0;i<PERIOD;i++) bexponent[i]=0;
    cexponent=0;
}

void SymTerm::set_cExponent(const int32_t aw) {
    if (aw < 0) {
        printf("\nerror. negative c-exponents not supported.\n");
        exit(99);
    }
    cexponent=aw;
}

void SymTerm::set_bExponent_id_w(
    const int32_t pos,
    const int32_t ex
) {
    if ( (pos < 0) || (pos >= PERIOD) ) {
        printf("\nerror. set exponent out of var range: symterm%i\n",pos);
        exit(99);
    }

    if (ex < 0) {
        printf("\nerror. negative b-exponents not supported.\n");
        exit(99);
    }

    bexponent[pos]=ex;

}

void SymTerm::sort_bExponents_dec(void) {
    // sort exponents in decreasing order
    int8_t changed=1;
    while (changed > 0) {
        changed=0;

        for(int32_t a=0;a<PERIOD;a++) {
            for(int32_t b=(a+1);b<PERIOD;b++) {
                if (bexponent[a] < bexponent[b]) {
                    changed=1;
                    int32_t ex=bexponent[a];
                    bexponent[a]=bexponent[b];
                    bexponent[b]=ex;
                }
            } // b
        } // a

    } // while

}

void SymTerm::ausgabe(FILE* f) {
    factor.ausgabe(f);
    if (cexponent == 1) fprintf(f,"*c");
    else if (cexponent > 1) fprintf(f,"*c^%i",cexponent);

    for(int32_t k=0;k<PERIOD;k++) {
        if (bexponent[k] > 1) {
            fprintf(f,"*b%i^%i",k,bexponent[k]);
        } else
        if (bexponent[k] == 1) {
            fprintf(f,"*b%i",k);
        }
    } // k

}

void SymTerm::ausgabeB(FILE* f) {

    for(int32_t k=0;k<PERIOD;k++) {
        if (bexponent[k] > 1) {
            fprintf(f,"*b%i^%i",k,bexponent[k]);
        } else
        if (bexponent[k] == 1) {
            fprintf(f,"*b%i",k);
        }
    } // k

}

// AZRational
int8_t AZRational::save(FILE* f){
    if (!f) return 0;

    if (zaehler.save(f) <= 0) return 0;
    if (nenner.save(f) <= 0) return 0;

    return 1;

}

int8_t AZRational::load(FILE* f){
    if (!f) return 0;

    if (zaehler.load(f) <= 0) return 0;
    if (nenner.load(f) <= 0) return 0;

    return 1;

}

void AZRational::setToZero(void) {
    zaehler.setToZero();
    zaehler.setToOne();
}

void AZRational::copyFrom( AZRational& A ) {
    zaehler.copyFrom( A.zaehler );
    nenner.copyFrom( A.nenner );
}

void AZRational::ausgabe(FILE* f,const int32_t maxt) {
    fprintf(f,"( ");
    zaehler.ausgabe(f,maxt);
    fprintf(f,") / (");
    nenner.ausgabe(f,maxt);
    fprintf(f,")");
}

// AZpoly
AZPoly* AZPoly::get_coeffList_TV(int32_t& hip,const int32_t VARIDX) {
    // returns a newlay allocated array of coeffiicent
    // polynomials for variable VARIDX (A--Z, not seed c)

    if ( (VARIDX < 0) || (VARIDX >= ANZSYMVARS) ) return NULL;

    hip=get_hipow(VARIDX);
    AZPoly* co=new AZPoly[ hip+1 ];

    for(int32_t k=0;k<=hip;k++) {
        if (coeff_TID(co[k],VARIDX,k) != 0) {
            delete[] co;
            return NULL;
        }
    }

    return co;

}

int8_t AZPoly::save(FILE* f) {
    if (!f) return 0;

    WRITE32(f,anzterms)
    for(int32_t k=0;k<anzterms;k++) {
        if (terms[k].save(f) <= 0) return 0;
    } // k

    return 1;

}

int8_t AZPoly::load(FILE* f) {
    if (!f) return 0;

    setToZero();
    useLex=0; // stabndard: off

    int32_t w32;
    READ32(f,w32)
    AZterm azt;

    for(int32_t k=0;k<w32;k++) {
        if (azt.load(f) <= 0) return 0;
        if (addTerm( azt ) != 0) return 0;
    } // k

    return 1;

}

int32_t AZPoly::getLexicographicMax_TI(
    int32_t& erg,
    const int32_t lvaridx // this is the lexicaographically MOST important variable
) {
    if (
        (lvaridx < 0) ||
        (lvaridx >= ANZSYMVARS)
    ) {
        printf("\nerror getLex variable out of rnage\n");
        exit(99);
    }

    // first c-exponent then az
    if (anzterms <= 0) return -1; // error

    erg=0; // first term
    for(int32_t k=1;k<anzterms;k++) {
        if (terms[erg].azexponent[lvaridx] < terms[k].azexponent[lvaridx]) {
            erg=k;
            continue;
        }

        if (terms[erg].azexponent[lvaridx] > terms[k].azexponent[lvaridx]) {
            continue;
        }

        if (azterm_azcExponentVgl_AB( terms[erg], terms[k] ) < 0) erg=k;
    } // k

    return 0;

}

char* AZPoly::getVariables(char* erg) {
    // erg must provied enpough mnemory for AZSYMVARS+1
    for(int32_t k=0;k<(1+ANZSYMVARS);k++) {
        erg[k]='-';
    }
    erg[ANZSYMVARS + 1]=0;

    for(int32_t k=0;k<anzterms;k++) {
        if (terms[k].cexponent != 0) erg[ANZSYMVARS]='c';

        for(int32_t v=0;v<ANZSYMVARS;v++) {
            if (terms[k].azexponent[v] != 0) erg[v]=elemsymname[v];
        }
    }

    return erg;

}

int8_t AZPoly::predSimple(void) {
    // ret <= 0 => no, > 0 => yes
    // all coefficients are numbers and or c-terms (seed)
    // or A but no other variables
    for(int32_t k=0;k<anzterms;k++) {
        for(int32_t v=1;v<ANZSYMVARS;v++) {
            // start at 1 as A is allowed
            if (terms[k].azexponent[v] != 0) return 0;
        } // v
    } // k

    return 1;

}

int8_t AZPoly::nohigherthan(const int32_t IDX) {
    if ( (IDX < 0) || (IDX >= ANZSYMVARS) ) {
        printf("\nerror. nohigher for unknown variable IDX=%i\n",IDX);
        exit(99);
    }

    // ret: bool
    for(int32_t k=0;k<anzterms;k++) {

        for(int32_t v=(IDX+1);v<PERIOD;v++) {
            if (terms[k].azexponent[v] != 0) return 0; // higher present
        } // v

    } // k

    return 1;
}

int32_t AZPoly::coeff_TID(
    AZPoly& erg,
    const int32_t VAR,
    const int32_t DEG
) {
    erg.setToZero();

    // ret < 0 => error
    // ret = 0 => success
    // currently NOT working for c-coefficients
    if ( (VAR < 0) || (VAR >= ANZSYMVARS) ) return -1; // error

    // all terms that have azexponent[VAR]=DEG as a polynomial
    // but that exponent set to 0
    int32_t error=0;
    for(int32_t k=0;k<anzterms;k++) {
        if (terms[k].azexponent[VAR] != DEG) continue;

        AZterm azt;
        azt.copyFrom( terms[k] );
        azt.azexponent[VAR]=0; // remove tha texponent
        error += erg.addTerm(azt);
    } // k

    return error;

}

int32_t AZPoly::get_hipow(const int32_t varidx) {
    if ( (varidx < 0) || (varidx >= ANZSYMVARS) ) {
        printf("\nerror. var symbol %i not valid.\n",varidx);
        exit(99);
    }

    // varidx=0 => A, B=1 etc
    int32_t hi=0;
    for(int32_t k=0;k<anzterms;k++) {
        if (terms[k].azexponent[varidx] > hi) {
            hi=terms[k].azexponent[varidx];
        }
    } // k

    return hi;
}

int8_t AZPoly::predZero(void) {
    if (anzterms > 1) return 0;
    if (anzterms <= 0) return 1;

    if (terms[0].factor.vorz != 0) return 0;

    return 1;
}

void AZPoly::invertSign(void) {
    for(int32_t k=0;k<anzterms;k++) {
        terms[k].invertSign();
    } // k
}

AZPoly::AZPoly() {
    anzterms=0;
    anzmem=0;
    terms=NULL;
    useLex=1;
    sortidx=NULL;
};

AZPoly::~AZPoly() {
    if (terms) delete[] terms;
    if (useLex > 0) if (sortidx) delete[] sortidx;
}

void AZPoly::setToZero(void) {
    anzterms=0;
}

void AZPoly::setToOne(void) {
    anzterms=0;
    AZterm T;
    T.setToZero();
    T.setFactor(1);

    if (addTerm(T) != 0) {
        printf("\nerror. az set to one error.\n");
        exit(99);
    }
}

void AZPoly::ausgabe(FILE* f,const int32_t maxt) {
    if (anzterms == 0) {
        fprintf(f,"0");
        return;
    }

    int32_t te=anzterms;
    if (maxt > 0) if (te > maxt) te=maxt;

    for(int32_t k=0;k<te;k++) {
        terms[k].ausgabe(f);
        if (k != (anzterms-1)) fprintf(f,"+");
    } // k

    if (te != anzterms) {
        fprintf(f," ... (%i more) ",anzterms-te);
    }

}

int8_t AZPoly::provideMemory1more(void) {
    // ret <= 0 => error
    // > 0 => success
    // enough memory for one more term ?

    if (
        (!terms) ||
        ( (anzterms+1) >= anzmem )
    ) {
        // not enough memory
        int32_t newmem=anzmem + AZTERMSPERBLOCK;
        if (newmem < anzmem) {
            // overflow
            return 0;
        }

        AZterm* t2=new AZterm[newmem];
        if (!t2) {
            return 0;
        }

        // copy if present
        if ( (terms) && (anzterms > 0) )  {
            for(int32_t k=0;k<anzterms;k++) {
                t2[k].copyFrom(terms[k]);
            }
            delete[] terms;
        }

        terms=t2;

        if (useLex > 0) {
            int32_t* t2=new int32_t[newmem];
            if (!t2) {
                return 0;
            }

            // copy if present
            if ( (sortidx) && (anzterms > 0) )  {
                for(int32_t k=0;k<anzterms;k++) {
                    t2[k]=sortidx[k];
                }
                delete[] sortidx;
            }

            sortidx=t2;

        } // useLex

        // anzterms remains unchangted
        anzmem=newmem;

    } // new memory

    return 1;

}

void AZPoly::setMemory(const int32_t a) {
    anzterms=0;
    if (a < anzmem) return; // enough memory

    delete[] terms;
    terms=new AZterm[a];
    if (!terms) {
        fprintf(stderr,"\nerror. memory. az poly set\n");
        exit(99);
    }
    anzmem=a;

    if (useLex > 0) {
        delete[] sortidx;
        sortidx=new int32_t[a];
        if (!sortidx) {
            fprintf(stderr,"\nerror. memory2. az poly set\n");
            exit(99);
        }
    }

}

int32_t AZPoly::subTerm(AZterm& A) {
    // ret 0 => success, !=0 = >error
    AZterm B;
    B.copyFrom(A);
    B.factor.vorz *= -1;

    return addTerm(B);
}

int32_t AZPoly::addTerm_FCVE(
    const int64_t afactor,
    const int32_t acexp,
    const int32_t avar,
    const int32_t aexp
) {
    if ( (avar < 0) || (avar >= ANZSYMVARS) ) {
        printf("\nerror. variable out of A-Z\n");
        exit(99);
    }

    AZterm T;
    T.setToZero();
    T.setFactor( afactor );
    T.set_cExponent( acexp );
    T.set_azExponent_id_w( avar,aexp );

    return addTerm(T);
}

int32_t AZPoly::addPoly(AZPoly& A) {
    // ret 0 => success, !=0 = >error
    // zero terms are handled in addTerm
    for(int32_t a=0;a<A.anzterms;a++) {
        if (addTerm(A.terms[a]) != 0) return -1;
    } // a

    return 0;

}

int32_t AZPoly::subPoly(AZPoly& A) {
    // ret <= 0 => error
    // > 0 => success
    // subTemr handles zro terms

    for(int32_t a=0;a<A.anzterms;a++) {
        if (subTerm(A.terms[a]) != 0) return -1;
    } // a

    return 0;

}

int8_t AZPoly::findTermLexPosition_TA(
    int32_t& pos,
    AZterm& A
) {
    if (anzterms <= 0) {
        pos=0;
        return 0; // not present
    }

    // ret <= 0 => not present, pos is the position of the
    // term A, so that all before that are SMALLER according to
    // azterm_azcExponentVgl_AB

    // ret > 0 => pos is the position of the CAZ-exponents of A

    // binary seearch
    if (azterm_azcExponentVgl_AB(A,terms[sortidx[0]]) < 0) {
        // left outside
        pos=0;
        return 0;
    }

    // right outside ?
    if (azterm_azcExponentVgl_AB(A,terms[sortidx[anzterms-1]]) > 0) {
        pos=anzterms; // outside new
        return 0;
    }

    // now may be inside
    int32_t l=0,r=anzterms-1;
    while (1) {
        int32_t diff=r-l+1;
        if (diff > 4) {
            // binary search at midpoint
            int32_t m=l+(diff >> 1);
            int8_t vgl=azterm_azcExponentVgl_AB(A,terms[sortidx[m]]);
            if (vgl == 0) {
                // found
                pos=m;
                return 1;
            }

            if (vgl < 0) {
                // A < m => in left half
                r=m;
                continue;
            }

            if (vgl > 0) {
                // A > m => right half
                l=m;
                continue;
            }

        } else {
            // a few elements, look linearly
            for(int32_t k=l;k<=r;k++) {
                int8_t vgl=azterm_azcExponentVgl_AB(A,terms[sortidx[k]]);
                if (vgl == 0) {
                    pos=k;
                    return 1; // present
                }

                if (vgl < 0) {
                    // A < term[ortidx[k]] for first time
                    pos=k;
                    return 0; // not present
                    // new term should be put at k and everyhting moved right
                }
            } // k

            // error
            printf("\nerror binary search.\n");
            exit(99);

        } //

    } // while

    // error

    printf("\nerror binary search end of func.\n");
    exit(99);

}

int32_t AZPoly::lexRemoveTermAtSort(const int32_t insort) {
    // ret < 0 => error, 0 => success

    int32_t interm=sortidx[insort];

    if (interm == (anzterms-1)) {
        // move sortidx one to the left
        for(int32_t k=insort;k<(anzterms-1);k++) {
            sortidx[k]=sortidx[k+1];
        } // k

        // in terms only the lastg element must be removed
        // this is done by just decreasing anzterms
        anzterms--;

        return 0;

    }

    // the term in terms is inside
    // so terms must be shuffled too
    // swap the last terms[anzterms-1] onto terms[interm]
    // and adjust sortidx accordingly
    for(int32_t k=insort;k<(anzterms-1);k++) {
        sortidx[k]=sortidx[k+1];
    } // k

    terms[interm].copyFrom(terms[anzterms-1]);

    for(int32_t k=0;k<(anzterms-1);k++) {
        if (sortidx[k] == (anzterms-1)) {
            sortidx[k]=interm;
            break;
        }
    } // k

    anzterms--;

    return 0; // success

}

int32_t AZPoly::addTerm(AZterm& A) {

    if (A.factor.isZero() > 0) return 0; // nothing to add, but no error

    // ret 0 => success, !=0 = >error
    if (provideMemory1more() <= 0) return -1;

    int32_t error=0;

    if (useLex > 0) {
        int32_t insort=-1;
        int8_t present=findTermLexPosition_TA(insort,A);

        // present > 0 => add A to term[pos]
        // and remove zero if necessary
        if (present > 0) {
            BigInt w;
            if (bigintAdd_TAB(w,terms[sortidx[insort]].factor,A.factor) != 0) {
                printf("\n  error. bigint. probably overflow.");
                exit(99);
            }

            if (w.isZero() > 0) {
                error += lexRemoveTermAtSort(insort);
            } else {
                terms[sortidx[insort]].factor.copyFrom( w );
            }

            return error;

        } // present
        else {
            // present <= 0
            // new term: set at position pos
            // (pos=0..anzterm. if anzterm => add at end
            // and move pos+1 on to the right
            // memory enough to store an additional element

            for(int32_t k=anzterms;k>insort;k--) {
                sortidx[k]=sortidx[k-1];
            } // k
            terms[anzterms].copyFrom( A );
            sortidx[insort]=anzterms;
            anzterms++;

        } // not present

        return 0; // succes

    } // useLex

    // does the term itself already exist ?
    int32_t idx=-1;
    for(int32_t k=0;k<anzterms;k++) {
        if (azterm_azcExponentVgl_AB(terms[k],A) == 0) {
            idx=k;
            break;
        }
    } // k

    if (idx < 0) {
        // add new term, memory sufficient
        // term not zero
        terms[anzterms].copyFrom(A);
        anzterms++;
    } else {
        // add the factors
        BigInt w;
        if (bigintAdd_TAB(w,terms[idx].factor,A.factor) != 0) {
            printf("\n  error. bigint. probably overflow.");
            exit(99);
        }

        if (w.isZero() > 0) {
            // term be removed
            if (idx != (anzterms-1) ) {
                terms[idx].copyFrom(terms[anzterms-1]);
            }

            anzterms--;

        } else {
            terms[idx].factor.copyFrom(w);
        }

    }

    return 0;

}

void AZPoly::copyFrom(AZPoly& A) {
    setToZero();
    setMemory(A.anzterms+AZTERMSPERBLOCK);

    if ( (useLex > 0) && (A.useLex > 0) ) {
        // copy
        for(int32_t k=0;k<A.anzterms;k++) {
            terms[k].copyFrom( A.terms[k] );
            sortidx[k]=A.sortidx[k];
        } // k

        anzterms=A.anzterms;

        return;
    }

    if ( (useLex > 0) && (A.useLex <= 0) ) {
        setToZero();
        for(int32_t k=0;k<A.anzterms;k++) {
            // use addTerm to generate sortidx
            if (addTerm(A.terms[k] ) != 0) {
                printf("\nerror. copy azpoly\n");
                exit(99);
            }
        } // k

        return;

    }

    for(int32_t k=0;k<A.anzterms;k++) {
        terms[k].copyFrom( A.terms[k] );
    } // k

    anzterms=A.anzterms;

}

// SymPoly lex
int8_t SymPoly::predZero(void) {
    if (anzterms > 1) return 0;
    if (anzterms <= 0) return 1;

    if (terms[0].factor.vorz != 0) return 0;

    return 1;
}

int8_t SymPoly::save(FILE* f) {
    if (!f) return 0;

    WRITE32(f,anzterms)
    for(int32_t k=0;k<anzterms;k++) {
        if (terms[k].save(f) <= 0) return 0;
    } // k

    return 1;

}

int8_t SymPoly::load(FILE* f) {
    if (!f) return 0;

    setToZero();

    int32_t w32;
    READ32(f,w32)
    SymTerm azt;

    for(int32_t k=0;k<w32;k++) {
        if (azt.load(f) <= 0) return 0;
        if (addTerm( azt ) != 0) return 0;
    } // k

    return 1;

}

int32_t SymPoly::getLexicographicMax_TI(
    int32_t& erg,
    const int32_t lvaridx // this is the lexicaographically MOST important variable
) {
    if (
        (lvaridx < 0) ||
        (lvaridx >= ANZSYMVARS)
    ) {
        printf("\nerror getLex variable out of rnage\n");
        exit(99);
    }

    // first c-exponent then az
    if (anzterms <= 0) return -1; // error

    erg=0; // first term
    for(int32_t k=1;k<anzterms;k++) {
        if (terms[erg].bexponent[lvaridx] < terms[k].bexponent[lvaridx]) {
            erg=k;
            continue;
        }

        if (terms[erg].bexponent[lvaridx] > terms[k].bexponent[lvaridx]) {
            continue;
        }

        if (symterm_bcExponentVgl_AB( terms[erg], terms[k] ) < 0) erg=k;
    } // k

    return 0;

}

SymPoly::SymPoly() {
    anzterms=0;
    anzmem=0;
    terms=NULL;
    useLex=1; // use sorted list for terms w.r.t c and b-exponents
    sortidx=NULL;
};

SymPoly::~SymPoly() {
    if (terms) delete[] terms;
    if (useLex > 0) if (sortidx) delete[] sortidx;
}

void SymPoly::setToZero(void) {
    anzterms=0;
}

void SymPoly::setToOne(void) {
    anzterms=0;
    SymTerm T;
    T.setToZero();
    T.setFactor(1);

    if (addTerm(T) != 0) {
        printf("\nerror. az set to one error.\n");
        exit(99);
    }
}
void SymPoly::ausgabe(FILE* f,const int32_t maxt) {
    if (anzterms == 0) {
        fprintf(f,"0");
        return;
    }

    int32_t te=anzterms;
    if (maxt > 0) if (te > maxt) te=maxt;

    for(int32_t k=0;k<te;k++) {
        terms[k].ausgabe(f);
        if (k != (anzterms-1)) fprintf(f,"+");
    } // k

    if (te != anzterms) {
        fprintf(f," ... (%i more) ",anzterms-te);
    }

}

int8_t SymPoly::provideMemory1more(void) {
    // ret <= 0 => error
    // > 0 => success
    // enough memory for one more term ?

    if (
        (!terms) ||
        ( (anzterms+1) >= anzmem )
    ) {
        // not enough memory
        int32_t newmem=anzmem + SYMTERMSPERBLOCK;
        if (newmem < anzmem) {
            // overflow
            return 0;
        }

        SymTerm* t2=new SymTerm[newmem];
        if (!t2) {
            return 0;
        }

        // copy if present
        if ( (terms) && (anzterms > 0) )  {
            for(int32_t k=0;k<anzterms;k++) {
                t2[k].copyFrom(terms[k]);
            }
            delete[] terms;
        }

        terms=t2;

        if (useLex > 0) {
            int32_t* t2=new int32_t[newmem];
            if (!t2) {
                return 0;
            }

            // copy if present
            if ( (sortidx) && (anzterms > 0) )  {
                for(int32_t k=0;k<anzterms;k++) {
                    t2[k]=sortidx[k];
                }
                delete[] sortidx;
            }

            sortidx=t2;

        } // useLex

        // anzterms remains unchangted
        anzmem=newmem;

    } // new memory

    return 1;

}

void SymPoly::setMemory(const int32_t a) {
    anzterms=0;
    if (a < anzmem) return; // enough memory

    delete[] terms;
    terms=new SymTerm[a];
    if (!terms) {
        fprintf(stderr,"\nerror. memory. az poly set\n");
        exit(99);
    }
    anzmem=a;

    if (useLex > 0) {
        delete[] sortidx;
        sortidx=new int32_t[a];
        if (!sortidx) {
            fprintf(stderr,"\nerror. memory2. az poly set\n");
            exit(99);
        }
    }

}

int32_t SymPoly::subTerm(SymTerm& A) {
    // ret 0 => success, !=0 = >error
    SymTerm B;
    B.copyFrom(A);
    B.factor.vorz *= -1;

    return addTerm(B);
}

int32_t SymPoly::addPoly(SymPoly& A) {
    // ret 0 => success, !=0 = >error
    // zero terms are handled in addTerm
    for(int32_t a=0;a<A.anzterms;a++) {
        if (addTerm(A.terms[a]) != 0) return -1;
    } // a

    return 0;

}

int32_t SymPoly::subPoly(SymPoly& A) {
    // ret <= 0 => error
    // > 0 => success
    // subTemr handles zro terms

    for(int32_t a=0;a<A.anzterms;a++) {
        if (subTerm(A.terms[a]) != 0) return -1;
    } // a

    return 0;

}

int8_t SymPoly::findTermLexPosition_TA(
    int32_t& pos,
    SymTerm& A
) {
    if (anzterms <= 0) {
        pos=0;
        return 0; // not present
    }

    // ret <= 0 => not present, pos is the position of the
    // term A, so that all before that are SMALLER according to
    // azterm_azcExponentVgl_AB

    // ret > 0 => pos is the position of the CAZ-exponents of A

    // binary seearch
    if (symterm_bcExponentVgl_AB(A,terms[sortidx[0]]) < 0) {
        // left outside
        pos=0;
        return 0;
    }

    // right outside ?
    if (symterm_bcExponentVgl_AB(A,terms[sortidx[anzterms-1]]) > 0) {
        pos=anzterms; // outside new
        return 0;
    }

    // now may be inside
    int32_t l=0,r=anzterms-1;
    while (1) {
        int32_t diff=r-l+1;
        if (diff > 4) {
            // binary search at midpoint
            int32_t m=l+(diff >> 1);
            int8_t vgl=symterm_bcExponentVgl_AB(A,terms[sortidx[m]]);
            if (vgl == 0) {
                // found
                pos=m;
                return 1;
            }

            if (vgl < 0) {
                // A < m => in left half
                r=m;
                continue;
            }

            if (vgl > 0) {
                // A > m => right half
                l=m;
                continue;
            }

        } else {
            // a few elements, look linearly
            for(int32_t k=l;k<=r;k++) {
                int8_t vgl=symterm_bcExponentVgl_AB(A,terms[sortidx[k]]);
                if (vgl == 0) {
                    pos=k;
                    return 1; // present
                }

                if (vgl < 0) {
                    // A < term[ortidx[k]] for first time
                    pos=k;
                    return 0; // not present
                    // new term should be put at k and everyhting moved right
                }
            } // k

            // error
            printf("\nerror binary search.\n");
            exit(99);

        } //

    } // while

    // error

    printf("\nerror binary search end of func.\n");
    exit(99);

}

int32_t SymPoly::lexRemoveTermAtSort(const int32_t insort) {
    // ret < 0 => error, 0 => success

    int32_t interm=sortidx[insort];

    if (interm == (anzterms-1)) {
        // move sortidx one to the left
        for(int32_t k=insort;k<(anzterms-1);k++) {
            sortidx[k]=sortidx[k+1];
        } // k

        // in terms only the lastg element must be removed
        // this is done by just decreasing anzterms
        anzterms--;

        return 0;

    }

    // the term in terms is inside
    // so terms must be shuffled too
    // swap the last terms[anzterms-1] onto terms[interm]
    // and adjust sortidx accordingly
    for(int32_t k=insort;k<(anzterms-1);k++) {
        sortidx[k]=sortidx[k+1];
    } // k

    terms[interm].copyFrom(terms[anzterms-1]);

    for(int32_t k=0;k<(anzterms-1);k++) {
        if (sortidx[k] == (anzterms-1)) {
            sortidx[k]=interm;
            break;
        }
    } // k

    anzterms--;

    return 0; // success

}

int32_t SymPoly::addTerm(SymTerm& A) {

    if (A.factor.isZero() > 0) return 0; // nothing to add, but no error

    // ret 0 => success, !=0 = >error
    if (provideMemory1more() <= 0) return -1;

    int32_t error=0;

    if (useLex > 0) {
        int32_t insort=-1;
        int8_t present=findTermLexPosition_TA(insort,A);

        // present > 0 => add A to term[pos]
        // and remove zero if necessary
        if (present > 0) {
            BigInt w;
            if (bigintAdd_TAB(w,terms[sortidx[insort]].factor,A.factor) != 0) {
                printf("\n  error. bigint. probably overflow.");
                exit(99);
            }

            if (w.isZero() > 0) {
                error += lexRemoveTermAtSort(insort);
            } else {
                terms[sortidx[insort]].factor.copyFrom( w );
            }

            return error;

        } // present
        else {
            // present <= 0
            // new term: set at position pos
            // (pos=0..anzterm. if anzterm => add at end
            // and move pos+1 on to the right
            // memory enough to store an additional element

            for(int32_t k=anzterms;k>insort;k--) {
                sortidx[k]=sortidx[k-1];
            } // k
            terms[anzterms].copyFrom( A );
            sortidx[insort]=anzterms;
            anzterms++;

        } // not present

        return 0; // succes

    } // useLex

    // does the term itself already exist ?
    int32_t idx=-1;
    for(int32_t k=0;k<anzterms;k++) {
        if (symterm_bcExponentVgl_AB(terms[k],A) == 0) {
            idx=k;
            break;
        }
    } // k

    if (idx < 0) {
        // add new term, memory sufficient
        // term not zero
        terms[anzterms].copyFrom(A);
        anzterms++;
    } else {
        // add the factors
        BigInt w;
        if (bigintAdd_TAB(w,terms[idx].factor,A.factor) != 0) {
            printf("\n  error. bigint. probably overflow.");
            exit(99);
        }

        if (w.isZero() > 0) {
            // term be removed
            if (idx != (anzterms-1) ) {
                terms[idx].copyFrom(terms[anzterms-1]);
            }

            anzterms--;

        } else {
            terms[idx].factor.copyFrom(w);
        }

    }

    return 0;

}

void SymPoly::copyFrom(SymPoly& A) {
    setToZero();
    setMemory(A.anzterms+SYMTERMSPERBLOCK);

    if ( (useLex > 0) && (A.useLex > 0) ) {
        // copy
        for(int32_t k=0;k<A.anzterms;k++) {
            terms[k].copyFrom( A.terms[k] );
            sortidx[k]=A.sortidx[k];
        } // k

        anzterms=A.anzterms;

        return;
    }

    if ( (useLex > 0) && (A.useLex <= 0) ) {
        setToZero();
        for(int32_t k=0;k<A.anzterms;k++) {
            // use addTerm to generate sortidx
            if (addTerm(A.terms[k] ) != 0) {
                printf("\nerror. copy azpoly\n");
                exit(99);
            }
        } // k

        return;

    }

    for(int32_t k=0;k<A.anzterms;k++) {
        terms[k].copyFrom( A.terms[k] );
    } // k

    anzterms=A.anzterms;

}

int8_t azterm_azcExponentVgl_AB(
    AZterm& A,
    AZterm& B
) {
    // -1: exponents A < B

    if (A.cexponent < B.cexponent) return -1;
    if (A.cexponent > B.cexponent) return +1;

    for(int32_t k=0;k<ANZSYMVARS;k++) {
        if (A.azexponent[k] < B.azexponent[k]) return -1;
        if (A.azexponent[k] > B.azexponent[k]) return +1;
    } // k

    return 0;

}

int8_t symtermAllVgl_AB(SymTerm& A,SymTerm& B) {
    int8_t vgl=bigintVgl_AB(A.factor,B.factor);

    if (vgl != 0) return vgl;

    return symterm_bcExponentVgl_AB(A,B);
}

int8_t symterm_bcExponentVgl_AB(
    SymTerm& A,
    SymTerm& B
) {
    // -1: exponents A < B

    if (A.cexponent < B.cexponent) return -1;
    if (A.cexponent > B.cexponent) return +1;

    for(int32_t k=0;k<PERIOD;k++) {
        if (A.bexponent[k] < B.bexponent[k]) return -1;
        if (A.bexponent[k] > B.bexponent[k]) return +1;
    } // k

    return 0;

}

int32_t termMul_TAB(
    SymTerm& erg,
    SymTerm& A,
    SymTerm& B
) {
    // ret 0 => success, !=0 = >error
    if (bigintMul_TAB(erg.factor,A.factor,B.factor) != 0) {
        printf("\n  error. bigint. probably overflow.");
        exit(99);
    }
    erg.cexponent = sum_int32t(	A.cexponent,B.cexponent );
    for(int32_t k=0;k<PERIOD;k++) {
        erg.bexponent[k] = sum_int32t( A.bexponent[k], B.bexponent[k] );
    }

    return 0;

}

int32_t polyAdd_TAB(
    SymPoly& erg,
    SymPoly& A,
    SymPoly& B
) {
    // ret 0 => success, !=0 = >error
    erg.setToZero();
    int32_t m=maximumI(A.anzterms,B.anzterms);
    erg.setMemory( m+SYMTERMSPERBLOCK );
    erg.copyFrom(A);

    for(int32_t b=0;b<B.anzterms;b++) {
        if (erg.addTerm(B.terms[b]) != 0) return -1;
    }

    return 0;

}

int32_t polySub_TAB(
    SymPoly& erg,
    SymPoly& A,
    SymPoly& B
) {
    erg.setToZero();
    int32_t m=maximumI(A.anzterms,B.anzterms);
    erg.setMemory( m+SYMTERMSPERBLOCK );
    erg.copyFrom(A);

    for(int32_t b=0;b<B.anzterms;b++) {
        if (erg.subTerm(B.terms[b]) != 0) return -1;
    }

    return 0;

}

int32_t polyMul_TAB(
    SymPoly& erg,
    SymPoly& A,
    SymPoly& B
) {
    // ret 0 => success, !=0 = >error
    erg.setToZero();
    SymPoly tmp;
    tmp.setMemory( B.anzterms+SYMTERMSPERBLOCK );

    for(int32_t a=0;a<A.anzterms;a++) {
        if (polyMul_TTermA(tmp,A.terms[a],B) != 0) return -1; // error
        erg.addPoly(tmp);
    } // a

    return 0;

}

int32_t polyMul_TTermA(
    SymPoly& erg,
    SymTerm& A,
    SymPoly& B
) {
    // ret 0 => success, !=0 = >error
    erg.setToZero();
    erg.setMemory( B.anzterms+SYMTERMSPERBLOCK );
    SymTerm T;

    for(int32_t b=0;b<B.anzterms;b++) {
        if (termMul_TAB(T,A,B.terms[b]) != 0) return -1;
        if (erg.addTerm(T) != 0) return -1;
    } // b

    return 0;

}

// AZpoly operations
int32_t aztermMul_TAB(
    AZterm& erg,
    AZterm& A,
    AZterm& B
) {
    // ret 0 => success, !=0 = >error
    if (bigintMul_TAB(erg.factor,A.factor,B.factor) != 0) {
        printf("\n  error. bigint. probably overflow.");
        exit(99);
    }
    erg.cexponent = sum_int32t(A.cexponent,B.cexponent);
    for(int32_t k=0;k<ANZSYMVARS;k++) {
        erg.azexponent[k] = sum_int32t(A.azexponent[k],B.azexponent[k]);
    }

    return 0;

}

int32_t azpolyAdd_TAB(
    AZPoly& erg,
    AZPoly& A,
    AZPoly& B
) {
    // ret 0 => success, !=0 = >error
    erg.setToZero();
    int32_t m=maximumI(A.anzterms,B.anzterms);
    erg.setMemory( m+SYMTERMSPERBLOCK );
    erg.copyFrom(A);

    for(int32_t b=0;b<B.anzterms;b++) {
        if (erg.addTerm(B.terms[b]) != 0) return -1;
    }

    return 0;

}

int32_t azpolySub_TAB(
    AZPoly& erg,
    AZPoly& A,
    AZPoly& B
) {
    erg.setToZero();
    int32_t m=maximumI(A.anzterms,B.anzterms);
    erg.setMemory( m+SYMTERMSPERBLOCK );
    erg.copyFrom(A);

    for(int32_t b=0;b<B.anzterms;b++) {
        if (erg.subTerm(B.terms[b]) != 0) return -1;
    }

    return 0;

}

int32_t azpolyMul_TAB(
    AZPoly& erg,
    AZPoly& A,
    AZPoly& B
) {
    // ret 0 => success, !=0 = >error
    erg.setToZero();

    if (
        (A.anzterms <= 0) ||
        (B.anzterms <= 0)
    ) return 0; // success

    AZPoly tmp;
    int32_t sum=sum_int32t(A.anzterms, B.anzterms);
    int32_t sum2=sum_int32t(sum,AZTERMSPERBLOCK);
    tmp.setMemory( sum2 );

    for(int32_t a=0;a<A.anzterms;a++) {
        if (azpolyMul_TTermA(tmp,A.terms[a],B) != 0) return -1; // error
        erg.addPoly(tmp);
    } // a

    return 0;

}

int32_t azpolyMul_TTermA(
    AZPoly& erg,
    AZterm& A,
    AZPoly& B
) {
    // ret 0 => success, !=0 = >error
    erg.setToZero();

    if (A.factor.isZero() > 0) return 0; // zero

    erg.setMemory( B.anzterms+AZTERMSPERBLOCK );
    AZterm T;

    for(int32_t b=0;b<B.anzterms;b++) {
        if (aztermMul_TAB(T,A,B.terms[b]) != 0) return -1;
        // as terms[b] is distinct from all other, after
        // multiplying with T it still is, so
        // it can be set directly
        erg.terms[b].copyFrom( T );
        if ( (erg.useLex > 0) && (B.useLex > 0) ) {
            erg.sortidx[b]=B.sortidx[b]; // sorting does not change if after multiplication of ALL terms in a polynomial with the same AZterm
        }

    } // b
    erg.anzterms=B.anzterms;

    return 0;

}

int32_t createFullySymmetric_TA(
    SymPoly& erg,
    SymTerm& A
) {
    erg.setToZero();

    // order b-exponents of T
    SymTerm Tordered;
    Tordered.copyFrom(A);
    Tordered.sort_bExponents_dec();
    // how many non-zero exponents ?
    int32_t nze=0;
    for(int32_t k=0;k<PERIOD;k++) {
        if (Tordered.bexponent[k] != 0) nze++;
    }

    if (nze <= 0) {
        // constant 0-polynom
        erg.setToZero();
        return 0; // szccess
    }

    int32_t error=0;
    SymTerm one;

    // iterate
    int32_t loop[ANZSYMVARS];
    for(int32_t k=0;k<PERIOD;k++) loop[k]=0;
    int8_t first=1;

    while (1) {
        if (first <= 0) {
            // increment
            int8_t inc=0;
            for(int32_t k=0;k<nze;k++) {
                loop[k]++;
                if (loop[k] < PERIOD) {
                    inc=1;
                    break; // successfully set
                }

                // carry-over
                loop[k]=0;
            } // k

            if (inc <= 0) break; // all done

        } // first

        first=0;

        // is it a valid permutation ?
        int8_t valid=1;
        for(int32_t k1=0;k1<nze;k1++) {
            for(int32_t k2=(k1+1);k2<nze;k2++) {
                if (loop[k1] == loop[k2]) {
                    valid=0;
                    break;
                }
            } // k2

            if (valid <= 0) break;
        } // k1

        if (valid <= 0) continue;

        // build term
        one.setToZero();
        one.setFactor(Tordered.factor);
        one.cexponent=Tordered.cexponent;

        for(int32_t k=0;k<nze;k++) {
            one.set_bExponent_id_w(loop[k],Tordered.bexponent[k]);
        }

        if (erg.addTerm(one) != 0) return -1;

    } // while

    return error;

}

int8_t getFile_TF(DynSlowString& erg,const char* afn) {
    FILE* f=fopen(afn,"rt");
    if (!f) return 0;

    erg.setEmpty();

    while (!feof(f)) {
        char ch;
        ch=fgetc(f);
        if (
            (ch != ' ') &&
            (ch >= 32)
        ) erg.add(ch);

    } // while

    fclose(f);

    return 1;

}

void delfile(const char* afn) {
    char tt[4096];
    sprintf(tt,"del %s",afn);
    system(tt);
}

int8_t maximaexpandfactor_TA(
    DynSlowString& erg,
    DynSlowString& A

) {
    erg.setEmpty();
    erg.add(A);

    const char callfn[]="ef2.bat";

    FILE *f=fopen(callfn,"wt");
    if (!f) return 0;

    fprintf(f,"kill(all)$\nnumer:false$\ndisplay2d:false$\nttyoff:true$\n");
    fprintf(f,"poly:expand(%s)$\n",A.text);

    /*
    printf("\n\n\nATTN remove\n\n");
    fprintf(f,"poly:-poly$\n");
    printf("\n\n\nATTN remove\n\n");
    */

    fprintf(f,"\nttyoff:false$\ndisplay(poly);\n");

    fclose(f);

    const char ergfn[]="ef2.erg";

    // MAXIMAPATH
    char tt[2048];
    sprintf(tt,"call c:\\maxima-5.42.2\\bin\\maxima.bat <%s >%s",callfn,ergfn);
    system(tt);

    DynSlowString allfile;

    if (getFile_TF(allfile,ergfn) <= 0) {
        printf("\nerror. file2.\N");
        exit(99);
    }

    // case-sensitive
    // tstart at poly=
    // end at (% thereafter
    char* poly=strstr(allfile.text,"poly=");

    if (!poly) return 0;

    char* eofptr=strstr(poly,"(%");

    if (!eofptr) return 0;

    eofptr[0]=0;

    erg.setEmpty();
    int32_t pos=5; // after poly=
    // if there is a ( or ), fdiscard as this encloses beginning numbers
    while (poly[pos] != 0) {
        if (
            (poly[pos] == '(') ||
            (poly[pos] == ')')

        ) {
            pos++;
            continue;
        }

        erg.add(poly[pos]);
        pos++;
    }

    delfile(ergfn);
    delfile(callfn);

    return 1;

}

int8_t getTermFromString_TA(
    AZterm& erg,
    DynSlowString& A
) {
    // split at *
    int32_t last=-1;
    int32_t len=strlen(A.text);

    erg.setToZero();
    erg.setFactor(1); // in case of no number e.g. c^2

    DynSlowString sub;

    // <= as needs tor ead OVER the end for ik-1 to be a valöid right end
    for(int32_t k=0;k<=len;k++) {
        if (
            (A.text[k] == '*') ||
            (A.text[k] == 0)
        ) {
            // from last+1..
            if (last < 0) last=-1;

            sub.setEmpty();
            for(int32_t a=(last+1);a<=(k-1);a++) {
                sub.add(A.text[a]);
            }

            // what is sub ?
            if (
                (sub.text[0] == '-') ||
                (
                    (sub.text[0] >= '0') &&
                    (sub.text[0] <= '9')
                )
            ) {
                // number
                erg.factor.setstr(sub.text);
            } else {
                // variable
                char* p=strchr(sub.text,'^');
                int32_t expo=1;
                if (!p) expo=1; else expo=atoi(p+1);

                // CASE -sensitive
                if (sub.text[0] == 'c') {
                    erg.set_cExponent(expo);
                } else
                if (
                    (sub.text[0] >= 'A') &&
                    (sub.text[0] <= 'Z')
                ) {
                    erg.set_azExponent_id_w(sub.text[0]-'A',expo);
                } else {
                    printf("\nerror. variable.\n");
                    exit(99);
                }
            }

            last=k;
        }

    } // k

    return 1;
}

int8_t getPolyFromString_TA(
    AZPoly& erg,
    DynSlowString& inA2
) {
    erg.setToZero();

    DynSlowString inA;
    // rfemove spaces
    inA.setEmpty();
    if (!inA2.text) return 0; // error
    for(int32_t k=0;k<strlen(inA2.text);k++) {
        if (inA2.text[k] != ' ') inA.add(inA2.text[k]);
    }

    DynSlowString A;
    A.setEmpty();
    // remove superflupous signs
    int32_t al=strlen(inA.text);
    if (al <= 0) return 0; // error

    int32_t k=0;
    // [al] = 0 as string end, soi can be read

    // leading + be removed
    while (k < al) if (inA.text[k] == '+') k++; else break;

    // k is a t the first useable chaacter
    k--;

    while (k < al) {
        k++;

        // look ahead
        if (
            (inA.text[k] == '+') &&
            (inA.text[k+1] == '+')
        ) {
            A.add('+');
            k++; // jump over 2nd +
            continue;
        }

        if (
            (inA.text[k] == '+') &&
            (inA.text[k+1] == '-')
        ) {
            A.add('-');
            k++; // jump over 2nd
            continue;
        }

        if (
            (inA.text[k] == '-') &&
            (inA.text[k+1] == '+')
        ) {
            A.add('-');
            k++; // jump over 2nd
            continue;
        }

        if (
            (inA.text[k] == '-') &&
            (inA.text[k+1] == '-')
        ) {
            A.add('+');
            k++; // jump over 2nd
            continue;
        }

        A.add(inA.text[k]);

    }

    // split at + to get terms
    int32_t last=-1;
    int32_t len=strlen(A.text);
    DynSlowString term;
    AZterm T;

    // <= as needs tor ead OVER the end
    for(int32_t k=0;k<=len;k++) {
        if (
            (A.text[k] == '+') ||
            (A.text[k] == '-') ||
            (A.text[k] == 0)
        ) {

            if (k==0) {
                continue; // unary sign
            }

            // term from last+1..k-1
            term.setEmpty();
            if (last < 0) last=-1; // start at 0

            for(int32_t a=(last+1);a<=(k-1);a++) {
                term.add(A.text[a]);
            }

            if (getTermFromString_TA(T,term) <= 0) {
                printf("\nerror. cannot parse term\n");
                exit(99);
            }

            int32_t error=0;

            if (last < 0) error += erg.addTerm(T);
            else {
                if (A.text[last] == '+') error += erg.addTerm(T);
                else if (A.text[last] == '-') error += erg.subTerm(T);
            }

            if (error != 0) {
                printf("\nerror adding term\N");
                exit(99);
            }

            last=k;

        } // at +
    } // k

    return 1;
}

int32_t createHcPoly_TA(
    AZPoly& erg,
    const int32_t APER
) {
    // polynoimial in c for the centers of hyperbolic components
    // opf exact period APER
    int32_t error=0;

    AZPoly hcs[ANZSYMVARS+8];
    for(int32_t k=0;k<=APER;k++) {
        hcs[k].setToZero();
    }

    // as divisoon only works on the A-Z
    // varialbes but not on the seed LOWERCASE c
    // here UPPERCASE first variable is taken

    int32_t varidx=0;

    AZterm T;
    T.setToZero();
    T.setFactor(1);
    T.set_azExponent_id_w(varidx,1);

    error += hcs[1].addTerm(T); // c

    // build hc's of period PER and divisiors
    for(int32_t k=2;k<=APER;k++) {
        // hcs[k]=(hcs[k-1]^2+c)
        AZPoly cp;
        cp.copyFrom(hcs[k-1]);
        error += azpolyMul_TAB(hcs[k],hcs[k-1],cp);
        error += hcs[k].addTerm(T);

    } // k

    if (error != 0) return -1;

    struct Lokal {
        DivHlp hlp;
    };
    NEWLOKAL

    // divide out divisor-hcs
    for(int32_t k=2;k<=APER;k++) {
        // all below k have been constructed
        for(int32_t kv=1;kv<k;kv++) {
            if ( (k % kv) != 0) continue;

            AZPoly div,rem;
            error += azpolyDiv_TRABLH(div,rem,hcs[k],hcs[kv],varidx,var->hlp);
            if (error != 0) { RETL(-1) }
            if (rem.anzterms != 0) { printf("\nerror. expected remainder-.free"); exit(99); }

            hcs[k].copyFrom( div );

        } // kv

    } // kv

    // now the varialbe VARIXD has to be moved to LOWERCASE c
    erg.copyFrom( hcs[APER] );
    for(int32_t k=0;k<erg.anzterms;k++) {
        erg.terms[k].cexponent = erg.terms[k].azexponent[varidx];
        erg.terms[k].azexponent[varidx]=0;
    }

    return error;
}

int32_t getNumberOfCycles(const int32_t APER) {
    // for f)=z^2+c in APER-fold iteration
    // and f[N]=f(f(...(z)) N-fold
    // determines the z-degree of the polynomial which describes
    // the exact period APER cycles

    // strict1: z^2+c-z=0                   => degree 2
    // strict2: (f[2]-z) / strict1          => degree = degree 2-fold - degree strict1
    // strict3: (f[3]-z) / strict1
    // strict4: (f[4]-z) / (strict1*strict2)

    //d egree N-fold = 2^N
    int32_t deg[ANZSYMVARS+8]; // indexed by period
    deg[0]=0;
    deg[1]=2;

    for(int32_t k=2;k<=APER;k++) {
        int32_t D=(int64_t)1 << k; // nfold-degree
        for(int32_t div=1;div<k;div++) {
            if ( (k % div) == 0) {
                D = D-deg[div];
            }
        }
        deg[k]=D;

    }

    // then the number of cycles is DEG[APER] / APER
    if ( (deg[APER] % APER) != 0) {
        printf("\nimplementation error. cycles.");
        exit(99);
    }

    return (deg[APER] / APER);

}

void factorizeExternally(
    EquationList& unclassified,
    const int32_t tgt0,
    Equation0& A
) {
    if (A.nonright.anzterms != 0) {
        printf("\n  error. not able to mfactor if rhs is non-zero\N");
        return;
    }

    const char ergfn[]="factor.erg";

    printf("\n  externally factoring id=%i (factor exponents are disregarded) ... ",A.id);

    // factorize A.azleft with maxima
    // factors are stored at ids tgt...
    // factors are stored singularily
    // (iof a factor F1 occurs F1^6, the 6 will be discarded)

    const char fn2[]="_factor.bat";
    const char callfn[]="_factormeta.bat";

    FILE *f=fopen(callfn,"wt");
    fprintf(f,"batch(\"%s\");\n",fn2);
    fclose(f);

    f=fopen(fn2,"wt");
    fprintf(f,"kill(all)$\nnumer:false$\ndisplay2d:false$\nttyoff:true$\n");
    fprintf(f,"poly:expand(\n");

    A.azleft.ausgabe(f,-1);

    fprintf(f,")$\nfa:factor( poly )$\nttyoff:false$\n");
    //fprintf(f,"writefile(\"%s\");\n",ergfn);
    fprintf(f,"display(fa);\n");
    //fprintf(f,"closefile();\n");

    fclose(f);
    char tt[2048];
    // MAXIMAPATH
    sprintf(tt,"call c:\\maxima-5.42.2\\bin\\maxima.bat <%s >%s",callfn,ergfn);
    system(tt);

    // fa =
    // parse
    DynSlowString allfile;
    if (getFile_TF(allfile,ergfn) <= 0) {
        printf("\nerror. file.\N");
        exit(99);
    }

    int32_t id=tgt0;

    // now look for "fa="
    char* fa=strstr(allfile.text,"fa=");
    if (!fa) {
        printf("\nerror. no fa");
        exit(99);
    }
    char *eofptr=strstr(fa,"(%");
    if (!eofptr) {
        printf("\nerror. no eofptr");
        exit(99);
    }

    eofptr[0]=0;

    Equation0 oneeq;
    SymPoly nullpoly;
    nullpoly.setToZero();

    //printf("\nfound: %s",fa);
    int32_t klammer=0;
    int32_t pos=0;
    int32_t startpos=-1;
    AZPoly az;

    DynSlowString factor,factorexpanded;

    while (fa[pos] != 0) {
        pos++;

        if (startpos < 0) {
            // look for (
            if (fa[pos] == '(') {
                klammer++;
                if (klammer == 1) {
                    // new factor
                    startpos=pos;
                }
            }
        } else {
            // in factor, look for )
            if (fa[pos] == ')') {
                klammer--;
                if (klammer == 0) {
                    // and of factor
                    // factor [startpos..pos]
                    factor.setEmpty();
                    for(int32_t k=startpos;k<=pos;k++) {
                        factor.add(fa[k]);
                    }

                    // expand factor
                    if (maximaexpandfactor_TA(factorexpanded,factor) <= 0) {
                        printf("\nerror. expanding factor.\n");
                        exit(99);
                    }

                    //printf("\nexpanded factor= %s",factorexpanded.text);

                    // contains only variables *, ^numbers and +
                    // split at ecvery + to get a term
                    // then split tersm at * to get numbers or variable exponent
                    if (getPolyFromString_TA(az,factorexpanded) <= 0) {
                        printf("\nerror. cannot convert string to poly\N");
                        exit(99);
                    }

                    //printf("\n  fexp= |%s|",factorexpanded.text);
                    printf("\n    |= ");
                    az.ausgabe(stdout,12);

                    //exit(99);

                    oneeq.setAll_ABC(id,az,nullpoly);
                    int32_t error=0;
                    error += oneeq.fullSimplify();
                    error += unclassified.addEquation(oneeq);

                    if (error != 0) {
                        printf("\n  error. adding.\n");
                        return;
                    }

                    id++;

                    startpos=-1;
                }
            }
        }
    } // while


    delfile(ergfn);
    delfile(callfn);
    delfile(fn2);

}

void literature(FILE* fdescr) {
    fprintf(fdescr,"\n\nLiterature:\n");
    fprintf(fdescr,"[1] B Blum Smith, S Coskey. The fundamental theorem on symmetric polynomials. History's first whiff of Galois tehory. 2010.\n");
    fprintf(fdescr,"[2] K Conrad. Symmetric functions.\n");
    fprintf(fdescr,"[3] D Giarrusso, Y Fisher. A parameterization of the period 3 hyperbolic component of the Mandelbrot set. Proc Am Math Soc 123(12) 1995.\n");
    fprintf(fdescr,"[4] A Brown. Equations for periodic solutions of a logisitc difference equation. J Austral Math Soc 23 1981.\n");
    fprintf(fdescr,"[5] J Stephenson, T Ridgway. Formulae for cycles in the MandelBrot set. I, II, III. Physics A 1992.\n");
    fprintf(fdescr,"[6] A Healy. Resultants, resolvents and the computation of Galois group.\n");
    fprintf(fdescr,"[7] EH Bareiss. Sylvester's identity and multistep integer-preserving Gaussean elimination. 1968\n");
}

int32_t periodAll_interactive(void) {

    #define EXITW(W) \
    {\
        if (error != 0) { fprintf(stderr,"\nerror. %i\n",W); return (W); }\
    }

    int32_t error=0;
    #define SWAPLIST(AA,BB) \
    {\
        EquationList* tmplist=(AA);\
        (AA)=(BB);\
        (BB)=tmplist;\
    }

    #define CREATELIST(VAR) \
        (VAR)=new EquationList;\
        if ( !(VAR) ) { fprintf(stderr,"\nerror. memory 3.\n"); exit(99); }\
        (VAR)->fastClear();

    #define FREELIST(VAR) \
        if (VAR) delete (VAR);\
        (VAR)=NULL;

    // build a system of equations
    EquationList *unclassified;
    CREATELIST(unclassified)

    int8_t donotprint=0;
    Solution *sol=new Solution[ANZSYMVARS];
    if (!sol) { fprintf(stderr,"\nerror. memory solutions.\n"); exit(99); }
    for(int32_t k=0;k<ANZSYMVARS;k++) sol[k].setToUnsolved();

    FILE *fskript=NULL;
    fdescr=NULL;
    char ltmp[2048],utmp[2048];
    int8_t addfss=0; // flag for fdescr
    int8_t addcycbin2=0,addps=0;

    sprintf(fnbase,"_md2_p%i",
        PERIOD
    );

    // interactive loop
    // menu
    while (1) {
        if (donotprint > 0) {
            donotprint=0;
        } else {
            printf("\n\nperiod=%i: SET of equations",PERIOD);
            printf(" (#=%i)",unclassified->anzeq);
            printf("\n");
            unclassified->ausgabeVars(stdout,OUTPUTLEN);
        }

        if (fskript) {
            if (fgets(ltmp,2000,fskript) <= 0) {
                // close it, EOF
                fclose(fskript);
                fskript=NULL;
                donotprint=1;
                continue;
            }

        } else {
            printf("\n> ");
            scanf("%2000s",ltmp);
        }

        chomp(ltmp);
        strcpy(utmp,ltmp); upper(utmp);
        if (strstr(utmp,"QUIT")) {
            if (fdescr) {
                literature(fdescr);

                fclose(fdescr);
            }

            fdescr=NULL;

            break;
        }

        if (utmp[0] == '/') {
            donotprint=1;
            continue; // comment, in skript mainly
        }

        if (strlen(utmp) <= 0) {
            donotprint=1;
            continue;
        }


        if (strstr(utmp,"SIMPLIFYOFF") == utmp) {
            AUTOMATICSIMPLIFY=0;
            printf("\n  automatic simplifications switched OFF");
            donotprint=1;
        } else
        if (strstr(utmp,"SIMPLIFYON") == utmp) {
            AUTOMATICSIMPLIFY=1;
            printf("\n  automatic simplifications switched ON");
            donotprint=1;
        } else
        if (strstr(utmp,"REPEATITER(INF") == utmp) {
            REPEATITER=-1; // infinity
        } else
        if (strstr(utmp,"REPEATITER(ONCE") == utmp) {
            REPEATITER=1; // once
        } else
        if (!strcmp(utmp,"DESCRIBEOFF")) {
            if (fdescr) {
                fprintf(fdescr,"\n\nFinal: The system consists of equations");
                int8_t allzero=1;
                for(int32_t k=0;k<unclassified->anzeq;k++) {
                    if (k != 0) fprintf(fdescr,",");
                    fprintf(fdescr," %i",unclassified->aseq[k].id);
                    if (unclassified->aseq[k].nonright.anzterms != 0) allzero=0;
                }
                if (allzero > 0) {
                    fprintf(fdescr,"\nwhere all right hand sides = 0 ");
                }
                fprintf(fdescr,"\nin the variables");
                char varl[ANZSYMVARS+8];
                unclassified->get_VariablesLeft(varl);
                int8_t virst=1;
                for(int32_t k=0;k<strlen(varl);k++) {
                    if (varl[k] == '-') continue;
                    if (virst <= 0) fprintf(fdescr,",");
                    virst=0;
                    fprintf(fdescr," %c",varl[k]);
                }

                literature(fdescr);

                fprintf(fdescr,"\n\n");

                fclose(fdescr);
            }
            fdescr=NULL;
            donotprint=1;
            printf("\n description file closed.\n");
            continue;
        } else
        if (!strcmp(utmp,"DESCRCRLF")) {
            if (fdescr) fprintf(fdescr,"\n");
            donotprint=1;
        }
        if (strstr(utmp,"DESCR(") == utmp) {
            if (fdescr) fprintf(fdescr,"%s\n",&ltmp[6]);
            donotprint=0;

        } else
        if (strstr(utmp,"CONSTRBINDIV(") == utmp) {
            // tgtid,a0,a1,adiff,BETASTRING
            int32_t tgtid,a0,a1,adiff;
            char cycorfull;
            char term[2048];
            if (sscanf(&utmp[13],"%i,%i,%i,%i,%c,%s",&tgtid,&a0,&a1,&adiff,&cycorfull,term) != 6) {
                printf("\n  wrong. entry. error.\n");
                donotprint=1;
                continue;
            }

            int8_t fully=0;
            if (cycorfull == 'F') fully=1;
            //
            Equation0 A;
            if (constructBetaBin_TI_prodA0A1B_FsumT(A,tgtid,a0,a1,adiff,fully,term) != 0) {
                printf("\n  wrong. construction.");
                donotprint=1;
                continue;
            }

            if (unclassified->addEquation(A) != 0) {
                printf("\n  wrong. adding equation.");
                donotprint=1;
                continue;
            }

            if (fdescr) {
                fprintf(fdescr,"\n%i) Equation %i is based on a binomial reduction of the cyclic function\n",tgtid,tgtid);
                fprintf(fdescr,"\n  %i",a1);
                fprintf(fdescr,"\n  PROD [ b_k^2 - b_(k+delta mod p)^2) ]");
                fprintf(fdescr,"\n  k=%i",a0);
                fprintf(fdescr,"\n\nwith delta=%i in this case and ",adiff);
                fprintf(fdescr,"which can be expressed in two ways");
                fprintf(fdescr,"\n\n  PROD [ (b_k+b_(k+delta mod p))*(b_k - b_(k+delta mod p)) ]\n");
                fprintf(fdescr,"\nand on the other hand after using the iterational relatiionship b_k^2=b_(k+1 mod p)-c\n");
                fprintf(fdescr,"to lower the degrees\n");
                fprintf(fdescr,"\n  PROD [ b_(k+1 mod p)-c - b_(k+delta+1 mod p) + c ]");
                fprintf(fdescr,"\n\nwhich equals as a cyclic function PROD (b_k-b_(k+delta mod p)");
                fprintf(fdescr,"\n\nHence both sides can be divided by the same product ");
                fprintf(fdescr,"which results in\n");
                fprintf(fdescr,"\n  1 = PROD [ b_k + b_(k+delta mod p) ]");


                if ( (cycorfull == 'F') ||(cycorfull == 'f') ) {
                    DynSlowString tstr;
                    tstr.setEmpty();
                    int32_t indanz=0;
                    DynSlowString one;
                    one.setEmpty();
                    for(int32_t k=0;k<=strlen(term);k++) {
                        if ( (term[k] == ',') || (term[k] == 0) ) {
                            if (one.text) {
                                if (strlen(one.text) > 0) {
                                    if (one.text[0] != '0') {
                                        char tt[200];
                                        sprintf(tt,"b_i%i^",indanz);
                                        tstr.add(tt);
                                        tstr.add(one);
                                        tstr.add("*");
                                    }
                                    indanz++;
                                    one.setEmpty();
                                }
                            }
                        } else one.add(term[k]);

                    }

                    int32_t sl=0;
                    if (tstr.text) sl=strlen(tstr.text);

                    if (sl > 0) {
                        tstr.text[sl-1]=0;
                        fprintf(fdescr,"\n\nThis equation is then multiplied by ");
                        fprintf(fdescr,"the fully symmetric function\n");
                        fprintf(fdescr,"\n  %i",PERIOD-1);

                        fprintf(fdescr,"\n  SUM %s",tstr.text);
                        fprintf(fdescr,"\n  0 <= i0");
                        for(int32_t k=1;k<indanz;k++) {
                            fprintf(fdescr," <> i%i",k);
                        }
                        fprintf(fdescr," <= %i",PERIOD-1);
                    }

                } else
                if ( (cycorfull == 'C') ||(cycorfull == 'c') ) {
                    DynSlowString tstr;
                    tstr.setEmpty();
                    int32_t indanz=0;
                    DynSlowString one;
                    one.setEmpty();
                    for(int32_t k=0;k<=strlen(term);k++) {

                        if ( (term[k] == ',') || (term[k] == 0) ) {
                            if (strlen(one.text) > 0) {
                                if (one.text) {
                                    if (one.text[0] != '0') {
                                        char tt[200];
                                        sprintf(tt,"b_(k+%i mod p)^",indanz);
                                        tstr.add(tt);
                                        tstr.add(one);
                                        tstr.add("*");
                                    }
                                    indanz++;
                                    one.setEmpty();
                                }
                            }
                        } else one.add(term[k]);

                    }

                    int32_t sl=0;
                    if (tstr.text) sl=strlen(tstr.text);

                    if (sl > 0) {
                        tstr.text[sl-1]=0;

                        fprintf(fdescr,"\n\nThis equation is then multiplied by ");
                        fprintf(fdescr,"the cyclic function\n");
                        fprintf(fdescr,"\n  %i",PERIOD-1);


                        fprintf(fdescr,"\n  SUM %s",tstr.text);
                        fprintf(fdescr,"\n  0 <= k <= %i",PERIOD-1);
                    }

                } else {
                    printf("\nerror cycsymbol.");
                    exit(99);
                }


                fprintf(fdescr,"\n\nPowers of roots are reduced by the iterational relation and the");
                fprintf(fdescr,"\nequation is expressed in elementary symmetric functions as far as possible, which leads to\n\n");
                A.ausgabe(fdescr,-1);
                fprintf(fdescr,"\n");

            }

        } else
        if (strstr(utmp,"DESCRIBEON(") == utmp) {
            if (fdescr) fclose(fdescr);
            fdescr=fopen(&utmp[11],"wt");
            if (!fdescr) {
                printf("\n  error file open.");
            } else {
                fprintf(fdescr,"Notation:\n- Roots of cycles are denoted by b_0, b_1, .. b_(period-1) and are taken to be distinct.\n");
                fprintf(fdescr,"- Elementary symmetric functions are denoted by capital letters A (sum b_i), B (sum b_i*b_k with i <> k) etc).\n");
                fprintf(fdescr,"- The seed in the iteration z^2+c is denoted by the lowcase letter 'c'.\n\n");
                fprintf(fdescr,"- simplifcation of equation involve dividing by the greatest common divisor and reducing the degree of\n  the seed if it occurs in all terms.\n\n");
            }
            donotprint=0;

        } else
        if (!strcmp(utmp,"EMPTY")) {
            unclassified->fastClear();
            for(int32_t k=0;k<ANZSYMVARS;k++) sol[k].setToUnsolved();

            continue;
        } else

        if (!strcmp(utmp,"SKRIPTENDE")) {
            fclose(fskript);
            fskript=NULL;
            donotprint=1;
            continue;
        } else
        if (strstr(utmp,"SKRIPT(") == utmp) {
            fskript=fopen(&utmp[7],"rt");
            // can be NULL, then no skript opened
            donotprint=1;
            continue;
        } else
        if (strstr(utmp,"OUTLEN(") == utmp) {
            int32_t a;

            if (sscanf(&utmp[7],"%i",&a) != 1) {
                printf("\n  wrong. entry. error.");
                donotprint=1;
                continue;
            }

            OUTPUTLEN=a;
            donotprint=0;

        } else
        if (strstr(utmp,"PERIOD(") == utmp) {
            int32_t aper;

            if (sscanf(&utmp[7],"%i",&aper) != 1) {
                printf("\n  wrong. entry. error.");
                donotprint=1;
                continue;
            }

            if ( (aper < MINPERIOD) || (aper > ANZSYMVARS) ) {
                printf("\n  error. currently only periods %i..%i supported.\n",
                    MINPERIOD,ANZSYMVARS);
                donotprint=1;
                continue;
            }

            // delete all equations
            unclassified->fastClear();
            for(int32_t k=0;k<ANZSYMVARS;k++) sol[k].setToUnsolved();

            PERIOD=aper;

            printf("\ninitializing elementary symmetric functions ... ");
            initSigma(PERIOD);

            if (fdescr) {
                fprintf(fdescr,"Goal:\nFlow of algebraic steps to construct relationships between A (sum of periodic points),\n");
                fprintf(fdescr,"the seed c and %c (product of periodic points) in the iteration case z^2+c for period %i\n",
                    elemsymname[PERIOD-1],PERIOD);
            }

            sprintf(fnbase,"_md2_p%i",
                PERIOD
            );

            continue;

        } else
        if (strstr(utmp,"PRINT(") == utmp) {
            int32_t aid;

            if (sscanf(&utmp[6],"%i",&aid) != 1) {
                printf("\nwrong entry error\n");
                donotprint=1;
                continue;
            }

            Equation0* A=unclassified->getEquationPtr_id(aid);

            if (!A) {
                printf("\nwrong entry error\n");
                donotprint=1;
                continue;
            }

            fprintf(flog,"\n\nprint equation id=%i\n",aid);
            A->ausgabe(flog,-1);
            fprintf(flog,"\n");
            fflush(flog);
            donotprint=1;

        } else
        if (
            (strstr(utmp,"ADDEQ(") == utmp) ||
            (strstr(utmp,"SUBEQ(") == utmp) ||
            (strstr(utmp,"MULEQ(") == utmp) ||
            (strstr(utmp,"DEFEQ(") == utmp)
        ) {
            int8_t add=0,sub=0,mul=0,def=0;
            if (strstr(utmp,"ADDEQ(") == utmp) add=1;
            else if (strstr(utmp,"SUBEQ(") == utmp) sub=1;
            else if (strstr(utmp,"MULEQ(") == utmp) mul=1;
            else if (strstr(utmp,"DEFEQ(") == utmp) def=1;

            // ADDEQ(tgtid,id0,factor0,cexp0,id1,factor1,cexp1");
            int32_t tgtid,id0,id1,factor0,factor1,cexp0,cexp1;
            int32_t pow=-1;
            if (def > 0) {
                if (sscanf(&utmp[6],"%i,%i,%i,%i,%i,%i,%i,%i",
                    &tgtid,&id0,&factor0,&cexp0,
                    &id1,&factor1,&cexp1,&pow) != 8
                ) {
                    if (sscanf(&utmp[6],"%i,%i,%i,%i,%i,%i,%i",
                        &tgtid,&id0,&factor0,&cexp0,
                        &id1,&factor1,&cexp1) != 7
                    ) {
                        printf("\n  wrong entry error\n");
                        donotprint=1;
                        continue;
                    } else pow=1;
                }
            } else {
                if (sscanf(&utmp[6],"%i,%i,%i,%i,%i,%i,%i",
                    &tgtid,&id0,&factor0,&cexp0,
                    &id1,&factor1,&cexp1) != 7
                ) {
                    printf("\n  wrong entry error\n");
                    donotprint=1;
                    continue;
                }
            }

            if (
                (cexp0 < 0) ||
                (cexp1 < 0)
            ) {
                printf("\n  wrong. entry. error. negative c-exponents not allowed.\n");
                donotprint=1;
                continue;
            }

            Equation0 tgt,*A,*B;
            A=unclassified->getEquationPtr_id(id0);
            B=unclassified->getEquationPtr_id(id1);

            if ( (!A) || (!B) || (id0 == id1) ) {
                printf("\n  unknown or identical equations\n");
                donotprint=1;
                continue;
            }

            AZterm azt0,azt1;
            SymTerm st0,st1;
            st0.setToZero(); st0.setFactor( factor0 ); st0.set_cExponent( cexp0 );
            st1.setToZero(); st1.setFactor( factor1 ); st1.set_cExponent( cexp1 );
            azt0.setToZero(); azt0.setFactor( factor0 ); azt0.set_cExponent( cexp0 );
            azt1.setToZero(); azt1.setFactor( factor1 ); azt1.set_cExponent( cexp1 );

            AZPoly az0,az1,az2;
            SymPoly sp0,sp1,sp2;

            // azt0*A.left + azt1*B.left = st0*A.right + st1*B.right
            error=0;
            error += azpolyMul_TTermA(az0,azt0,A->azleft);
            error += azpolyMul_TTermA(az1,azt1,B->azleft);
            error += polyMul_TTermA(sp0,st0,A->nonright);
            error += polyMul_TTermA(sp1,st1,B->nonright);

            if (add > 0) {
                error += azpolyAdd_TAB(az2,az0,az1);
                error += polyAdd_TAB(sp2,sp0,sp1);
            } else if (sub > 0) {
                error += azpolySub_TAB(az2,az0,az1);
                error += polySub_TAB(sp2,sp0,sp1);
            } else if (mul > 0) {
                error += azpolyMul_TAB(az2,az0,az1);
                error += polyMul_TAB(sp2,sp0,sp1);
            } else if (def > 0) {
                //printf("\npow is= %i (rhs %i,%i)",pow,sp0.anzterms,sp1.anzterms);
                // rhs must be zeroi currenlty
                if (
                    (sp0.anzterms != 0) ||
                    (sp1.anzterms != 0)
                ) error--;
                else {
                    AZPoly divisor;
                    if (azpolyPow_TNA(divisor,pow,az1) != 0) {
                        printf("\n  error. not able to compute power of equation");
                        donotprint=1;
                        continue;
                    }
                    AZPoly rem;
                    DivHlp h;
                    // only working if remainder-free, so leading variable not relevant
                    error += azpolyDiv_TRABLH(az2,rem,az0,divisor,0,h);
                    if (rem.anzterms != 0) {
                        printf("\nerror. not remainder-free");
                        error--;
                    }
                }
                sp2.setToZero();
            }

            if (error == 0) {
                tgt.setAll_ABC(tgtid,az2,sp2);
                error += tgt.fullSimplify();
                error += unclassified->addEquation( tgt );
            }

            if (error != 0) {
                printf("\n  error.\n");
                donotprint=1;
                continue;
            }

            if (fdescr) {
                fprintf(fdescr,"\n\n%i) For Equation %i, the following operation is performed\n",tgtid,tgtid);
                fprintf(fdescr,"\n  eq%i*(",id0);
                azt0.ausgabe(fdescr);
                fprintf(fdescr,"  )\n    ");
                if (add > 0) fprintf(fdescr,"+");
                else if (sub > 0) fprintf(fdescr,"-");
                else if (mul > 0) fprintf(fdescr,"*");
                else if (def > 0) fprintf(fdescr,"/ (if remainder-free and right equations sides are 0)");
                fprintf(fdescr,"\n  [eq%i*(",id1);
                azt1.ausgabe(fdescr);
                fprintf(fdescr,")]");

                if (def > 0) fprintf(fdescr,"^%i\n",pow);
                else fprintf(fdescr,"\n");

                fprintf(fdescr,"\nwhich leads (after possible simplification) to\n\n  eq");
                tgt.ausgabe(fdescr,-1);

            }

        } else
        if (strstr(utmp,"HIPOW(") == utmp) {
            int32_t aid;
            char varsym;

            if (sscanf(&utmp[6],"%i,%c",&aid,&varsym) != 2) {
                printf("\n  wrong entry error.");
                continue;
            }

            int32_t var=getElemSymNameIdx(varsym);
            Equation0* A=unclassified->getEquationPtr_id(aid);

            if (
                (!A) ||
                (var < 0) || (var >= ANZSYMVARS)
            ) {
                printf("\n  wrong entry error. equation number out of range.");
                continue;
            }

            int32_t hid=A->azleft.get_hipow(var);
            printf("\n  => degree(%c)= %i\n",elemsymname[var],hid);
            donotprint=1;
            continue;

        } else
        if (strstr(utmp,"HIPOWS(") == utmp) {
            char varsym;

            if (sscanf(&utmp[7],"%c",&varsym) != 1) {
                printf("\n  wrong entry error.");
                continue;
            }

            int32_t var=getElemSymNameIdx(varsym);

            if (
                (var < 0) || (var >= ANZSYMVARS)
            ) {
                printf("\n  wrong entry error. equation number out of range.");
                continue;
            }

            for(int32_t k=0;k<unclassified->anzeq;k++) {
                int32_t hid=unclassified->aseq[k].azleft.get_hipow(var);
                printf("\n  id=%i : %i",unclassified->aseq[k].id,hid);
            } // k

            donotprint=1;
            continue;

        } else
        if (strstr(utmp,"ADDFSS(") == utmp) {
            char t2[2048];
            int32_t aid,afactor,cexp;
            if (sscanf(&utmp[7],"%i,%i,%i,%s",&aid,&afactor,&cexp,t2) != 4) {
                printf("\n  wrong entry error.\n");
                donotprint=1;
                continue;
            }

            int32_t err=0;
            if (fdescr) {
                if (addfss > 0) {
                    fprintf(fdescr,"\n\n%i) For equation %i, construct the following function and express it.\n",aid,aid);
                    fprintf(fdescr,"in two ways as above.\n");
                } else {
                    fprintf(fdescr,"\n\n%i) For equation %i, construct the following function and first express\n",aid,aid);
                    fprintf(fdescr,"it in elementary symmetric functions as is [2] and secondly, apply\n");
                    fprintf(fdescr,"the iterational relationship (b_(i+1 mod %i)-c)^n=b_i^2n\n");
                    if (REPEATITER > 0) {
                        fprintf(fdescr,"%i-times ",REPEATITER);
                    } else {
                        fprintf(fdescr,"as often as possible ");
                    }

                    fprintf(fdescr,"express it as far as possible by",PERIOD);
                    fprintf(fdescr,"elementary symmetric functions\nand equate both expressions to arrive at the final equation.\n");
                    addfss=1;
                }
            }

            err += unclassified->makeAddEquation_fullSymmetric_IFCBs(aid,afactor,cexp,t2);

            Equation0* A=NULL;
            if (err == 0) {
                A=unclassified->getEquationPtr_id(aid);

                if (A) {
                    err += A->fullSimplify();
                    if (fdescr) {
                        fprintf(fdescr,"\n\n  eq");
                        A->ausgabe(fdescr,-1);
                        fprintf(fdescr,"\n");
                    }
                } else err--;

            }

            if (err != 0) {
                printf("\n  error. adding.\n");
                donotprint=1;
                continue;
            }

            fflush(fdescr);

        } else
        if (strstr(utmp,"VALIDATE(") == utmp) {
            int32_t aid;
            if (sscanf(&utmp[9],"%i",&aid) != 1) {
                printf("\n  wrong. entry. error.");
                donotprint=1;
                continue;
            }

            char varsym=elemsymname[PERIOD-1];
            int32_t var=getElemSymNameIdx(varsym);

            if ( (var < 0) || (var >= ANZSYMVARS) ) {
                printf("\n  error.s ymbol.");
                donotprint=1;
                continue;
            }

            Equation0* A=unclassified->getEquationPtr_id(aid);
            if (!A) {
                printf("\n  wrong. entry. error.");
                donotprint=1;
                continue;
            }

            if (A->nonright.anzterms != 0) {
                printf("\n  error. right-hand side msut be 0");
                donotprint=1;
                continue;
            }

            int32_t deg=A->azleft.get_hipow(var);
            // all other degrees must be 0
            int8_t other0=1;
            for(int32_t k=0;k<(PERIOD-1);k++) {
                int32_t kv=getElemSymNameIdx(elemsymname[k]);
                if (A->azleft.get_hipow(kv) != 0) {
                    other0=0;
                    break;
                }
            } // k

            if (other0 <= 0) {
                printf("\n  error. only variable %c is allowed to have a non-zero degree\n",varsym);
                donotprint=1;
                continue;
            }

            int32_t numberCyclesPerSeed=getNumberOfCycles(PERIOD);

            if (deg != numberCyclesPerSeed) {
                LOGMSG3("\n\n  FAILED. Number of cycles per seed in period %i is %i\n",PERIOD,numberCyclesPerSeed);
                LOGMSG3("  but %c-degree of equation is %i\n",varsym,deg);
                donotprint=1;
                continue;
            }

            AZPoly hc;
            AZPoly coeff0;
            AZPoly id;
            int32_t error=0;
            error += createHcPoly_TA(hc,PERIOD);
            if (error == 0) error += A->azleft.coeff_TID(coeff0,var,0);
            if (error == 0) error += azpolySub_TAB(id,hc,coeff0);

            if (error != 0) {
                printf("\n  error creating hyperbolic center polynomial");
                donotprint=1;
                continue;
            }

            if (id.anzterms != 0) {
                printf("\n  FAILED. %c^0-coefficient does NOT equal hc-polynomial",var);
                donotprint=1;
                continue;
            }

            donotprint=1;

            // all passed
            LOGMSG4("\n\nValidation passed:\n  - degree of %c equals the maximum number of cycles per seed\n  - the %c^0-coefficient equals the formula for the hyperbolic centers of period %i\n",
                varsym,varsym,PERIOD);

            if (fdescr) {
                fprintf(fdescr,"\n\n\nVALIDATION of equation %i:\n\n",aid);
                A->ausgabe(fdescr,-1);
                fprintf(fdescr,"\n\nFor a given seed c, the equation results in a polynomial only in %c.",varsym);
                fprintf(fdescr,"\nThat polynomial's %c-degree of %i equals the maximum number of cycles of ",varsym,deg);
                fprintf(fdescr,"length %i.\n",PERIOD);

                fprintf(fdescr,"\nThe hyperbolic centers of exact period %i are described by\n\n",PERIOD),

                hc.ausgabe(fdescr,-1);

                fprintf(fdescr,"\n\nwhich equals the coefficient of the %c^0-term, i.e. after setting the\n",varsym);
                fprintf(fdescr,"product of periodic points %c to 0 (= superattracting).\n",varsym);

            }
        }
        if (strstr(utmp,"SETSTR(") == utmp) {
            int32_t aid;
            char formel[4096];

            // CASE-sensitive
            if (sscanf(&ltmp[7],"%i,%s",&aid,formel) != 2) {
                printf("\nwrong entry error");
                donotprint=1;
                continue;
            }

            AZPoly az;
            DynSlowString d;
            d.setEmpty();
            d.add(formel);

            if (getPolyFromString_TA(az,d) <= 0) {
                printf("\n  error. formula.\n");
                donotprint=1;
                continue;
            }

            Equation0 eq2;
            SymPoly nullpoly;
            AZPoly aznull;
            aznull.setToZero();
            nullpoly.setToZero();
            // no simplification at setting
            eq2.setAll_ABC(aid,aznull,nullpoly);
            int32_t err=0;
            err += unclassified->addEquation(eq2);
            if (err != 0) {
                printf("\n  error. adding.");
                donotprint=1;
                continue;
            }
            Equation0* A=unclassified->getEquationPtr_id(aid);
            if (!A) {
                printf("\n  error. adding.");
                donotprint=1;
                continue;
            }
            // now COPY
            A->azleft.copyFrom( az );

            if (fdescr) {
                fprintf(fdescr,"\n\n%i) Equation %i is directly set to\n\n",aid,aid);
                A->ausgabe(fdescr,-1);
                fprintf(fdescr,"\n");
            }

            donotprint=0;


        } else
        if (strstr(utmp,"DEL(") == utmp) {
            int32_t aid,aid2;

            if (sscanf(&utmp[4],"%i,%i",&aid,&aid2) == 2) {
            } else
            if (sscanf(&utmp[4],"%i",&aid) == 1) {
                aid2=aid;
            } else {
                printf("\n  wrong entry error.");
                continue;
            }

            for(int32_t k=aid;k<=aid2;k++) {
                // return discard as not all id's must be there
                (void)unclassified->removeEquationById(k);
            }

        } else
        if (strstr(utmp,"RESTABV(") == utmp) {
            int32_t tgtid,srcid,fin;
            char varsym;

            if (sscanf(&utmp[8],"%i,%i,%i,%c",&fin,&tgtid,&srcid,&varsym) != 4) {
                printf("\n  wrong entry error.");
                continue;
            }

            int32_t var=getElemSymNameIdx(varsym);
            Equation0* tgt=unclassified->getEquationPtr_id(tgtid);
            Equation0* src=unclassified->getEquationPtr_id(srcid);

            if (
                (!tgt) || (!src) ||
                (tgtid == srcid) ||
                (var < 0) || (var >= ANZSYMVARS)
            ) {
                printf("\n  wrong entry error. equation number out of range.");
                continue;
            }

            if (
                (tgt->nonright.anzterms != 0) ||
                (src->nonright.anzterms != 0)
            ) {
                printf("\n  error. only empty rhs is allowed.");
                donotprint=1;
                continue;
            }

            AZPoly res;
            int32_t clock0=clock();
            int32_t ret;

            ret=fractionFreeGaussMultistep2_TABI(
                res,
                tgt->azleft,
                src->azleft,
                var);

            if (ret != 0) {
                printf("\nerror. resultant\n");
                donotprint=1;
                continue;
            }

            int32_t clock1=clock();
            double dtime=clock1-clock0;
            dtime /= CLOCKS_PER_SEC;
            int32_t dclock=(int32_t)ceil(dtime);
            LOGMSG2("\n  %i sec",dclock);

            Equation0 eq;
            eq.setAll_ABC(fin,res,tgt->nonright);
            if (
                (eq.fullSimplify() != 0)
            ) {
                printf("\nerror reducing\n");
                donotprint=1;
                continue;
            }
            error += unclassified->addEquation(eq);

            donotprint=1;

            if (fdescr) {
                fprintf(fdescr,"\n\n%i) For equation %i, the resultant of eq%i and eq%i with respect to symbol %c\n",
                    fin,fin,tgtid,srcid,elemsymname[var]);
                fprintf(fdescr,"is computed using Bareiss' %i-step algorithm [7].\n\n  eq",BAREISS_STEP);
                eq.ausgabe(fdescr,-1);
            }

            continue;
        } else
        if (strstr(utmp,"XFACTOR2(") == utmp) {
            int32_t tgt,src;
            if (sscanf(&utmp[9],"%i,%i",&tgt,&src) != 2) {
                printf("\n  wrong. entry. error1.");
                donotprint=1;
                continue;
            }

            Equation0* A=unclassified->getEquationPtr_id(src);

            if (
                (!A)
            ) {
                printf("\n  wrong. entry. error2.");
                donotprint=1;
                continue;
            }

            factorizeExternally(*unclassified,tgt,*A);
            donotprint=0;

        } else
        if (strstr(utmp,"RESTABVALL(") == utmp) {
            int32_t offset,aid;
            char varsym;

            if (sscanf(&utmp[11],"%c,%i,%i",&varsym,&aid,&offset) != 3) {
                printf("\n  wrong entry error.");
                continue;
            }

            int32_t var=getElemSymNameIdx(varsym);
            Equation0* wrt=unclassified->getEquationPtr_id(aid);

            if (
                (!wrt) ||
                (var < 0) || (var >= ANZSYMVARS)
            ) {
                printf("\n  wrong entry error. equation number out of range.");
                continue;
            }

            if (wrt->nonright.anzterms != 0) {
                printf("\n  error. right hand side must be zero.");
                donotprint=1;
                continue;
            }


            AZPoly res;
            int32_t ret;
            int32_t anz=unclassified->anzeq;
            int32_t* idx=new int32_t[anz];

            for(int32_t k=0;k<anz;k++) {
                idx[k]=unclassified->aseq[k].id;
            }

            for(int32_t k=0;k<anz;k++) {
                if (idx[k] == aid) continue;

                Equation0* A=unclassified->getEquationPtr_id(idx[k]);
                if (!A) continue;

                if (A->nonright.anzterms != 0) continue;

                int32_t clock0=clock();

                ret=fractionFreeGaussMultistep2_TABI(
                    res,
                    A->azleft,
                    wrt->azleft,
                    var);

                if (ret != 0) {
                    printf("\nerror. resultant\n");
                    donotprint=1;
                    continue;
                }

                int32_t clock1=clock();
                double dtime=clock1-clock0;
                dtime /= CLOCKS_PER_SEC;
                int32_t dclock=(int32_t)ceil(dtime);
                LOGMSG2("\n  %i sec",dclock);

                Equation0 eq;
                eq.setAll_ABC(idx[k]+offset,res,A->nonright);

                if (
                    (eq.fullSimplify() != 0)
                ) {
                    printf("\nerror reducing\n");
                    donotprint=1;
                    continue;
                }
                error += unclassified->addEquation(eq);

                if (fdescr) {
                    fprintf(fdescr,"\n\n%i) For equation %i, the resultant of eq%i and eq%i with respect to symbol %c\n",
                        idx[k]+offset,idx[k]+offset,idx[k],aid,elemsymname[var]);
                    fprintf(fdescr,"is computed using Bareiss' %i-step algorithm [7].\n\n  eq",BAREISS_STEP);

                    eq.ausgabe(fdescr,-1);
                }

                if (unclassified->removeEquationById(idx[k]) != 0) {
                    printf("\n  error deleting %i.\n",idx[k]);
                    donotprint=1;
                    continue;

                }

            } // k

            donotprint=0;
            delete[] idx;

        } else
        if (strstr(utmp,"SAVE(") == utmp) {
            FILE *f=fopen(&utmp[5],"wb");
            if (!f) { printf("\nerror.\n"); donotprint=1; continue; }

            unclassified->save(f);
            WRITE32(f,ANZSYMVARS)
            for(int32_t k=0;k<ANZSYMVARS;k++) {
                sol[k].save(f);
            }

            fclose(f);
            printf("\n  saved.");
            donotprint=1;

        } else
        if (strstr(utmp,"LOAD(") == utmp) {
            FILE *f=fopen(&utmp[5],"rb");
            if (!f) { printf("\nerror.\n"); donotprint=1; continue; }

            unclassified->load(f);
            int32_t w32;
            READ32(f,w32)
            if (w32 != ANZSYMVARS) {
                printf("\nerror. symbol number not consistent from file stored to  compiled\n");
                exit(99);
            }

            for(int32_t k=0;k<w32;k++) {
                sol[k].load(f);
            }

            fclose(f);
        } else
        if (strstr(utmp,"LOADEQ(") == utmp) {
            int32_t aid;
            char fn[1024];

            if (sscanf(&utmp[7],"%i,%s",&aid,fn) != 2) {
                printf("\n  wrong entry error.");
                continue;
            }

            FILE *f=fopen(fn,"rb");
            if (!f) {
                printf("\nerror. file.\n");
                donotprint=1;
                continue;
            }

            Equation0 E;
            if (E.load(f) <= 0) {
                printf("\nerror. file.\n");
                donotprint=1;
                continue;
            }
            E.id=aid;

            error=0;
            error += unclassified->addEquation( E );
            if (error != 0) {
                printf("\nerror. adding.\n");
                donotprint=1;
                continue;
            }

            fclose(f);
            donotprint=0;

            if (fdescr) {
                fprintf(fdescr,"\n\n%i) Equation %i has been directly set to\n\n",aid,aid);
                E.ausgabe(fdescr,-1);
                fprintf(fdescr,"\n");
            }

        } else
        if (strstr(utmp,"SAVEEQ(") == utmp) {
            int32_t aid;
            char fn[1024];

            if (sscanf(&utmp[7],"%i,%s",&aid,fn) != 2) {
                printf("\n  wrong entry error.");
                continue;
            }

            Equation0* A=unclassified->getEquationPtr_id(aid);
            if (!A) {
                printf("\n  error. equation.\n");
                donotprint=1;
                continue;
            }

            FILE *f=fopen(fn,"wb");
            if (!f) {
                printf("\nerror. file.\n");
                donotprint=1;
                continue;
            }

            if (A->save(f) <= 0) {
                printf("\nerror. file.\n");
                donotprint=1;
                continue;
            }

            fclose(f);
            donotprint=1;
            printf("\n  saved.");

        } else
        if (strstr(utmp,"SOLT(") == utmp) {
            int32_t tgtid,srcid;
            char varsym;

            if (sscanf(&utmp[5],"%i,%i,%c",&tgtid,&srcid,&varsym) != 3) {
                printf("\n  wrong entry error.");
                continue;
            }

            int32_t var=getElemSymNameIdx(varsym);
            Equation0* A=unclassified->getEquationPtr_id(srcid);

            if (
                (!A) || (tgtid == srcid) ||
                (var < 0) ||
                (var >= ANZSYMVARS)
            ) {
                printf("\n  wrong entry error. equation number out of range.");
                continue;
            }

            if (A->nonright.anzterms != 0) {
                printf("\n  error. only empty rhs is allowed.");
                donotprint=1;
                continue;
            }

            donotprint=1;

            AZRational erg;
            if (substitute_TAB(erg,A->azleft,sol[var]) != 0) {
                printf("\nerror substituing\n");
                continue;
            }

            // if the original equation is = 0
            // NENNER can be ignored
            // now replace equation in list with the
            // NUMERATOR, as equations are = 0 => so denominator
            // will be removed

            Equation0 E;
            SymPoly sym;
            sym.setToZero();

            E.setAll_ABC(tgtid,erg.zaehler,sym);
            error=0;
            error += E.fullSimplify();
            error += unclassified->addEquation( E );

            if (error != 0) {
                printf("\nerror. reduction.\n");
                exit(99);
            }

            if (fdescr) {
                fprintf(fdescr,"\n\n%i) For equation %i, the solution rational for symbol %c is substituted into equation %i\nwhich, after bringing every term to the same denominator and multiplying to get a polynomial, leads to\n\n  eq",tgtid,tgtid,elemsymname[var],srcid);
                E.ausgabe(fdescr,-1);
            }

        } else
        if (strstr(utmp,"SOLTALL(") == utmp) {
            int32_t offset,srcid;
            char varsym;

            if (sscanf(&utmp[8],"%c,%i",&varsym,&offset) != 2) {
                printf("\n  wrong entry error.");
                continue;
            }

            int32_t var=getElemSymNameIdx(varsym);

            if (
                (offset <= 0) ||
                (var < 0) ||
                (var >= ANZSYMVARS)
            ) {
                printf("\n  wrong entry error. equation number out of range.");
                continue;
            }

            donotprint=1;

            int32_t anz=unclassified->anzeq;
            int32_t* idx=new int32_t[anz];

            for(int32_t k=0;k<anz;k++) {
                idx[k]=unclassified->aseq[k].id;
            }

            for(int32_t k=0;k<anz;k++) {
                Equation0* A=unclassified->getEquationPtr_id(idx[k]);
                if (!A) continue;

                if (A->nonright.anzterms != 0) {
                    continue;
                }

                AZRational erg;
                if (substitute_TAB(erg,A->azleft,sol[var]) != 0) {
                    printf("\nerror substituing\n");
                    continue;
                }

                Equation0 E;
                SymPoly sym;
                sym.setToZero();

                E.setAll_ABC(idx[k]+offset,erg.zaehler,sym);
                error=0;
                error += E.fullSimplify();
                error += unclassified->addEquation( E );

                if (error != 0) {
                    printf("\nerror. reduction.\n");
                    exit(99);
                }

                if (fdescr) {
                    fprintf(fdescr,"\n\n%i) For equation %i, the solution rational for symbol %c is substituted into equation %i\nwhich, after bringing every term to the same denominator and multiplying to get a polynomial, leads to\n\n  eq",idx[k]+offset,idx[k]+offset,elemsymname[var],srcid);
                    E.ausgabe(fdescr,-1);
                }

                // now delete A
                if (unclassified->removeEquationById(idx[k]) != 0) {
                    printf("\n  error. deleting\n");
                    exit(99);
                }

            } // k

            delete[] idx;
            donotprint=0;

        } else
        // solving linear
        if (strstr(utmp,"SOLVE(") == utmp) {
            int32_t aid;
            char varsym;
            if (sscanf(&utmp[6],"%i,%c",&aid,&varsym) != 2) {
                printf("\n  wrong entry error.");
                continue;
            }

            int32_t var=getElemSymNameIdx(varsym);
            Equation0* A=unclassified->getEquationPtr_id(aid);

            if (
                (!A)||
                (var < 0) ||
                (var >= ANZSYMVARS)
            ) {
                printf("\n  wrong entry error. equation number out of range.");
                continue;
            }

            if (A->nonright.anzterms != 0) {
                printf("\n  error. only empty rhs is allowed.");
                donotprint=1;
                continue;
            }

            donotprint=1;

            AZRational solve;
            if (A->solveForLinear_TI(solve,var) != 0) {
                printf("\n  error solving.\n");
                continue;
            }

            sol[var].setSolution_id(var,solve);
            //sol[var].ausgabe(stdout,OUTLEN);

            if (fdescr) {
                fprintf(fdescr,"\n\n--) Equation %i is solved for symbol %c which leads to\n\n  ",aid,elemsymname[var]);
                sol[var].ausgabe(fdescr,-1);
            }

        } else
        // simple or perfect simple variable
        if (strstr(utmp,"SIMP(") == utmp) {
            int32_t aid;
            if (sscanf(&utmp[5],"%i",&aid) != 1) {
                printf("\n  wrong entry error.");
                continue;
            }

            Equation0* E=unclassified->getEquationPtr_id(aid);
            if (!E) {
                printf("\n  wrong entry error. equation number out of range.");
                continue;
            }

            // indexed by elemsymname
            for(int32_t var=0;var<PERIOD;var++) {
                int32_t hid=E->azleft.get_hipow(var);
                int8_t degree1=0,coeffsimple=0,perfect=0;
                if (hid != 1) continue;

                degree1=1;

                #define OUTCHARACT \
                {\
                    printf("\n  %c: ",elemsymname[var]);\
                    if (perfect > 0) printf("perfect");\
                    else if (coeffsimple > 0) printf("simple");\
                    else if (degree1 > 0) printf("linear");\
                }

                AZPoly coeff;
                if (E->azleft.coeff_TID(coeff,var,hid) != 0) {
                    printf("\nerror. classification.\n");
                    break;
                }

                // characterize coefficients
                if (coeff.predSimple() <= 0) {
                    OUTCHARACT
                    continue;
                }

                coeffsimple=1;

                if (E->azleft.nohigherthan(var) > 0) perfect=1;

                OUTCHARACT

            } // var
            donotprint=1;

        } else
        // coeff
        if (strstr(utmp,"COEFF(") == utmp) {
            int32_t aid,varidx,varexp;
            char varsymbol;
            if (sscanf(&utmp[6],"%i,%c,%i",&aid,&varsymbol,&varexp) != 3) {
                printf("  wrong entry error.\n");
                continue;
            }
            varidx=getElemSymNameIdx(varsymbol);
            Equation0* A=unclassified->getEquationPtr_id(aid);

            if (
                (!A) ||
                (varidx < 0)
            ) {
                printf("\nwrong entry error.\n");
                donotprint=1;
                continue;
            }

            AZPoly coeff;
            if (A->azleft.coeff_TID(coeff,varidx,varexp) != 0) {
                printf("\nerror. coefficient.\n");
                donotprint=1;
                continue;
            }

            printf("\ncoefficient /");
            char var[ANZSYMVARS + 1 +8];
            coeff.getVariables(var);
            printf("%s/ = ",var);
            coeff.ausgabe(stdout,32);
            donotprint=1;
        } else
        if (strstr(utmp,"COEFFS(") == utmp) {
            int32_t varidx,varexp;
            char varsymbol;
            if (sscanf(&utmp[7],"%c,%i",&varsymbol,&varexp) != 2) {
                printf("  wrong entry error.\n");
                continue;
            }
            varidx=getElemSymNameIdx(varsymbol);

            if (
                (varidx < 0)
            ) {
                printf("\nwrong entry error.\n");
                donotprint=1;
                continue;
            }

            for(int32_t k=0;k<unclassified->anzeq;k++) {
                Equation0* A=&unclassified->aseq[k];;

                AZPoly coeff;
                if (A->azleft.coeff_TID(coeff,varidx,varexp) != 0) {
                    printf("\nerror. coefficient.\n");
                    donotprint=1;
                    continue;
                }

                printf("\n[id=%i] => /",A->id);
                char var[ANZSYMVARS + 1 +8];
                coeff.getVariables(var);
                printf("%s/ = ",var);
                coeff.ausgabe(stdout,32);
            } // k

            donotprint=1;
        } else
        if (strstr(utmp,"ELMT(") == utmp) {
            int32_t tgtid,var,varexp,id0,id1;
            char varsym;
            if (sscanf(&utmp[5],"%i,%i,%c,%i,%i",&tgtid,&id0,&varsym,&varexp,&id1) != 5) {
                printf("  wrong entry error.\n");
                continue;
            }

            Equation0* A0=unclassified->getEquationPtr_id(id0);
            Equation0* A1=unclassified->getEquationPtr_id(id1);

            var=getElemSymNameIdx(varsym);

            if (
                (!A0) || (!A1) || (id0 == id1) ||
                (id0 == tgtid) || (id1 == tgtid) ||
                (var < 0)
            ) {
                printf("\nwrong entry error.\n");
                donotprint=1;
                continue;
            }

            // tgt=coeff(id0,sym^exp)*id1 - coeff(id1,sym^exp)*id0
            // reduce
            if (
                (A0->nonright.anzterms > 0) ||
                (A1->nonright.anzterms > 0)
            ) {
                donotprint=1;
                printf("\n  error. eliminate only works with em,empty rrhs.\n");
                continue;
            }

            AZPoly eliminate,coeffA,coeffB;
            if (eliminate_TTTIEAB(
                    eliminate,
                    coeffA,coeffB,
                    var,varexp,
                    A0->azleft,
                    A1->azleft
                ) != 0
            ) {
                printf("\nerror. eliminating.\n");
                exit(99);
            }

            Equation0 E;
            SymPoly sym;
            sym.setToZero();
            E.setAll_ABC(tgtid,eliminate,sym);

            error=0;
            error += E.fullSimplify();
            error += unclassified->addEquation(E);

            if (error != 0) {
                printf("\n  error. adding.\n");
                donotprint=1;
                continue;
            }

            if (fdescr) {
                fprintf(fdescr,"\n\n%i) For equation %i, terms harbouring %c^%i are eliminated by subtracting multiples of\n",tgtid,tgtid,elemsymname[var],varexp);
                fprintf(fdescr,"equations %i and %i:\n",id0,id1);
                fprintf(fdescr,"\n  eq%i*(",id0);
                coeffB.ausgabe(fdescr,-1);
                fprintf(fdescr,")\n    -\n  eq%i(",id1);
                coeffA.ausgabe(fdescr,-1);
                fprintf(fdescr,")\n\nwhich leads to\n\n  eq");
                E.ausgabe(fdescr,-1);
            }

        } // eliminate
        else
        if (strstr(utmp,"ELMTALL(") == utmp) {
            int32_t offset,var,varexp,aid;
            char varsym;
            if (sscanf(&utmp[8],"%c,%i,%i,%i",&varsym,&varexp,&aid,&offset) != 4) {
                printf("  wrong entry error.\n");
                continue;
            }

            Equation0* A0=unclassified->getEquationPtr_id(aid);
            var=getElemSymNameIdx(varsym);

            if (
                (offset <= 0) ||
                (!A0) ||
                (var < 0)
            ) {
                printf("\nwrong entry error.\n");
                donotprint=1;
                continue;
            }

            // tgt=coeff(id0,sym^exp)*id1 - coeff(id1,sym^exp)*id0
            // reduce
            if (
                (A0->nonright.anzterms > 0)
            ) {
                donotprint=1;
                printf("\n  error. eliminate only works with em,empty rrhs.\n");
                continue;
            }

            int32_t anz=unclassified->anzeq;
            int32_t* idx=new int32_t[anz];

            for(int32_t k=0;k<anz;k++) {
                idx[k]=unclassified->aseq[k].id;
            }

            for(int32_t k=0;k<anz;k++) {
                if (idx[k] == aid) continue;

                Equation0* A1=unclassified->getEquationPtr_id(idx[k]);

                if (A1->nonright.anzterms != 0) {
                    continue;
                }

                AZPoly eliminate,coeffA,coeffB;
                if (eliminate_TTTIEAB(
                        eliminate,
                        coeffA,coeffB,
                        var,varexp,
                        A0->azleft,
                        A1->azleft
                    ) != 0
                ) {
                    printf("\nerror. eliminating.\n");
                    exit(99);
                }

                Equation0 E;
                SymPoly sym;
                sym.setToZero();
                E.setAll_ABC(idx[k]+offset,eliminate,sym);

                error=0;
                error += E.fullSimplify();
                error += unclassified->addEquation(E);

                if (error != 0) {
                    printf("\n  error. adding.\n");
                    donotprint=1;
                    continue;
                }

                if (fdescr) {
                    fprintf(fdescr,"\n\n%i) For equation %i, terms harbouring %c^%i are eliminated by subtracting multiples of\n",idx[k]+offset,idx[k]+offset,elemsymname[var],varexp);
                    fprintf(fdescr,"equations %i and %i:\n",aid,idx[k]);
                    fprintf(fdescr,"\n  eq%i*(",aid);
                    coeffB.ausgabe(fdescr,-1);
                    fprintf(fdescr,")\n    -\n  eq%i(",idx[k]);
                    coeffA.ausgabe(fdescr,-1);
                    fprintf(fdescr,")\n\nwhich leads to\n\n  eq");
                    E.ausgabe(fdescr,-1);
                }

                // delete A1
                if (unclassified->removeEquationById(idx[k]) != 0) {
                    printf("\n  error deleting %i.\n",idx[k]);
                    donotprint=1;
                    continue;

                }

            } // k

            delete[] idx;

        } // eliminate


    } // while menu

    FREELIST(unclassified)

    return 0; // no error

}

int32_t replace_repeated_iter(
    SymPoly& A
) {
    SymPoly cpA;
    cpA.copyFrom(A);

    int32_t error=0;
    error += replace_repeated_iter_TA(A,cpA);
    if (error == 0) return 0; // success

    A.copyFrom(cpA); // leave unchanged
    return -1;

}

int32_t replace_repeated_iter_TA(
    SymPoly& erg,
    SymPoly& A
) {
    erg.setToZero();
    SymPoly terg;
    erg.copyFrom(A);

    int32_t error=0;
    int32_t todoiter=REPEATITER;

    while (erg.get_maxBexponent() > 1) {
        //printf("\nerg%i= ",todoiter);
        //erg.ausgabe(stdout,-1);

        if (REPEATITER > 0) {
            todoiter--;
            if (todoiter < 0) {
                printf("\nout");
                break;
            }
        }

        error += replace_iter_TA(terg,erg);
        if (error != 0) return -1;
        erg.copyFrom(terg);

    } // while

    return error;
}

int32_t replace_iter_TA(
    SymPoly& erg,
    SymPoly& A
) {
    // ret != 0 => error
    // = 0 => success

    // replace all b_i^2 with b_(i+1)-c (i+1 mod ANZVAR)
    // bNext = b^2+c

    erg.setToZero();
    erg.setMemory( A.anzmem );
    int32_t error=0;

    SymPoly tmp;
    for(int32_t a=0;a<A.anzterms;a++) {
        tmp.setToZero();
        if (replace_iter_TTerm(tmp,A.terms[a]) != 0) return -1;

        if (erg.addPoly(tmp) != 0) return -1;

    } // a

    return error;
}

int32_t replace_iter_TTerm(
    SymPoly& erg,
    SymTerm& A
) {
    // ret != 0 => error, else success
    // pure factor*c^cexponent
    SymTerm T;
    T.copyFrom(A);
    for(int32_t k=0;k<PERIOD;k++) T.set_bExponent_id_w(k,0);
    if (erg.addTerm(T) != 0) return -1;

    SymPoly tmp,tmp2;
    for(int32_t k=0;k<PERIOD;k++) {
        if (A.bexponent[k] <= 0) continue; // nbothing to do

        if (replace_iter_T_bid_expo(tmp,k,A.bexponent[k]) != 0) return -1;

        tmp2.copyFrom(erg);
        if (polyMul_TAB(erg,tmp2,tmp) != 0) return -1;
    } // k

    return 0;

}

int32_t replace_iter_T_bid_expo(
    SymPoly& erg,
    const int32_t bnr,
    const int32_t bexpo
) {
    // e.g. even exponent b_bnr^bexpo
    // will be replaced by (b_bnr+1 mod p-c)^(bexpo/2)
    // (with odd exponent, b_bnr^1 remains
    // this is DIFFERENT from e.g.
    // (b_bnr+1 mod p)*b_bnr^(bexpo-2), which, if performed
    // to the full extent leads to the same ifnal result
    // butr at one step, it's idfferent
    erg.setToZero();

    // ret != 0 => error, success else
    if (bexpo < 0) {
        printf("\nerror. no netfgauvioe exponent allowed.\n");
        exit(99);
    }

    if (bexpo == 0) {
        erg.setToOne();
        return 0;
    }

    SymTerm T;

    if (bexpo == 1) {
        T.setToZero();
        T.setFactor(1);
        T.set_bExponent_id_w(bnr,1);

        return erg.addTerm(T);
    }

    int32_t mod=bexpo % 2;
    int32_t div=(bexpo-mod) >> 1;
    int32_t bnext=(bnr+1) % PERIOD;

    // term=b_nr^bexpo
    // erg=(b_nr^2)^div * b_nr^mod
    // erg=(b_next-c)^div*b_nr^mod
    T.setToZero();
    T.setFactor(1);
    T.set_bExponent_id_w(bnr,mod);
    erg.setToZero();
    if (erg.addTerm(T) != 0) return -1;
    // erg=(b_next-c)^div*b_nr^mod
    SymPoly bnextc;
    bnextc.setToZero();
    T.setToZero();
    T.setFactor(1);
    T.set_bExponent_id_w(bnext,1);
    if (bnextc.addTerm(T) != 0) return -1;
    T.setToZero();
    T.setFactor(1);
    T.set_cExponent(1);
    if (bnextc.subTerm(T) != 0) return -1;

    SymPoly tmp;
    for(int32_t k=1;k<=div;k++) {
        tmp.copyFrom(erg);
        if (polyMul_TAB(erg,tmp,bnextc) != 0) return -1;
    } // k

    return 0; // success

}

int32_t split_TsymTnonA(
    SymPoly& ergsym,
    SymPoly& ergnon,
    SymPoly& inA
) {
    // ret != 0 => error, else success
    // splits the polynomial A into a symmetric and a non-symmetric
    // part

    ergsym.setToZero();
    ergnon.setToZero();

    // go over A and pick an arbitrary term
    // remove all terms with the same c-exponent and
    // a permutation in the b-exponents
    // CAVE: take cxare of the FACTOR (ggT necessary)

    SymPoly Acurr;
    Acurr.copyFrom(inA);
    SymPoly split;
    split.setMemory( inA.anzmem );

    while (1) {
        if (Acurr.anzterms <= 0) break; // done

        // take the first term
        SymTerm major;
        major.copyFrom(Acurr.terms[0]);
        split.setToZero();
        split.addTerm( major );
        int32_t nze=major.get_nzb();

        // remove out all terms with identical c-.exponent
        // and a permutation of the b-exponents
        for(int32_t k=1;k<Acurr.anzterms;k++) {
            if (Acurr.terms[k].cexponent != major.cexponent) continue;

            if (predPermutation_AB(major.bexponent,Acurr.terms[k].bexponent) <= 0) continue;

            // found one => set to split
            split.addTerm( Acurr.terms[k] );
        } // k

        if (nze == 0) {
            // constant term (pure number or c-term8s) */
            // is symmwtric, but cannot be expressed by elementary
            // symmetric functions
            // => treat as ASYMMETRIC
            if (ergnon.addPoly(split) != 0) {
                return -1;
            }
            if (Acurr.subPoly(split) != 0) {
                return -1;
            }

            continue; // next round
        }

        // now it be symmetrical
        int64_t expected=getbPermutationCount(major);

        if (split.anzterms != expected) {
            // non-symmetrical
            if (ergnon.addPoly(split) != 0) {
                return -1;
            }
            if (Acurr.subPoly(split) != 0) {
                return -1;
            }

            continue; // next round
        }

        // could be symmetriclal, now check the factors
        // if they havbe difering signs => non-symnmetrical
        int8_t signeq=1;
        int8_t sign=split.terms[0].factor.vorz;
        for(int32_t k=1;k<split.anzterms;k++) {
            if (split.terms[k].factor.vorz != sign) {
                signeq=0;
                break;
            }
        } // k

        if (signeq <= 0) {
            // non-symmetrical

            if (ergnon.addPoly(split) != 0) {
                return -1;
            }
            if (Acurr.subPoly(split) != 0) {
                return -1;
            }

            continue; // next round
        }


        // if sign < 0 => setf=max of all factors
        // else min
        BigInt setf;
        setf.copyFrom(split.terms[0].factor);
        for(int32_t k=1;k<split.anzterms;k++) {
            if (sign > 0) {
                if (bigintVgl_AB(split.terms[k].factor,setf) < 0) {
                    setf.copyFrom(split.terms[k].factor);
                }
            } else {
                if (bigintVgl_AB(split.terms[k].factor,setf) > 0) {
                    setf.copyFrom(split.terms[k].factor);
                }
            }
        } // k

        // now set setf to all factors in split
        // this part is symmetrical
        for(int32_t k=0;k<split.anzterms;k++) {
            split.terms[k].setFactor(setf);
        } // k

        if (ergsym.addPoly(split) != 0) {
            printf("\nsplit: "); split.ausgabe(stdout,-1);
            printf("\nsetf: "); setf.ausgabe(stdout);
            return -1;
        }
        if (Acurr.subPoly(split) != 0) {
            return -1;
        }

    } // while

    return 0; // success

}

int8_t predPermutation_AB(
    int32_t* A,
    int32_t* B
) {
    // are the entries of A a permutation of B ?=
    // ret > 0 => yes, else no
    // A,B are of size PERIOD
    int8_t used[ANZSYMVARS];
    for(int32_t k=0;k<PERIOD;k++) used[k]=0;

    for(int32_t a=0;a<PERIOD;a++) {
        // look for A.bexponent[a] in the unused-parts of B.bexponent

        int32_t found=-1;
        for(int32_t b=0;b<PERIOD;b++) {
            if (used[b] > 0) continue;

            if (B[b] == A[a]) {
                found=b;
                break;
            }
        } // b

        if (found < 0) return 0; // no permutation
        used[found]=1;

    } // a

    return 1; // permutation

}

int64_t fakultaet(const int64_t w) {
    if (w > 12) {
        printf("\nn! not supported that high\n");
        exit(99);
    }

    int64_t erg=1;

    for(int64_t k=2;k<=w;k++) erg *= k;

    return erg;
}

int64_t getbPermutationCount(SymTerm& A) {
    // b_i^3 => 5
    // b_i^3*b_k => 20
    // b_i^3*b_k*b_l => 30 ? as multiplication is symmetric
    int64_t nze=A.get_nzb();
    int64_t num = 1;
    int64_t w0=PERIOD;
    for(int32_t k=1;k<=nze;k++) {
        num = num * w0;
        w0--;
    }

    // let {a,b,c..} be the set of DIFFERENT exponents
    // and N(a), N(b) etc their occurance number
    // then erg = erg / ( N(A)! * N(B)! * ...
    int32_t what[ANZSYMVARS];
    int32_t howoften[ANZSYMVARS];
    int32_t anz=0;

    for(int32_t k=0;k<PERIOD;k++) {
        if (A.bexponent[k] <= 0) continue;

        int32_t idx=-1;
        for(int32_t m=0;m<anz;m++) {
            if (what[m] == A.bexponent[k]) {
                idx=m;
                break;
            }
        }

        if (idx >= 0) {
            howoften[idx]++;
        } else {
            howoften[anz]=1;
            what[anz]=A.bexponent[k];
            anz++;
        }
    } // k

    int64_t div=1;
    for(int32_t k=0;k<anz;k++) {
        div = div * fakultaet(howoften[k]);
    }

    if ( (num % div) != 0) {
        printf("\nerror. expected formula wrong\n");
        exit(99);
    }

    return (num / div); // exact
}

int8_t AZterm::save(FILE* f) {
    if (!f) return 0;

    factor.save(f);
    WRITE32(f,cexponent)
    WRITE32(f,ANZSYMVARS)
    for(int32_t k=0;k<ANZSYMVARS;k++) {
        WRITE32(f,azexponent[k])
    }

    return 1;

}

int8_t AZterm::load(FILE* f) {
    if (!f) return 0;

    factor.load(f);
    READ32(f,cexponent)
    int32_t sym;
    READ32(f,sym)
    if (sym != ANZSYMVARS) {
        printf("\nerror. file stored different symbol number\n ");
        exit(99);
    }
    for(int32_t k=0;k<ANZSYMVARS;k++) {
        READ32(f,azexponent[k])
    }

    return 1;

}

int8_t SymTerm::save(FILE* f) {
    if (!f) return 0;

    factor.save(f);
    WRITE32(f,cexponent)
    WRITE32(f,PERIOD)
    for(int32_t k=0;k<PERIOD;k++) {
        WRITE32(f,bexponent[k])
    }

    return 1;

}

int8_t SymTerm::load(FILE* f) {
    if (!f) return 0;

    factor.load(f);
    READ32(f,cexponent)
    int32_t sym;
    READ32(f,sym)
    if (sym != PERIOD) {
        printf("\nerror. file stored different root number\n ");
        exit(99);
    }
    for(int32_t k=0;k<PERIOD;k++) {
        READ32(f,bexponent[k])
    }

    return 1;

}

void AZterm::invertSign(void) {
    factor.vorz *= -1;
}

void AZterm::setToZero(void) {
    factor.setToZero();
    cexponent=0;
    for(int32_t k=0;k<ANZSYMVARS;k++) azexponent[k]=0;
}

void AZterm::setFactor(BigInt& A) {
    factor.copyFrom(A);
}

int32_t AZterm::setFactor(const int64_t a) {
    char tt[256];
    sprintf(tt,"%I64d",a);
    return factor.setstr(tt);
}

int32_t AZterm::setFactor(const char* as) {
    return factor.setstr(as);
}

void AZterm::set_azExponent_id_w(const int32_t nr,const int32_t w) {
    if (
        (nr < 0) || (nr >= ANZSYMVARS)
    ) {
        fprintf(stderr,"error. set azterm out of index %i [0..%i]\n",nr,ANZSYMVARS-1);
        exit(99);
    }

    if (
        (w < 0)
    ) {
        fprintf(stderr,"error. negative az-exponent.\n");
        exit(99);
    }

    azexponent[nr]=w;

}

void AZterm::set_azExponent_ch_w(const char ch,const int32_t w) {
    int nr=-1;
    if ( (ch >= 'A') && (ch <= 'Z') ) nr=ch-'A';
    else
    if ( (ch >= 'a') && (ch <= 'z') ) nr=ch-'a';

    set_azExponent_id_w(nr,w);
}

void AZterm::set_cExponent(const int32_t w) {
    if (w < 0) {
        printf("\nerror. negative c-exponent not allowed.\n");
        exit(99);
    }
    cexponent=w;
}

void AZterm::ausgabe(FILE* f) {
    factor.ausgabe(f);
    if (cexponent > 0) {
        fprintf(f,"*c");
        if (cexponent > 1) fprintf(f,"^%i",cexponent);
    }

    for(int32_t k=0;k<ANZSYMVARS;k++) {
        if (azexponent[k] > 0) {
            fprintf(f,"*%c",elemsymname[k]);
            if (azexponent[k] > 1) fprintf(f,"^%i",azexponent[k]);
        }
    }

}

void AZterm::copyFrom(AZterm& A) {
    factor.copyFrom(A.factor);
    cexponent=A.cexponent;
    for(int32_t k=0;k<ANZSYMVARS;k++) azexponent[k]=A.azexponent[k];
}

int32_t expressByElem_Texp_Tnon_A(
    AZPoly& erg,
    SymPoly& non,

    SymPoly& inpoly
) {
	erg.setToZero();
	non.setToZero();

    int32_t error=0;

	// f3curr:expand(f3f),
	SymPoly curr;
	error += split_TsymTnonA(curr,non,inpoly);

	if (error != 0) return -1;
	if (curr.anzterms <= 0) {
        erg.setToZero();
        return 0; // success
	}

	// f3ret:0,
	erg.setToZero();

	SymTerm lt,sterm;
	SymPoly sprunning,sphelp,spcopy;
	AZPoly aprunning,aphelp;
	AZterm aterm;

	// while f3curr # 0 and f3goon # 0 do (
	while (curr.anzterms > 0) {
        error += polyLtbexponent_TA(lt,curr);
        if (error != 0) return -1;

        sterm.setToZero();
        sterm.setFactor(lt.factor);
        sterm.set_cExponent(lt.cexponent);

        sprunning.setToZero();
        error += sprunning.addTerm(sterm);

        aterm.setToZero();
        aterm.setFactor(lt.factor);
        aterm.set_cExponent(lt.cexponent);
        aprunning.setToZero();
        error += aprunning.addTerm(aterm);

        for(int32_t i=0;i<PERIOD;i++) {
            int32_t diffexpo=0;
            if (i < (PERIOD-1)) diffexpo=lt.bexponent[i] - lt.bexponent[i+1];
            else diffexpo=lt.bexponent[i];

            if (diffexpo != 0) {
                aterm.setToZero();
                aterm.setFactor(1);
                aterm.set_azExponent_id_w(i,diffexpo);
                aphelp.copyFrom(aprunning);
                error += azpolyMul_TTermA(aprunning,aterm,aphelp);

                int32_t e1=error;
                error += polyPow_TAE(sphelp,sigma[i],diffexpo);
                spcopy.copyFrom(sprunning);
                error += polyMul_TAB(sprunning,sphelp,spcopy);

                if (error != 0) {
                    return -1;
                }

            }
        }

        spcopy.copyFrom(curr);
        error += polySub_TAB(curr,spcopy,sprunning);

        aphelp.copyFrom(erg);
        error += azpolyAdd_TAB(erg,aphelp,aprunning);

	}// while

	if (error != 0) return -1;

	// if the non-symmetric part contains e.g. only c-terms but no roots
	// => it will be added to the wsymmetric part as it is
	if (
        (non.anzterms > 0) &&
        (non.get_maxBexponent() <= 0)
	) {
        AZPoly azh;
        if (convertSymToAZ_TA(azh,non) != 0) return -1;
        if (erg.addPoly(azh) != 0) return -1;
        non.setToZero(); // non-symmetric part now empty
	}

    return error;
}

int32_t polyLtbexponent_TA(
    SymTerm& erg,
    SymPoly& A
) {
    // ret <= 0 => error, > 0 => success
    // leading term in b-Exponents
    // ignore c-exponent,as thoise will be deltat with by
    // repeated call
    if (A.anzterms <= 0) {
        erg.setToZero();
        return 0; // success
    }

    erg.copyFrom( A.terms[0] );
    for(int32_t a=1;a<A.anzterms;a++) {
        for(int32_t b=0;b<PERIOD;b++) {
            if (A.terms[a].bexponent[b] > erg.bexponent[b]) {
                // new better found
                erg.copyFrom(A.terms[a]);
                break;
            } else if (A.terms[a].bexponent[b] < erg.bexponent[b]) {
                break;
            }
        } // b
    } // a


    return 0;

}

void initSigma(const int32_t APER) {
    // constructs the elementary sym functions in the roots b0,b1...
    // sigma[0] = sum roots (sigma1, represents A) etc

    int32_t error=0;
    SymTerm T;
    T.setToZero();
    T.setFactor(1);
    // SymPoly
    for(int32_t k=0;k<APER;k++) {
        T.setToZero();
        T.setFactor(1);
        for(int32_t i=0;i<=k;i++) {
            T.set_bExponent_id_w(i,1);
        }

        error += createFullySymmetric_TA(sigma[k],T);
        // as only the elementary part is needed, allf actors
        // can be set to 1
        for(int32_t i=0;i<sigma[k].anzterms;i++) {
            sigma[k].terms[i].setFactor(1);
        }

    } // k

    if (error != 0) {
        printf("\nerror. cannot build sigma, powwersum\n");
        exit(99);
    }

}

int32_t polyPow_TAE(
    SymPoly& erg,
    SymPoly& A,
    const int32_t expo
) {
    // ret != 0 => error, success else
    int32_t error=0;

    if (expo < 0) return -1; // pure polynomials only

    // continued multiplication
    if (expo == 0) {
        erg.setToOne();
        return 0; // success
    }

    SymPoly help;
    erg.copyFrom(A);

    for(int32_t k=2;k<=expo;k++) {
        // erg <- erg*A
        help.copyFrom(erg);
        error += polyMul_TAB(erg,help,A);
    } // k

    return error;

}

int32_t convertSymToAZ_TA(
    AZPoly& ergaz,
    SymPoly& insym
) {
    // ret != 0 => error
    // if insym is free of b's => can be expressed as an AZpoly
    int32_t error=0;

    if (insym.get_maxBexponent() > 0) return -1; // error

    ergaz.setToZero();
    AZterm azt;
    for(int32_t k=0;k<insym.anzterms;k++) {
        azt.setToZero();
        azt.setFactor(insym.terms[k].factor);
        azt.set_cExponent(insym.terms[k].cexponent);
        error += ergaz.addTerm(azt);
    } // k

    return 0;
}

int32_t cyclicSum_TA(
    SymPoly& erg,
    SymTerm& inA
) {
    // constructs the cyclic sum of A where the indicies are SHIFTED to the
    // right by one for PERIOD times
    erg.setToZero();

    int32_t error=0;
    SymTerm A;
    A.copyFrom(inA);

    for(int32_t k=0;k<PERIOD;k++) {
        error += erg.addTerm(A);

        // now shift exponents one to the right
        int32_t sic=A.bexponent[0 ];

        for(int32_t m=1;m<PERIOD;m++) {
            A.bexponent[m-1]=A.bexponent[m];
        } // m

        A.bexponent[PERIOD-1]=sic;

    } // k

    return error;
}

int32_t Equation0::fullSimplify(void) {
    if (AUTOMATICSIMPLIFY <= 0) return 0; // no error

    int32_t error=0;

    // reduce degree of roots ojn right side
    // express right side by elem sym funcs and move to the left
    // move pure c-terms to the left
    // reduce gcd, variabvles

    if (REPEATITER <= 0) {
        if (nonright.anzterms > 0) {
            if (REPEATITER <= 0) {
                error += replace_repeated_iter(nonright);
            }
        }
    }

    if (nonright.anzterms > 0) {
        error += expressElemRight();
    }

    if (nonright.anzterms > 0) {
        error += shiftPureCtermsToLeft();
    }

    error += reduce_gcd_c();

    return error;

}

int8_t Equation0::load(FILE* f) {
    if (!f) return 0;

    int32_t w32;
    READ32(f,w32)
    id=w32;

    if (azleft.load(f) <= 0) return 0;
    if (nonright.load(f) <= 0) return 0;

    return 1;

}

int8_t Equation0::save(FILE* f) {
    if (!f) return 0;

    WRITE32(f,id)
    if (azleft.save(f) <= 0) return 0;
    if (nonright.save(f) <= 0) return 0;

    return 1;

}

int32_t Equation0::expressElemRight(void) {
    // tries to express the NONRIGHT part in elem sym func and move those
    // to the left

    int32_t error=0;
    AZPoly az;
    SymPoly non;

    error += expressByElem_Texp_Tnon_A(az,non,nonright);
    if (error != 0) return -1;
    error += azleft.subPoly(az);
    nonright.copyFrom(non);
    error += shiftPureCtermsToLeft();

    return error;

}

int32_t Equation0::solveForLinear_TI(
    AZRational& solve,
    const int32_t VARIDX
) {
    // ret < 0 => error, =0 => success
    // right-side must be zero
    if (nonright.anzterms > 0) return -1; // error
    if ( (VARIDX < 0) || (VARIDX >= ANZSYMVARS) )  return -1;

    int32_t hid=azleft.get_hipow(VARIDX);
    if (hid != 1) return -1; // onnly linear can be solved

    AZPoly coeff;
    if (azleft.coeff_TID(coeff,VARIDX,1) != 0) return -1; // error

    // now copy all terms that do not have VARIDX to ZAEHLER
    solve.zaehler.setToZero();
    int32_t error=0;
    for(int32_t k=0;k<azleft.anzterms;k++) {
        if (azleft.terms[k].azexponent[VARIDX] == 0) {
            error += solve.zaehler.subTerm( azleft.terms[k] );
        }
    } // k

    solve.nenner.copyFrom( coeff );

    return error;

}

char* Equation0::getVariablesLeft(char* erg) {
    // erg must provide enough emmeory for ANZSYMVARS + 1 (c)
    for(int32_t v=0;v<ANZSYMVARS;v++) erg[v]='-';
    erg[ANZSYMVARS]='-'; // c
    erg[ANZSYMVARS+1]=0;

    for(int32_t k=0;k<azleft.anzterms;k++) {
        for(int32_t v=0;v<ANZSYMVARS;v++) {
            if (azleft.terms[k].azexponent[v] != 0) erg[v]=elemsymname[v];
        } // v
        if (azleft.terms[k].cexponent != 0) erg[ANZSYMVARS]='c';
    } // k

    return erg;
}

int32_t Equation0::reduce_gcd_c(void) {
    if (azleft.anzterms <= 0) {
        return 0; // nothing to do
    }

    int32_t error=0;

    struct Lokal {
        BigInt gcd,t1,rem;
        int32_t lowest_cExponent;
    };
    Lokal *var=new Lokal;
    if (!var) return -1; // error
    int8_t first;

    #define RETURNRED(WW) \
    {\
        if (var) delete var;\
        return (WW);\
    }

    #define COLLECTGCD(PP) \
    {\
        for(int32_t k=0;k<(PP).anzterms;k++) {\
            if ((PP).terms[k].cexponent < var->lowest_cExponent) {\
                var->lowest_cExponent=(PP).terms[k].cexponent;\
            }\
            \
            var->t1.copyFrom(var->gcd);\
            /* sign of factor will be ignored by gcd algorithm */\
            error += bigintGcd_abs_TAB(\
                var->gcd,\
                var->t1,\
                (PP).terms[k].factor\
            );\
        } /* k */\
    }

    #define DIVIDEBYGCD(PP) \
    {\
        for(int32_t k=0;k<(PP).anzterms;k++) {\
            error += bigintDiv_TRAB(\
                var->t1,\
                var->rem,\
                (PP).terms[k].factor,\
                var->gcd\
            );\
            if (var->rem.isZero() == 0) error-=1;\
            if (error != 0) {\
                printf("\nerror bigint div\n");\
                printf("\n  division not remainder free or overflow\n");\
                printf("  "); (PP).terms[k].factor.ausgabe(stdout);\
                printf("\n  /\n  ");\
                var->gcd.ausgabe(stdout);\
                printf("\n  div= "); var->t1.ausgabe(stdout);\
                printf("\n  rem(z%i)= ",var->rem.isZero()); var->rem.ausgabe(stdout);\
                RETURNRED(error)\
            }\
            \
            (PP).terms[k].setFactor( var->t1 );\
            int32_t nc=(PP).terms[k].cexponent - var->lowest_cExponent;\
            if (nc < 0) {\
                printf("\nerror. c exponent: %i\n",nc);\
                exit(99);\
            }\
            (PP).terms[k].set_cExponent( nc );\
        } /* k */\
    }

    // case 1: nonright = 0
    //      go over all terms in azleft and find the
    //      gcd (of |factor|) and the lowest cExponent of all
    //      if gcd > 1 => divide all terms by the (positive) value
    //      if lowest cExponent > 0 => decrease cExponent in all terms
    //      (basically a division by c)

    if (nonright.anzterms <= 0) {
        // initialize
        var->gcd.copyFrom( azleft.terms[0].factor );
        var->gcd.absolute();
        var->lowest_cExponent=azleft.terms[0].cexponent;

        COLLECTGCD(azleft)
        var->gcd.absolute();

        // dividing ?
        if (
            (var->lowest_cExponent > 0) ||
            (bigintVgl_AB(var->gcd,consts->big1) > 0)
        ) {
            DIVIDEBYGCD(azleft)
        }

        RETURNRED(error)

    } // empty rhs

    // case 2: nonright != 0
    //      gcd of |factor| of both the left and right side
    //      lowest cExponent
    //      and accordinglyx divide both sides and/or substract cExponent
    if (nonright.anzterms > 0) {

        var->gcd.copyFrom( nonright.terms[0].factor );
        var->gcd.absolute();
        var->lowest_cExponent=nonright.terms[0].cexponent;

        COLLECTGCD(azleft)
        COLLECTGCD(nonright)
        var->gcd.absolute();

        if (
            (var->lowest_cExponent > 0) ||
            (bigintVgl_AB(var->gcd,consts->big1) > 0)
        ) {
            DIVIDEBYGCD(azleft)
            DIVIDEBYGCD(nonright)
        }

        RETURNRED(error)

    }

    RETURNRED(error)

}

int32_t Equation0::shiftPureCtermsToLeft(void) {
    if (AUTOMATICSIMPLIFY <= 0) return 0; // no error

    if (nonright.anzterms <= 0) return 0; // no error, but nothing to do

    int32_t error=0;
    SymPoly spure;
    spure.setToZero();

    for(int32_t k=0;k<nonright.anzterms;k++) {
        if (nonright.terms[k].get_nzb() <= 0) {
            error += spure.addTerm( nonright.terms[k] );
            if (error != 0) return -1;
        }
    } // k

    AZPoly azpure;
    error += convertSymToAZ_TA(azpure,spure);
    if (error != 0) return -1;

    error += nonright.subPoly( spure );
    error += azleft.subPoly( azpure );

    return error;

}

void Equation0::copyFrom(Equation0& A) {
    id=A.id;
    azleft.copyFrom(A.azleft);
    nonright.copyFrom(A.nonright);
}

Equation0::Equation0() {
    id=-1;
    azleft.setToZero();
    nonright.setToZero();
}

void Equation0::ausgabe(FILE* f,const int32_t maxt) {
    fprintf(f,"[id=%d] ",id);
    azleft.ausgabe(f,maxt);
    fprintf(f," = ");
    nonright.ausgabe(f,maxt);
}

void Equation0::setAll_ABC(
    const int32_t aid,
    AZPoly& inaz,
    SymPoly& insym
) {
    id=aid;
    azleft.copyFrom(inaz);
    nonright.copyFrom(insym);

}

int8_t EquationList::save(FILE* f) {
    // ret <= 0 => error, > 0 => correct
    if (!f) return 0; // error

    WRITE32(f,anzeq)
    for(int32_t k=0;k<anzeq;k++) {
        if (aseq[k].save(f) <= 0) return 0;
    } // k

    return 1;
}

int8_t EquationList::load(FILE* f) {
    // ret <= 0 => error, > 0 => correct
    if (!f) return 0; // error

    fastClear();

    int32_t w32;
    READ32(f,w32)
    Equation0 eq;
    for(int32_t k=0;k<w32;k++) {
        if (eq.load(f) <= 0) return 0;
        if (eq.fullSimplify() != 0) return 0;
        if (addEquation(eq) != 0) return 0;
    } // k

    return 1;
}

Equation0* EquationList::getEquationPtr_id(const int32_t aid) {
    for(int32_t k=0;k<anzeq;k++) {
        if (aseq[k].id == aid) return &aseq[k];
    } // k

    return NULL;
}

int32_t EquationList::makeAddEquation_fullSymmetric_IFCBs(
    const int32_t aid,
    const int64_t afactor,
    const int32_t acexponent,
    const char* bstr
) {
    SymTerm T;
    T.setToZero();
    T.setFactor( afactor );
    T.set_cExponent( acexponent );
    char* all=new char[strlen(bstr)+4];
    char *all0=all;
    sprintf(all,"%s,",bstr);

    if (fdescr) {
        char tt[2048];
        tt[0]=0;
        if (afactor != 1) sprintf(&tt[strlen(tt)],"%i*",afactor);
        if (acexponent == 1) sprintf(&tt[strlen(tt)],"c*");
        else if (acexponent > 1) sprintf(&tt[strlen(tt)],"c^%i*",acexponent);

        fprintf(fdescr,"\n  %ssum [ ",tt);
    }

    int32_t error=0,bctr=0;
    for(int32_t k=0;k<PERIOD;k++) {
        char *p=strchr(all,',');
        if (!p) {
            // too few
            // set rest of exponents to 0
            for(int32_t k2=k;k2<PERIOD;k2++) T.set_bExponent_id_w(k2,0);
            break; // too few
        }
        p[0]=0;
        int64_t w;
        if (sscanf(all,"%i",&w) != 1) { error=-1; break; }
        T.set_bExponent_id_w(k,w);
        if (fdescr) {
            if (w > 0) {
                if (bctr > 0) fprintf(fdescr,"*");
                fprintf(fdescr,"b_i%i",bctr);
                if (w > 1) fprintf(fdescr,"^%i",w);
                bctr++;
            }
        }
        all=p+1; // next

    } // k

    delete[] all0;

    if (error != 0) return -1;

    if (fdescr) {
        fprintf(fdescr,"]\n");
        fprintf(fdescr,"    for indices ");
        fprintf(fdescr,"0 <= ");
        for(int32_t k=0;k<bctr;k++) {
            if (k > 0) fprintf(fdescr," <> ");
            fprintf(fdescr,"i%i",k);
        }
        fprintf(fdescr," <= %i\n",PERIOD-1);
    }

    SymPoly sp;
    error += createFullySymmetric_TA(sp,T);
    if (error != 0) return -1;

    Equation0 eq0;
    error += makeEquationSymPoly_TAI(eq0,sp,aid);
    if (error != 0) return -1;

    error += eq0.fullSimplify();
    error += addEquation( eq0 );

    return error;

}

int32_t EquationList::copyFrom(EquationList& A) {
    int32_t error=0;
    fastClear();

    for(int32_t k=0;k<A.anzeq;k++) {
        error += addEquation(A.aseq[k]);
    } // k

    return error;
}

int32_t EquationList::get_VariablesLeft(char* erg) {
    // erg must provide enough space for
    // ANZSYMVARS plus 1 (small c) and plus 1 (EOL)
    erg[0]=0;

    int32_t anz=0;
    int8_t present[ANZSYMVARS];
    int8_t presentc=0;
    int32_t exc=ANZSYMVARS;
    for(int32_t k=0;k<ANZSYMVARS;k++) {
        erg[k]='-';
        present[k]=0;
    }
    erg[exc]='-'; // for small c
    erg[ANZSYMVARS+1]=0;

    for(int32_t eq=0;eq<anzeq;eq++) {
        for(int32_t azt=0;azt<aseq[eq].azleft.anzterms;azt++) {
            for(int32_t ex=0;ex<ANZSYMVARS;ex++) {
                if (aseq[eq].azleft.terms[azt].azexponent[ex] != 0) {
                    if (present[ex] <= 0) {
                        present[ex]=1;
                        erg[ex]='A'+ex;
                        anz++;
                    }
                }

                if (aseq[eq].azleft.terms[azt].cexponent != 0) {
                    if (presentc <= 0) {
                        anz++;
                        presentc=1;
                        erg[exc]='c';
                    }
                }
            } // ex
        } // azt
    } // eq

    return anz;
}


int32_t EquationList::removeEquationById(const int32_t aid) {
    int32_t pos=-1;
    for(int32_t k=0;k<anzeq;k++) {
        if (aseq[k].id == aid) { pos=k; break; }
    }
    if (pos < 0) return -1; // treat as error

    if (anzeq == 1) {
        fastClear();
        return 0; // success
    }

    if (pos != (anzeq-1)) {
        aseq[pos].copyFrom(aseq[anzeq-1]);
    }

    anzeq--;

    return 0; // success
}

EquationList::EquationList() {
    anzeq=0;
}

void EquationList::fastClear(void) {
    anzeq=0;
}

int32_t EquationList::addEquation(Equation0& A) {
    return addEquation(
        A.id,
        A.azleft,
        A.nonright
    );
}

int32_t EquationList::addEquation(
    const int32_t aid,
    AZPoly& A,
    SymPoly& B
) {
    // ret < 0 => error
    // =0 => success

    if (anzeq >= (MAXEQ-2)) {
        printf("\nerror. equationlist not large enough\n");
        exit(99);
        return -1;
    }

    // id must be unieuq
    if (getEquationPtr_id(aid)) return -1; // already exists

    aseq[anzeq].id=aid;
    aseq[anzeq].azleft.copyFrom(A);
    aseq[anzeq].nonright.copyFrom(B);

    anzeq++;

    return 0; // success
}

void EquationList::ausgabe(FILE* f,const int32_t maxt) {
    char vars[ANZSYMVARS+8]; // PERIOD + 2 (1 for c, 1 for EOL)

    int32_t anzv=get_VariablesLeft(vars);

    fprintf(f,"\n|- %i equations in %i variables: %s\n\n",
        anzeq,anzv,vars

    );
    for(int32_t k=0;k<anzeq;k++) {
        aseq[k].ausgabe(f,maxt);
        fprintf(f,"\n");
    }
}

void EquationList::ausgabeVars(FILE* f,const int32_t maxt) {
    char vars[ANZSYMVARS+8]; // PERIOD + 2 (1 for c, 1 for EOL)

    int32_t anzv=get_VariablesLeft(vars);

    fprintf(f,"\n|- %i equations in %i variables: %s\n\n",
        anzeq,anzv,vars
    );

    for(int32_t k=0;k<anzeq;k++) {
        aseq[k].getVariablesLeft(vars);
        fprintf(f," %s ",vars);
        int32_t hiidx=-1;
        for(int32_t p=0;p<strlen(vars);p++) {
            // not c
            if (vars[p] == 'c') continue;
            if (vars[p] != '-') hiidx=p;
        }

        if (hiidx >= 0) {
            int32_t hi=aseq[k].azleft.get_hipow(hiidx);
            fprintf(f,"_%c^%i_ ",
                vars[hiidx],
                hi);
        }

        aseq[k].ausgabe(f,maxt);
        fprintf(f,"\n");
    }
}


int8_t termVgl_only_bExponent_AB(
    SymTerm& A,
    SymTerm& B
) {
    // ret < 0 => A < B
    // = 0 => A=B
    // > 0 => A > B
    for(int32_t k=0;k<PERIOD;k++) {
        if (A.bexponent[k] < B.bexponent[k]) return -1;
        if (A.bexponent[k] > B.bexponent[k]) return +1;
    } // k

    return 0;

}

int32_t makeEquationSymPoly_TAI(
    Equation0& erg,
    SymPoly& A,
    const int32_t aid
) {
    // ret: < 0 => error
    // = 0 => success

    // A is a function (need not be fully symmetric)
    // 1) exrpess A as elem sym B + nonB
    // 2) replace_repeated A and thene xpress elem sym C + nonC
    // move constant (number and/or c) to left side

    int32_t error=0;
    AZPoly azB,azC;
    SymPoly nonB,nonC,Aiter;
    error += expressByElem_Texp_Tnon_A(azB,nonB,A);
    if (error != 0) return -1; // error

    // nonB can be NON-repeatedly and contain higher powers
    error += replace_repeated_iter(nonB);

    if (fdescr) {
        fprintf(fdescr,"\nwhich gives\n\n  ");
        azB.ausgabe(fdescr,-1);
        if (nonB.anzterms > 0) {
            fprintf(fdescr," + ");
            nonB.ausgabe(fdescr,-1);
        }
    }

    // A = azB + nonB

    // now repleaceiter A
    //printf("\n  A= "); A.ausgabe(stdout,-1);
    error += replace_repeated_iter_TA(Aiter,A);
    if (error != 0) return -1; // error
    error += expressByElem_Texp_Tnon_A(azC,nonC,Aiter);
    if (error != 0) return -1; // error

    if (fdescr) {
        fprintf(fdescr,"\n    =\n  ");
        azC.ausgabe(fdescr,-1);
        if (nonC.anzterms > 0) {
            fprintf(fdescr," + ");
            nonC.ausgabe(fdescr,-1);
        }
    }

    // A = azC + nonC
    // Equation: azB + nonB = azC + nonC
    // azB - azC = nonC - nonB
    // the shift constant (number, c-terms) from right to left
    error += azB.subPoly( azC );
    error += nonC.subPoly( nonB );
    if (error != 0) return -1; // error

    erg.setAll_ABC(aid,azB,nonC);

    // shiftPureCtermsToLeft
    error += erg.shiftPureCtermsToLeft();

    return error;

}

int32_t substitute_TAB(
    AZRational& erg,
    AZPoly& into,
    Solution& sol
) {
    erg.setToZero();

    // ret < 0 => error
    // ret = 0 => success
    int32_t error=0;
    AZterm azt;
    azt.setToZero();
    azt.setFactor(1);

    AZPoly one;
    one.setToZero();

    if (one.addTerm(azt) != 0) return -1;

    int32_t hi=into.get_hipow(sol.varidx);
    if (hi <= 0) {
        // variable doesn't occur => result is polynomial
        erg.setToZero();
        erg.nenner.copyFrom( one );
        erg.zaehler.copyFrom( into );

        return 0; // success

    }

    AZPoly t1;
    erg.setToZero();

    // now build eq.zaehler^k, eq.nenner^k up to hipow
    AZPoly *eqzae=new AZPoly[hi+1];
    if (!eqzae) return -1;

    AZPoly *eqnen=new AZPoly[hi+1];
    if (!eqnen) {
        delete[] eqzae;
        return -1;
    }

    #define RETSUBST(WW) \
    {\
        if (eqzae) delete[] eqzae;\
        if (eqnen) delete[] eqnen;\
        \
        return (WW);\
    }

    eqzae[0].copyFrom(one);
    eqnen[0].copyFrom(one);

    error=0;
    for(int32_t k=1;k<=hi;k++) {
        error += azpolyMul_TAB(eqzae[k],eqzae[k-1],sol.eq.zaehler);
        error += azpolyMul_TAB(eqnen[k],eqnen[k-1],sol.eq.nenner);

        if (error != 0) { RETSUBST(-1) }

    } // k

    erg.setToZero();

    // erg = add over all terms in INTO
    // and add polynomial
    // T'=term*zaehler^[var-exponent]*nenner^[hi-varexo]
    AZPoly t2;
    error=0;
    erg.zaehler.setToZero();
    erg.nenner.copyFrom( eqnen[hi] );

    for(int32_t k=0;k<into.anzterms;k++) {
        // T'=term*zaehler^[var-exponent]*nenner^[hi-varexo]
        azt.copyFrom( into.terms[k] );
        int32_t currentexp=azt.azexponent[sol.varidx];

        // remove VAR-exponent
        azt.set_azExponent_id_w( sol.varidx,0 );

        // T'=term*zaehler^[var-exponent]*nenner^[hi-varexo]
        error += azpolyMul_TTermA(t1,azt,eqzae[currentexp]);
        // T'=t1*nenner^[hi-varexo]
        error += azpolyMul_TAB(t2,t1,eqnen[hi-currentexp]);
        // T'=t2
        if (error != 0) { RETSUBST(-1) }
        error += erg.zaehler.addPoly( t2 );

    } // k

    RETSUBST(error)

}

int8_t Solution::load(FILE* f) {
    if (!f) return 0;

    int32_t w32;
    READ32(f,w32)
    solved=w32;
    READ32(f,w32)
    varidx=w32;

    if (eq.load(f) <= 0) return 0;

    return 1;

}

int8_t Solution::save(FILE* f) {
    if (!f) return 0;

    WRITE32(f,solved)
    WRITE32(f,varidx);
    if (eq.save(f) <= 0) return 0;

    return 1;

}

void Solution::ausgabe(FILE* f,const int32_t maxt) {
    if (solved <= 0) return;

    fprintf(f,"%c = ",elemsymname[varidx]);
    eq.ausgabe(f,maxt);

}

Solution::Solution() {
    solved=0;
    varidx=-1;
    eq.setToZero();
}

void Solution::setToUnsolved(void) {
    solved=0;
    varidx=-1;
    eq.setToZero();
}

void Solution::setSolution_id(const int32_t VAR,AZRational& sol) {
    if ( (VAR < 0) || (VAR >= ANZSYMVARS) ) {
        printf("\nerror. var index too out of range at solution\n");
        exit(99);
    }

    varidx=VAR;
    solved=1;
    eq.copyFrom(sol);
}

int32_t eliminate_TTTIEAB(
    AZPoly& erg,
    AZPoly& coeffA,
    AZPoly& coeffB,
    const int32_t varidx,
    const int32_t varexp,
    AZPoly& A,
    AZPoly& B
) {
    // ret < 0 => error, = 0 => success
    // computes: coeff(A,varidx^vardxp)*B - coeff(B,variedx^varexp)*A
    erg.setToZero();
    AZPoly t1;
    int32_t error=0;
    error += A.coeff_TID(coeffA,varidx,varexp);
    error += B.coeff_TID(coeffB,varidx,varexp);
    if (error != 0) return -1;

    // computes: coeff(A,varidx^vardxp)*B - coeff(B,variedx^varexp)*A
    erg.setToZero();
    error += azpolyMul_TAB(erg,coeffA,B);
    // computes: erg - coeff(B,variedx^varexp)*A
    error += azpolyMul_TAB(t1,coeffB,A);
    error += erg.subPoly(t1);

    return error;

}

int32_t azpolyPow_TNA(
    AZPoly& erg,
    const int32_t ex,
    AZPoly& A
) {

    if (A.anzterms <= 0) {
        erg.setToZero();
        return -1; // error
    }

    if (ex < 0) {
        printf("\nerror. negative exponent not allowed.\n");
        exit(99);
    }

    // return as always < 0 => error
    // erg = A^ex
    erg.setToZero();

    if (ex == 0) {
        erg.setToOne();
        return 0;
    }

    if (ex == 1) {
        erg.copyFrom(A);
        return 0;
    }

    int32_t error=0;

    // version using squaring
    AZPoly* arr=new AZPoly[ex+1];
    if (!arr) { printf("\nerror. memory. pow\n"); exit(99); }
    for(int32_t k=0;k<=ex;k++) arr[k].setToZero();

    arr[1].copyFrom( A );
    int32_t hiex=-1;
    for(int32_t k=2;k<=ex;k<<=1) {
        error += azpolySqr_TA(arr[k],arr[k>>1]);
        hiex=k;
    } // k

    erg.copyFrom( arr[hiex] );
    AZPoly t1;

    int32_t todo=ex-hiex;
    while (todo > 0) {
        int32_t p=-1;
        for(int32_t k=1;k<=todo;k<<=1) {
            p=k;
        }

        //printf("\n  todo%i, binaer %i",todo,p);

        if(p <= 0) {
            printf("implementation error pow_tree.\n");
            exit(99);
        }

        if (arr[p].anzterms <= 0) {
            printf("\nerror. pow_tree. not-set polynomial at %i.\n",p);
            printf("\n  for ");
            A.ausgabe(stdout,10);
            printf("\n  ^%i\n",ex);
            exit(99);
        }

        t1.copyFrom(erg);
        error += azpolyMul_TAB(erg,t1,arr[p]);
        todo -= p;

    } // todo

    delete[] arr;

    return error;
}

int32_t azpolyDiv_TRABLH(
	AZPoly& dividing,
	AZPoly& remainder,
	AZPoly& A,
	AZPoly& B,
	const int32_t lvaridx, // leading var

	DivHlp& hlp // temporary objects, will be deleted outside
) {
	remainder.copyFrom(B); // standard return value: error

	if (B.predZero() > 0) {
		LOGMSG("\nError. Polynom division by zero 1\n");
		printf("\nA= "); A.ausgabe(stdout,-1);
		printf("\nB= "); B.ausgabe(stdout,-1);

		exit(99);
	}

	if (A.predZero() > 0) {
		// A / B = 0 / B => remainderfree with diviidend=0
		dividing.setToZero();
		remainder.setToZero();
		return 0; // success
	}

	// leadingTerm in B is: B.terms[anzterms-1] of sorted
	// or search for if array access is used

	int32_t ltb;
	if (B.getLexicographicMax_TI(ltb,lvaridx) != 0) return -1; // error

	if (ltb < 0) {
		LOGMSG("\nError. No leading term in B\n");
		exit(99);
	}

	AZterm* LTB=&B.terms[ltb];

	hlp.hlprest.copyFrom(A);
	dividing.setToZero();
	int32_t error=0;

	while (1) {

		if (hlp.hlprest.predZero() > 0) {
			remainder.setToZero();
			// dividing already set
			return error; // succ
		}

		int32_t ltrest;
		// use highest lexicographic term w.r.t to
		// variable to be eliminated, not the internal
		// sorting of the polynomial sorage
		if (hlp.hlprest.getLexicographicMax_TI(ltrest,lvaridx) != 0) return -1; // error

		if (ltrest < 0) {
			LOGMSG("\nError. PolynomDiv rest has no maximal element\n");
			exit(99);
		}

		// are the leading terms dividable
		AZterm *LTR=&hlp.hlprest.terms[ltrest];

		if (LTR->cexponent < LTB->cexponent) {
			remainder.copyFrom(hlp.hlprest);
			return error;
		}

		for(int32_t v=0;v<ANZSYMVARS;v++) {
            if (LTR->azexponent[v] < LTB->azexponent[v]) {
                remainder.copyFrom(hlp.hlprest);
                return error;
            }
		} // v

		// now variable monomial could be divided
		// what about the coefficient

		AZterm tone;

		error += bigintDiv_TRAB(
			tone.factor,
			hlp.irem,
			hlp.hlprest.terms[ltrest].factor,
			B.terms[ltb].factor
		);

		if (hlp.irem.vorz != 0) {
			// no remainder-free division possible
			remainder.copyFrom(hlp.hlprest);
			// dividing is not relevant the
			return error;
		}

		tone.cexponent=LTR->cexponent - LTB->cexponent;
		for(int32_t v=0;v<ANZSYMVARS;v++) {
            tone.azexponent[v]=LTR->azexponent[v] - LTB->azexponent[v];
		}

		error += azpolyMul_TTermA(hlp.hlpt2,tone,B);
		error += hlp.hlprest.subPoly(hlp.hlpt2);
		error += dividing.addTerm(tone);

		if (error != 0) return -1; // error

	} // while

	// no remainder free division possible
	// just return the last computed values

	remainder.copyFrom(hlp.hlprest);

	return error;

}

// symPolyDiv
int32_t polyDiv_TRABLH(
	SymPoly& dividing,
	SymPoly& remainder,
	SymPoly& A,
	SymPoly& B,
	const int32_t lvaridx, // leading var

	DivHlpSym& hlp // temporary objects, will be deleted outside
) {
	remainder.copyFrom(B); // standard return value: error

	if (B.predZero() > 0) {
		LOGMSG("\nError. Polynom division by zero 1\n");
		printf("\nA= "); A.ausgabe(stdout,-1);
		printf("\nB= "); B.ausgabe(stdout,-1);

		exit(99);
	}

	if (A.predZero() > 0) {
		// A / B = 0 / B => remainderfree with diviidend=0
		dividing.setToZero();
		remainder.setToZero();
		return 0; // success
	}

	// leadingTerm in B is: B.terms[anzterms-1] of sorted
	// or search for if array access is used

	int32_t ltb;
	if (B.getLexicographicMax_TI(ltb,lvaridx) != 0) return -1; // error

	if (ltb < 0) {
		LOGMSG("\nError. No leading term in B\n");
		exit(99);
	}

	SymTerm* LTB=&B.terms[ltb];

	hlp.hlprest.copyFrom(A);
	dividing.setToZero();
	int32_t error=0;

	while (1) {

		if (hlp.hlprest.predZero() > 0) {
			remainder.setToZero();
			// dividing already set
			return error; // succ
		}

		int32_t ltrest;
		// use highest lexicographic term w.r.t to
		// variable to be eliminated, not the internal
		// sorting of the polynomial sorage
		if (hlp.hlprest.getLexicographicMax_TI(ltrest,lvaridx) != 0) return -1; // error

		if (ltrest < 0) {
			LOGMSG("\nError. PolynomDiv rest has no maximal element\n");
			exit(99);
		}

		// are the leading terms dividable
		SymTerm *LTR=&hlp.hlprest.terms[ltrest];

		if (LTR->cexponent < LTB->cexponent) {
			remainder.copyFrom(hlp.hlprest);
			return error;
		}

		for(int32_t v=0;v<ANZSYMVARS;v++) {
            if (LTR->bexponent[v] < LTB->bexponent[v]) {
                remainder.copyFrom(hlp.hlprest);
                return error;
            }
		} // v

		// now variable monomial could be divided
		// what about the coefficient

		SymTerm tone;

		error += bigintDiv_TRAB(
			tone.factor,
			hlp.irem,
			hlp.hlprest.terms[ltrest].factor,
			B.terms[ltb].factor
		);

		if (hlp.irem.vorz != 0) {
			// no remainder-free division possible
			remainder.copyFrom(hlp.hlprest);
			// dividing is not relevant the
			return error;
		}

		tone.cexponent=LTR->cexponent - LTB->cexponent;
		for(int32_t v=0;v<ANZSYMVARS;v++) {
            tone.bexponent[v]=LTR->bexponent[v] - LTB->bexponent[v];
		}

		error += polyMul_TTermA(hlp.hlpt2,tone,B);
		error += hlp.hlprest.subPoly(hlp.hlpt2);
		error += dividing.addTerm(tone);

		if (error != 0) return -1; // error

	} // while

	// no remainder free division possible
	// just return the last computed values

	remainder.copyFrom(hlp.hlprest);

	return error;

}


int32_t azpolySqr_TA(
    AZPoly& erg,
    AZPoly& A
) {
    // idea later: split A into two parts B,C and
    // computed (B+C)^2 = B^2+2*B*C+C^2 = B*(B+2C)+C^2
    AZPoly t1;
    t1.copyFrom(A);

    return azpolyMul_TAB(erg,A,t1);

}

int32_t eqmultiply_TIAFC(
    Equation0& erg,
    const int32_t atgtid,
    Equation0& A,
    const int32_t afactor,
    const int32_t ac
) {
    // ret < 0 => error
    // multiplies the equation A on both sides with the
    // term AFACTOR*C^AC

    AZterm azt;
    azt.setToZero();
    azt.setFactor(afactor);
    azt.set_cExponent(ac);
    SymTerm st;
    st.setToZero();
    st.setFactor( azt.factor );
    st.set_cExponent( azt.cexponent );

    int32_t error=0;
    AZPoly az;
    SymPoly sp;
    error += azpolyMul_TTermA(az,azt,A.azleft);
    error += polyMul_TTermA(sp,st,A.nonright);

    erg.setAll_ABC(atgtid,az,sp);

    return error;

}

int32_t eqadd_TIAB(
    Equation0& erg,
    const int32_t tgtid,
    Equation0& A,
    Equation0& B
) {
    int32_t error=0;

    AZPoly az;
    SymPoly sp;

    error += azpolyAdd_TAB(az,A.azleft,B.azleft);
    error += polyAdd_TAB(sp,A.nonright,B.nonright);

    erg.setAll_ABC(tgtid,az,sp);

    return error;

}

int32_t getCorrectDimension(
	const int32_t degA,
	const int32_t degB,
	const int32_t step
) {
	int32_t s=sum_int32t(degA,degB);

	if (step == 1) return s;

	int32_t m=s % step;
	if ( (m == 1) && (s > step) ) return s; // all fits
	s=sum_int32t(s,step-m); // now MULTISTEP | s
	s=sum_int32t(s,1); // now dim % MULTISTEP = 1

	//LOGMSG3("\n  enlargening dimension to %i for stepsize %i\n",s,GLOBALMULTISTEP);

	return s;

}

int8_t setSylvesterMatrix_TABCDVS(
	MatrixPolynom& sylvester,
	AZPoly& A,
	AZPoly& B,
	AZPoly* Acoeff,
	AZPoly* Bcoeff,
	const int32_t VARIDX,
	const int32_t STEP
) {
    int32_t Ah=A.get_hipow(VARIDX);
    int32_t Bh=B.get_hipow(VARIDX);

    if ( (Ah <= 0) || (Bh <= 0) ) return 0; // error

	int32_t setdim=getCorrectDimension(Ah,Bh,BAREISS_STEP);

	LOGMSG3("\nSylvester matrix dimension %i x %i\n",setdim,setdim);
	sylvester.setDimension(setdim);

	int32_t add=setdim - sum_int32t(Ah,Bh);

	sylvester.setConstant0polynom();
	// diagonal elements in 0..add-1 is 1
	for(int32_t w=0;w<add;w++) {
		sylvester.entryYX[w*setdim+w].setToOne();
	} // w

	// matrix: first M rows: coefficients of F shifted
	int32_t x0=add-1;
	for(int32_t y=add;y<(add+Bh);y++) {
		x0++;

		for(int32_t i=0;i<=Ah;i++) {
			// entry address valid as all y,x0+i,dim< = 2^14
			sylvester.entryYX[y*sylvester.dim+(x0+i)].copyFrom(Acoeff[Ah-i]);
		} // i

	} // y

	// second N rows: coefficients of G shifted
	x0=add-1;
	for(int32_t y=(add+Bh);y<setdim;y++) {
		x0++;

		for(int32_t i=0;i<=Bh;i++) {
			sylvester.entryYX[y*sylvester.dim+(x0+i)].copyFrom(Bcoeff[Bh-i]);
		} // i

	} // y

	return 1; // success

}

int8_t setSylvesterMatrix_TABCDV_old(
	MatrixPolynom& sylvester,
	AZPoly& A,
	AZPoly& B,
	AZPoly* Acoeff,
	AZPoly* Bcoeff,
	const int32_t VARIDX
) {
    // ret <= 0 => werror
	// dimension of Sylvester matrix will be adjusted to
	// fit the used MULTISTEP, so that the last computed matrix
	// is dimension-1 and hold sexactly one entry as
	// the resultant

	if ( (VARIDX < 0) || (VARIDX >= ANZSYMVARS) ) return 0;
	int32_t Ah=A.get_hipow(VARIDX);
	int32_t Bh=B.get_hipow(VARIDX);

	if ( (Ah <= 0) || (Bh <= 0) ) return 0;

	int32_t setdim=sum_int32t(Ah,Bh);

	LOGMSG3("\n  Sylvester matrix dimension %i x %i\n",setdim,setdim);
	sylvester.setDimension(setdim);

	if ( (Ah >= MAXSYLDIMENSION) || (Bh >= MAXSYLDIMENSION) ) {
		LOGMSG2("\nError. Polynomial degree must not exceed %i\n",MAXSYLDIMENSION);
		exit(99);
	}

	sylvester.setConstant0polynom();

	// matrix: first M rows: coefficients of F shifted
	int32_t x0=-1;
	for(int32_t y=0;y<Bh;y++) {
		x0++;

		for(int32_t i=0;i<=Ah;i++) {
			// entry address valid as all y,x0+i,dim< = 2^14
			sylvester.entryYX[y*sylvester.dim+(x0+i)].copyFrom(
                Acoeff[Ah-i]);
		} // i

	} // y

	// second N rows: coefficients of G shifted
	x0=0-1;
	// setdim <add+G.Zdegree, so addition is fine
	for(int32_t y=Bh;y<setdim;y++) {
		x0++;

		for(int32_t i=0;i<=Bh;i++) {
			sylvester.entryYX[y*sylvester.dim+(x0+i)].copyFrom(
                Bcoeff[Bh-i]);
		} // i

	} // y

	return 1;

}

char* getMatrixFileName(char* erg,const int32_t k) {
	sprintf(erg,"_%s_k%i.matrix",fnbase,k);

	return erg;
}

// 2-step Bareiss mnethod
// multistep version
// support useMCarray

int32_t fractionFreeGaussMultistep2_TABI(
	AZPoly& res,
	AZPoly& A,
	AZPoly& B,
	const int32_t VARIDX // variable for which to compute resultant wih respect to
) {
	// multistep method by E Bareiss
	res.setToZero();

	/*

		indexing
		M0=Sylvester matrix
		Mk=target matrix to compute (note different from other resultant programs)
		M(dimension): end matrix whose n/n element is the resultant
		rows and columns: 1 to dim
		(access in memory is still C++based on 0)

		1. does file M0 exist ?
			no => create Sylvestermatrix and store
		2. Look for highest value h such that Mh exists
		3. h (1 to dim): h must be divisible by MULTISTEP
		4. k = h+MULTISTEP
		5. Load h and the necessary k',k' value from Mk'
			(k' = divisor,prob. h-2*MULTISTEP ?
		6. Compute matrix k directly onto harddrive

	*/

	int32_t Ah=A.get_hipow(VARIDX);
	int32_t Bh=B.get_hipow(VARIDX);

	if ( (Ah <= 0) || (Bh <= 0) ) return -1; // error

	int32_t expecteddim=getCorrectDimension(Ah,Bh,BAREISS_STEP);

	char fn[2048];
	int32_t klastfullcomputed=-1;
	// delete all existing file names
	char tt[2048];

	//klastfullcomputed=12;

	if (klastfullcomputed < 0) {
        for(int32_t k=expecteddim;k>=0;k--) {
            getMatrixFileName(fn,k);
            sprintf(tt,"del %s 2>nul",fn);
            system(tt);
        } // k
    }

	AZPoly t1,t2,polyConst1,c0,ci1,ci2,rem,nullpoly;
	DivHlp help;

	polyConst1.setToOne();
	nullpoly.setToZero();

	MatrixPolynom* Mkmd=new MatrixPolynom;
	Mkmd->setDimension(expecteddim);
	printf("\ncomputing the resultant with Bareiss' 2-step algorithm ...\n");

	#define RET2(WW) \
	{\
        if (Mkmd) delete Mkmd;\
        return (WW);\
	}

	if (klastfullcomputed < 0) {
		// create Sylvester matrix (with dimension divisible by multistep)
		printf("initial step with Sylvester matrix ... ");

		printf("\ncreating Sylvester matrix ... ");
		printf("\n  getting coefficients ... ");
		int32_t ah,bh;
		AZPoly *Acoeff=A.get_coeffList_TV(ah,VARIDX);
		AZPoly *Bcoeff=B.get_coeffList_TV(bh,VARIDX);

		if ( (ah != Ah) || (bh != Bh) ) {
            printf("\nerror. degress not consistenmt.\n");
            exit(99);
		}

        if (setSylvesterMatrix_TABCDVS(
			*Mkmd,
			A,B,
			Acoeff,Bcoeff,VARIDX,BAREISS_STEP
		) <= 0) {
            printf("\n  error setting Sylvester matirx.\n");
            exit(99);
		}

		char fn2[2048];
		getMatrixFileName(fn2,0);
		Mkmd->save(fn2);

		//Mkmd->ausgabe(stdout);

		klastfullcomputed=0; // M0

		delete[] Acoeff;
		delete[] Bcoeff;
	}

	int32_t ystart=1; // not zero

	// LOOP
	char fn2[2048];
	int32_t error=0;

	for(int32_t kloop=(klastfullcomputed+2);
		kloop < expecteddim;
		kloop += 2
	) {
		printf("\ncomputing %i/%i ... ",kloop,expecteddim-1);
		// free memory

		Mkmd->setDimension(expecteddim);

		getMatrixFileName(fn2,kloop-2);
		Mkmd->load(fn2,-1,-1); // full

		//printf("\n|%s|\n",fn2);
		//Mkmd->ausgabe(stdout);
		//printf("\n");

		AZPoly* adenom=NULL;

		if ( (kloop-2) > 0) {
			adenom=Mkmd->getPointer1based_YX(kloop-2,kloop-2);
		} else {
			adenom=&polyConst1;
		}

		printf(" |%i| ",adenom->anzterms);

		// c0=(	a(k-2)(k-1,k-1)*a(k-2)(k,k) -
		//		a(k-2)(k,k-1)*a(k)(k-1,k)
		//		) / a(k-2)(k-2,k-2)

		int8_t gerror1=0,gerror2=0;
		{
			{
				{
					if (azpolyMul_TAB(
						t1,
						*Mkmd->getPointer1based_YX(kloop-1,kloop-1),
						*Mkmd->getPointer1based_YX(kloop,kloop)
					) != 0) gerror1=1;
				}
				{
					if (azpolyMul_TAB(
						t2,
						*Mkmd->getPointer1based_YX(kloop,kloop-1),
						*Mkmd->getPointer1based_YX(kloop-1,kloop)
					) != 0) gerror2=1;

				}
			}
		}

		if (
            (gerror1 != 0) ||
            (gerror2 != 0)
        ) {
            RET2(-1)
        }

		if (t1.subPoly(t2) != 0) {
            RET2(-1)
        }

		error=0;
		error += azpolyDiv_TRABLH(c0,rem,t1,*adenom,VARIDX,help);
		if ( (error == 0) && (rem.anzterms != 0) ) {
            error--; // error
			LOGMSG("\nError. c0 not fraction free\n");
			printf("\nt1 / ad: |"); t1.ausgabe(stdout,-1);
			printf("|\n / |"); adenom->ausgabe(stdout,-1);
			printf("|\n");
			exit(99);
        }

		if (error != 0) {
            RET2(-1)
        }

		getMatrixFileName(fn2,kloop);
		FILE *fk=fopen(fn2,"wb");
		if (!fk) {
			LOGMSG("\nError. Cannot open target matrix for storage\n");
			exit(99);
		}

		// store dimension
		if (fwrite(&expecteddim,1,sizeof(expecteddim),fk) != sizeof(expecteddim)) {
			LOGMSG("\nError. Writing matrix. Probably invalid file. Deleting recommended.\n");
			exit(99);
		}

		for(int32_t iloop=ystart;iloop<=expecteddim;iloop++) {

			if (iloop >= kloop) {
				if (
					(expecteddim > 100) ||
					( (iloop & 0b11) == 0)
				) printf("%i ",expecteddim-iloop);
			}

			if (iloop > kloop) {
				// ci1=( a(k-1,k)*a(i,k-1)-a(k-1,k-1)*a(i,k) ) / adenom
				// set only usage flags for the polynomials used
				// for writing into here

				int8_t gerror1=0,gerror2=0;
				{
					{
						{
							if (azpolyMul_TAB(
								t1,
								*Mkmd->getPointer1based_YX(kloop-1,kloop),
								*Mkmd->getPointer1based_YX(iloop,kloop-1)
							) != 0) gerror1=1;

						}
						{
							if (azpolyMul_TAB(
								t2,
								*Mkmd->getPointer1based_YX(kloop-1,kloop-1),
								*Mkmd->getPointer1based_YX(iloop,kloop)
							) != 0) gerror2=1;

						}
					}
				}

				error=0;
				error += t1.subPoly(t2);

				if (
                    (error != 0) ||
                    (gerror1 != 0) ||
                    (gerror2 != 0)
                ) {
                    RET2(-1)
                }

                error += azpolyDiv_TRABLH(
                    ci1,rem,t1,*adenom,VARIDX,help);
                if ( (error == 0) && (rem.anzterms != 0) ) {
                    error--; // error
                    LOGMSG("\nError. c0 not fraction free\n");
                    printf("\nt1 / ad: |"); t1.ausgabe(stdout,-1);
                    printf("|\n / |"); adenom->ausgabe(stdout,-1);
                    printf("|\n");
                    exit(99);
                }

                if (error != 0) {
                    RET2(-1)
                }

				// ci2=( a(k,k-1)*a(i,k)-a(k,k)*a(i,k-1) ) / adenom

				gerror1=0;
				gerror2=0;
				{
					{
						{
							if (azpolyMul_TAB(
								t1,
								*Mkmd->getPointer1based_YX(kloop,kloop-1),
								*Mkmd->getPointer1based_YX(iloop,kloop)
							) != 0) gerror1=1;

						}
						{
							if (azpolyMul_TAB(
								t2,
								*Mkmd->getPointer1based_YX(kloop,kloop),
								*Mkmd->getPointer1based_YX(iloop,kloop-1)
							) != 0) gerror2=1;

						}
					}
				}

				error=0;
				error += t1.subPoly(t2);

				if (
                    (error != 0) ||
                    (gerror1 != 0) ||
                    (gerror2 != 0)
                ) {
                    RET2(-1)
                }

                error += azpolyDiv_TRABLH(
                    ci2,rem,t1,*adenom,VARIDX,help);
                if ( (error == 0) && (rem.anzterms != 0) ) {
                    error--; // error
                    LOGMSG("\nError. c0 not fraction free\n");
                    printf("\nt1 / ad: |"); t1.ausgabe(stdout,-1);
                    printf("|\n / |"); adenom->ausgabe(stdout,-1);
                    printf("|\n");
                    exit(99);
                }

                if (error != 0) {
                    RET2(-1)
                }

			} // compute ci1, ci2

			for(int32_t jloop=1;jloop<=expecteddim;jloop++) {

				// compute i,j >= kloop+1
				// and special element kloop,kloop

				// if too early => save null polynomial
				if ( !(
					( (iloop == kloop) && (jloop == kloop) ) ||
					( (iloop > kloop) && (jloop > kloop) )
				) ) {
                    if (nullpoly.save(fk) <= 0) {
                        printf("\n  error. saving matrix. recommended deleting stored files.\n");
                        RET2(-1);
                    }

					#ifdef _STOREROWS
					if (fwrite(&nullpolynom,1,sizeof(nullpolynom),frow) != sizeof(nullpolynom) ) {
						LOGMSG("\nError. Row1\n");
						exit(99);
					}
					#endif

					continue;
				}

				// compute k,k
				// k,kis actually an element of the previous
				// matrix, as the algorithm only compues
				// elements of the principal minors, hence
				// having coordinates k+1,k+1 or larger

				if ( (iloop == kloop) && (jloop == kloop) ) {
					// compute k,k
					// akk=(a(k-2)(k-1,k-1)*a(k-2)(k,k)-a(k-1,k)*a(k,k-1)) / a(k-2,k-2)

					gerror1=0;
					gerror2=0;
					{
						{
							{
								if (azpolyMul_TAB(
									t1,
									*Mkmd->getPointer1based_YX(kloop-1,kloop-1),
									*Mkmd->getPointer1based_YX(kloop,kloop)
								) != 0) gerror1=1;

							}
							{
								if (azpolyMul_TAB(
									t2,
									*Mkmd->getPointer1based_YX(kloop-1,kloop),
									*Mkmd->getPointer1based_YX(kloop,kloop-1)
								) != 0) gerror2=1;

							}
						}
					}

					error=0;
					error += t1.subPoly(t2);
					if (
                        (error != 0) ||
                        (gerror1 != 0) ||
                        (gerror2 != 0)
                    ) {
                        RET2(-1)
                    }

                    error=0;
                    error += azpolyDiv_TRABLH(
                        t2,rem,t1,*adenom,VARIDX,help);
                    if ( (error == 0) && (rem.anzterms != 0) ) {
                        error--; // error
                        LOGMSG("\nError. c0 not fraction free\n");
                        printf("\nt1 / ad: |"); t1.ausgabe(stdout,-1);
                        printf("|\n / |"); adenom->ausgabe(stdout,-1);
                        printf("|\n");
                        exit(99);
                    }

					if (error != 0) {
                        RET2(-1)
                    }

					if (t2.save(fk) <= 0) {
                        printf("\n  error. file. saving., recommended deleting.");
                        RET2(-1);
					}

					#ifdef _STOREROWS
					if (fwrite(&notnullpolynom,1,sizeof(nullpolynom),frow) != sizeof(nullpolynom)) {
						LOGMSG("\nError. Row2\n");
						exit(99);
					}
					t2.save(frow);
					#endif

					continue;
				}

				// entry to compute

				// a(k)ij = a(k-2)ij*c0+akj*ci1+a(k-1,j)*ci2
				//			--------------------------------------
				//				a(k-2)(k-2,k-2) = adenom

				gerror1=0;
				gerror2=0;
				{
					{
						{
							if (azpolyMul_TAB(
								t1,
								*Mkmd->getPointer1based_YX(iloop,jloop),
								c0
							) != 0) gerror1=1;

						}
						{
							if (azpolyMul_TAB(
								t2,
								*Mkmd->getPointer1based_YX(kloop,jloop),
								ci1
							) != 0) gerror2=1;

						}
					}
				}

				error=0;
				error += t1.addPoly(t2);

				if (
                    (error != 0) ||
                    (gerror1 != 0) ||
                    (gerror2 != 0)
                ) {
                    RET2(-1)
                }

				error += azpolyMul_TAB(
					t2,
					*Mkmd->getPointer1based_YX(kloop-1,jloop),
					ci2
				);
				error += t1.addPoly(t2);

                error += azpolyDiv_TRABLH(
                    t2,rem,t1,*adenom,VARIDX,help);
                if ( (error == 0) && (rem.anzterms != 0) ) {
                    error--; // error
                    LOGMSG("\nError. c0 not fraction free\n");
                    printf("\nt1 / ad: |"); t1.ausgabe(stdout,-1);
                    printf("|\n / |"); adenom->ausgabe(stdout,-1);
                    printf("|\n");
                    exit(99);
                }

				if (error != 0) {
                    RET2(-1)
                }

				if (t2.save(fk) <= 0) {
                    printf("\n  error. storing. recommended file delete.");
                    RET2(-1);
				}

			}  // xloop
		}  // yloop

		fclose(fk);

	} // kloop

	// computation ready => load last matrix
	// by dimension design the M(dim-1) was computed
	// and holds one entry at dim-1,dim-1 (so not dependent
	// on GLOBALMULTISTEP)

	getMatrixFileName(fn,expecteddim-1);
	Mkmd->load(fn,expecteddim-1,expecteddim-1);
	res.copyFrom(
		*Mkmd->getPointer1based_YX(expecteddim,expecteddim)
	);

	RET2(0) // success

}

// MatrixPolynom
AZPoly* MatrixPolynom::getPointer1based_YX(
	const int32_t ay1,
	const int32_t ax1
) {
	if (
		(ax1 >= 1) && (ax1 <= dim) &&
		(ay1 >= 1) && (ay1 <= dim)
	) return &entryYX[(ay1-1)*dim+(ax1-1)];

	LOGMSG3("\nError. Access to non-existent matrix entry y%i,x%i\n",ay1,ax1);
	exit(99);
}

int8_t MatrixPolynom::load(
	const char* afn,
	const int32_t loady,
	const int32_t loadx
) {
	/*

		if loady>=0 & loadx >= 0 => only that element is to
		be actually loaded, the rest of the materix is set to
		constant0

	*/

	FILE *f=fopen(afn,"rb");
	if (!f) {
		LOGMSG2("\nError. Not able to open file for matrix loading |%s|.\n",afn);
		exit(99);
	}

	int32_t d;
	READ32(f,d)

	// sets mgr to the same value, but this way calling
	// setDimension one is forced to provide a valid manager
	if (!entryYX) {
		setDimension(d); // sets dim and allocates memory
	} else {
		if (d != dim) {
			setDimension(d);
		} else {
			// reConstructor polynomial entries
			for(int32_t y=0;y<dim;y++) {
				for(int32_t x=0;x<dim;x++) {
					// index correct as dim <= 2^14
					entryYX[y*dim+x].setToZero();
				}
			}
		}
	}

	for(int32_t y=0;y<dim;y++) {
		for(int32_t x=0;x<dim;x++) {
            if (entryYX[y*dim+x].load(f) <= 0) {
                printf("\nerror reading matrix\n");
                exit(99);
            }
		} // x
	} // y

	fclose(f);

	return 1;
}

int8_t MatrixPolynom::save(const char* afn) {
	FILE *f=fopen(afn,"wb");
	if (!f) {
		LOGMSG("\nError. Not able to open file for matrix storage.\n");
		exit(99);
	}

	WRITE32(f,dim)

	for(int32_t y=0;y<dim;y++) {
		for(int32_t x=0;x<dim;x++) {
            if (entryYX[y*dim+x].save(f) <= 0) {
                printf("\nerror writing matrix. dfiscard fiel.\N");
                exit(99);
            }
		} // x
	} // y

	fclose(f);

	return 1;
}

int8_t MatrixPolynom::entryZero_YX(
	const int32_t ay,
	const int32_t ax
) {
	if ( (ax < 0) || (ax >= dim) || (ay < 0) || (ay >= dim) ) {
		LOGMSG("\nError. Access to non-existent matrix element\n");
		exit(99);
	}

	// ax,ay,dim <= 2^14
	if (entryYX[ay*dim + ax].anzterms <= 0) return 0;

	return 1;
}

void MatrixPolynom::ausgabe(FILE* f) {
	fprintf(f,"dim=%i\n",dim);
	if (!entryYX) {
		fprintf(f,"no poinzers.\n");
		return;
	}

	DynSlowString one;
	for(int32_t y=0;y<dim;y++) {
		fprintf(f,"[y%i]: ",y);
		for(int32_t x=0;x<dim;x++) {
			if ( (x % 8) == 0) {
				fprintf(f,"[x%i] ",x);
			}
			// x,y,dim <= 2^14
			entryYX[y*dim+x].ausgabe(f,-1);
			fprintf(f," / ");
		}

		fprintf(f,"\n");
	} // y
}

MatrixPolynom::MatrixPolynom() {
	dim=0;
	entryYX=NULL;
}

MatrixPolynom::~MatrixPolynom() {
	delete[] entryYX; // allocated with new
	// the polynomials themselves are released
	// when the appropirate termManager is destroyed
}

void MatrixPolynom::setDimension(
	const int32_t adim
) {
	if (adim >= MAXSYLDIMENSION) {
		LOGMSG("\nError. Matrix dimension must not exceed 2^14\n");
		exit(99);
	}

	if (entryYX) {
		if (dim != adim) {
            // delete always if not the correct size as
			delete[] entryYX;
			entryYX=NULL;
		}
	}
	dim=adim;
	if (!entryYX) {
		entryYX=new AZPoly[(int64_t)adim*adim];
	}
	if (!entryYX) {
		LOGMSG("\nError. Memory. MatrixPolynom.\n");
		exit(99);
	}

	setConstant0polynom();

}

void MatrixPolynom::setConstant0polynom(void) {
	if (!entryYX) return;

	for(int32_t y=0;y<dim;y++) {
		for(int32_t x=0;x<dim;x++) {
			entryYX[y*dim+x].setToZero();
		}
	}
}

int32_t SymPoly::get_maxBexponent(void) {
    int32_t erg=0;

    for(int32_t k=0;k<anzterms;k++) {
        for(int32_t m=0;m<PERIOD;m++) {
            if (terms[k].bexponent[m] > erg) erg=terms[k].bexponent[m];
        } // m
    } // k

    return erg;

}

int32_t convertStrToBTerm_TA(
    SymTerm& erg,
    const char* bstr
) {
    erg.setToZero();
    erg.setFactor( 1 );
    char* all=new char[strlen(bstr)+4];
    char *all0=all;
    sprintf(all,"%s,",bstr);

    int32_t error=0,bctr=0;
    for(int32_t k=0;k<PERIOD;k++) {
        char *p=strchr(all,',');
        if (!p) {
            // too few
            // set rest of exponents to 0
            for(int32_t k2=k;k2<PERIOD;k2++) erg.set_bExponent_id_w(k2,0);
            break; // too few
        }
        p[0]=0;
        int64_t w;
        if (sscanf(all,"%i",&w) != 1) { error=-1; break; }
        erg.set_bExponent_id_w(k,w);
        all=p+1; // next

    } // k

    delete[] all0;

    return 0; // success
}

int32_t constructBetaBin_TI_prodA0A1B_FsumT(
    Equation0& erg,
    // tgtid for equation
    const int32_t aid,

    // for binomial cyclic product
    const int32_t a0,
    const int32_t a1,
    const int32_t adiff,    // for cycluic product

    // for cyclic sum
    const int32_t fully,    // > 0 => fully symmetric function form term, else cyclic
    const char* aterm      // term for cyclic sum
) {
    int32_t error=0;

    // constructs the following SymPoly's
    // A: PROD[i=a0..a1} [ b_i^2-b_(i+diff % p)^2 ]
    // B: PROD[i=a0..a1} [ b_i  -b_(i+diff % p)   ]
    // C: PROD[i=a0..a1} [ b_i  +b_(i+diff % p)   ]
    // D: A repeated_itered
    // it holds A=B*C and A=D => D=B*C
    // now divide
    // D/B = C
    // if D/B is remainder-free, one has a new equation

    // then compute
    // E: sum{i=0..period-1} T(i) (T=term from TERMSTR)
    // equation: D/B * E = C*E
    // 0 = C*E - D/B*E
    // replace_iter and express by elem sym func

    // A: PROD[i=a0..a1} [ b_i^2-b_(i+diff % p)^2 ]
    // CAVE: if a0..a1=0..period-1 could use cyclicProduct_2sum function
    SymPoly A,sp1,cp;
    SymTerm T;
    A.setToOne();

    for(int32_t k=a0;k<=a1;k++) {
        int32_t kn=(k+adiff) % PERIOD;
        sp1.setToZero();
        T.setToZero();
        T.setFactor(1);
        T.set_bExponent_id_w(k,2);
        error += sp1.addTerm(T);

        T.setToZero();
        T.setFactor(1);
        T.set_bExponent_id_w(kn,2);
        error += sp1.subTerm(T);

        cp.copyFrom(A);
        error += polyMul_TAB(A,cp,sp1);

    } // k

    if (error != 0) return -1;

    // B: PROD[i=a0..a1} [ b_i  -b_(i+diff % p)   ]
    SymPoly B;
    B.setToOne();

    for(int32_t k=a0;k<=a1;k++) {
        int32_t kn=(k+adiff) % PERIOD;
        sp1.setToZero();
        T.setToZero();
        T.setFactor(1);
        T.set_bExponent_id_w(k,1);
        error += sp1.addTerm(T);

        T.setToZero();
        T.setFactor(1);
        T.set_bExponent_id_w(kn,1);
        error += sp1.subTerm(T);

        cp.copyFrom(B);
        error += polyMul_TAB(B,cp,sp1);

    } // k

    error += replace_repeated_iter(B); // as B is used as divisor of D

    if (error != 0) return -1;

    // C: PROD[i=a0..a1} [ b_i  +b_(i+diff % p)   ]
    // (not time-relevantr, so no need to shift in previous loop)
    SymPoly C;
    C.setToOne();

    for(int32_t k=a0;k<=a1;k++) {
        int32_t kn=(k+adiff) % PERIOD;
        sp1.setToZero();
        T.setToZero();
        T.setFactor(1);
        T.set_bExponent_id_w(k,1);
        error += sp1.addTerm(T);

        T.setToZero();
        T.setFactor(1);
        T.set_bExponent_id_w(kn,1);
        error += sp1.addTerm(T);

        cp.copyFrom(C);
        error += polyMul_TAB(C,cp,sp1);

    } // k

    error += replace_repeated_iter(C);

    if (error != 0) return -1;

    // D: A repeated_itered
    SymPoly D;
    error += replace_repeated_iter_TA(D,A);
    if (error != 0) return -1;

    // D/B = C
    // if D/B is remainder-free, one has a new equation
    SymPoly DB,rem;
    DivHlpSym dh;
    // in case of remainder-free => not relevant which variable is set the lexicographic highest, so just use a0
    // which occurs in B, so for remiander-freeness must also occur in D
    error += polyDiv_TRABLH(DB,rem,D,B,a0,dh);

    if (error != 0) return -1;
    if (rem.anzterms != 0) {
        printf("\n  error. division not remainder-free\n");
        return -1;
    }

    // temp equation: 0=C-D/B
    // E: sum{i=0..period-1} T(i) (T=term from TERMSTR)
    // full or cyclic sum 0..period-1 here
    SymPoly E;
    error += convertStrToBTerm_TA(T,aterm);
    if (error != 0) return -1;

    if (fully > 0) {
        //printf("\nfull sum "); T.ausgabe(stdout);
        error += createFullySymmetric_TA(E,T);
    } else {
        //printf("\ncycliic sum "); T.ausgabe(stdout);
        error += cyclicSum_TA(E,T);
    }
    //printf("\n  E => "); E.ausgabe(stdout,-1);

    if (error != 0) return -1;

    // equation: D/B * E = C*E
    // 0 = C*E - D/B*E
    // replace_iter and express by elem sym func
    SymPoly sp2;
    error += polyMul_TAB(sp1,C,E);
    // 0 = sp1 - D/B*E
    error += polyMul_TAB(sp2,DB,E);
    // 0 = sp1 - sp2
    error += sp1.subPoly(sp2);
    // 0 = sp1
    error += replace_repeated_iter(sp1);
    // express
    AZPoly az;
    SymPoly non;
    error += expressByElem_Texp_Tnon_A(az,non,sp1);
    // 0=AZ+NON
    // -AZ = NON
    az.invertSign();
    // AZ = NON
    erg.setAll_ABC(aid,az,non);
    erg.fullSimplify();

    return error;
}


// main
int32_t main(int32_t argc,char** argv) {
    // set
    initConstants();

    time0=clock();

    flog=fopen("elemsymfunc-core.log.txt","at");
    if (!flog) {
        printf("\nerror. open log file\n");
        exit(99);
    }
    fprintf(flog,"\n\n-------------------------\n");

    LOGMSG("\n\nSTATUS:\n");
    LOGMSG2("\n  Bareiss %i-step\n",BAREISS_STEP);

    initSigma(PERIOD);

    // principal function
    periodAll_interactive();

    // finalize
    double d=clock()-time0; d /= CLOCKS_PER_SEC;
    LOGMSG2("\n\nduration %.0lf sec\n",(double)d);

    fclose(flog);

    return 0;

}


