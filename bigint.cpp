#ifndef _BIGINT
#define _BIGINT

/*

	INFO: in case of an error, 0 is set to the resulting
	variable(s) to have highestusedidx in a valid range
	(so accumulation of errors can be performed)

	big integer on base 10^9

	it shall not be assumed that calling a function
	like bigIntMul etc. keeps the arguments constant
	at all times (it does at return, but signs might
	have been changed during computation)

	if using this in parallel applications, care must
	be taken if passing the same physical object in different
	parallel threads

*/


// consts

const int64_t UINT62MAX=((int64_t)1 << 62) - 1;
const int32_t INT31MAX=0b01111111111111111111111111111111;
const int32_t MAXBIGINTDIGITS=512;
// BASE*BASE must fit into 63 bits for multiplication purposes
const int32_t BIGINTBASE=1000000000; // 10^9 per 4-byte block
const uint64_t UINT64MAX=0xFFFFFFFFFFFFFFFF;

int32_t tenpower[] = {
	1,
	10,
	100,
	1000,
	10000,
	100000,
	1000000,
	10000000,
	100000000,
	1000000000
};


// structs

struct BigInt {
	int8_t vorz; // +1, -1, 0
	int32_t highestusedidx;
	int32_t digits[MAXBIGINTDIGITS]; // the higher the index
	// the higher the order

	BigInt();
	BigInt(const BigInt&);

	void save(FILE*);
	void load(FILE*);
	void set_int64(const int64_t);
	void setTo10power(const int32_t);
	int32_t addTo(BigInt&);
	int32_t subTo(BigInt&);
	void getStr(DynSlowString&);
	void getStrLZ(DynSlowString&);
	void setToZero(void);
	void copyFrom(BigInt&);
	void copyFrom(const BigInt&);
	int32_t shiftLeft10(const int32_t);
	int32_t shiftLeftBase(const int32_t);
	int32_t tendigitcount(void);
	void ausgabe(FILE*);
	int32_t inc(void);
	int32_t dec(void);
	int32_t setstr(const char*);
	int8_t isZero(void);
	int8_t isPositiveOne(void);
	void invertSign(void);
	int32_t getDig2(int64_t&);
	int8_t convertToInt64(int64_t&);
	void absolute(void);

	BigInt& operator=(const BigInt);
};


// forward

inline int32_t bigintDiv_TRAB(BigInt&,BigInt&,BigInt&,BigInt&);
inline int32_t bigintMul_TAB(BigInt&,BigInt&,BigInt&);
//int32_t  bigintMul_digit_TDA(BigInt&,int64_t&,BigInt&);
inline int32_t bigintMul_digit_TDA(BigInt&,int32_t&,BigInt&);
inline int32_t bigintAdd_TAB(BigInt&,BigInt&,BigInt&);
inline int32_t bigintAdd_abs_TAB(BigInt&,BigInt&,BigInt&);
inline int32_t bigintSub_TAB(BigInt&,BigInt&,BigInt&);
inline int32_t bigintSub_abs_TAB(BigInt&,BigInt&,BigInt&);
inline int32_t bigintSub_abs_ovgl_TAB(BigInt&,BigInt&,BigInt&);
inline int32_t bigintPow_TAE(BigInt&,BigInt&,const int32_t);
inline int32_t bigintGcd_abs_TAB(BigInt&,BigInt&,BigInt&);

inline int8_t bigintVgl_AB(BigInt&,BigInt&);
inline int8_t bigintVgl_abs_AB(BigInt&,BigInt&);


// globals

inline int32_t sum_int32t(const int32_t,const int32_t);
extern int64_t gmbiginthighestusedidx;


// routines

// BigInt
void BigInt::absolute(void) {
    if (vorz < 0) vorz=1;
}

int8_t BigInt::convertToInt64(int64_t& erg) {
    if (highestusedidx >= 2) return 0;
    // max value now BIGINTBASE^2,fits into int64_t

    erg=digits[0];
    if (highestusedidx == 1) erg=erg+digits[1]*BIGINTBASE;

    return 1;
}

int32_t BigInt::getDig2(int64_t& erg) {
    // return <=0 : error
    // return 1: fits

    if (highestusedidx >= 2) return -1; // to big
    // now: number is at most dig[1]*BASE+dig[0]
    if (BIGINTBASE != 1000000000) {
        printf("\nerror. bigint::getDig2 not possible with base different to 1000000000\n");
        exit(9);
    }

    uint64_t a=digits[0];
    if (highestusedidx == 1) a += (uint64_t)digits[1] * BIGINTBASE;

    if (a >= (uint64_t)1 << 62) {
        printf("\nerror. getDig2 too large.\n");
        exit(99);
    }

    // fits into int64_t

    erg=(int64_t)a;

    return 1;
}

int32_t BigInt::setstr(const char* as) {
    if (!as) return -1; // error

	// no white spaces
	if (BIGINTBASE != 1000000000) {
		LOGMSG("\nError. BigInt::setstr only working in base 10^9\n");
		exit(99);
	}

	// only at most 1 sign
    int32_t signs=0;
    int32_t kl=strlen(as);

    if (kl <= 0) return -1; // error

    for(int32_t k=0;k<kl;k++) {
        if (
            (as[k] == '-') ||
            (as[k] == '+')
        ) {
            signs++;
            if (signs > 1) break;
        }
    }

    if (signs > 1) {
        printf("\nbigint:: cannot set more than 1 sign\n");
        printf("|%s|\n",as);

        return -1; // erfror
    }

	//printf("\nbigint:setstr(%s)\n",as);
	highestusedidx=0;
	vorz=1;
	int32_t s0=0;
	if (as[0] == '-') {
		vorz=-1;
		s0=1;
	} else if (as[0] == '+') {
		s0=1;
	}

	//printf("1) (as|%s|) s0=%i (vorz%i)\n",as, s0,vorz);

	int32_t p=strlen(as)-1;

	if (p >= (9*MAXBIGINTDIGITS-4)) {
        printf("\nerror. bigint:setstr too many digits\n");
        exit(99);
	}
	int32_t setidx=0;
	highestusedidx=0;
	int8_t notzero=0;

	while (p >= s0) {
		// [p-8..p] is to be set into digit[setidx]
		int32_t pl=p-8;
		if (pl < s0) pl=s0;
		int64_t w;
		char tmp[256];
		for(int32_t i=pl;i<=p;i++) {
			tmp[i-pl]=as[i];
		}
		tmp[p-pl+1]=0;
		w=0;
		//printf("block |%s| ",tmp);
		int64_t DIGP=1;
		for(int64_t i=(strlen(tmp)-1);i>=0;i--) {
			int64_t dig=tmp[i]-'0';
			if ( (dig < 0) || (dig > 9) ) {
				return -1;
				//LOGMSG3("\nError. Did not recognize digits %s: %I64d\n",tmp,dig);
				//exit(99);
			}
			if (dig != 0) w += dig*DIGP;
			DIGP *= 10;
		};

		if ( (w < 0) || (w >= BIGINTBASE) ) {
			LOGMSG("\nError. BigInt. not recognized2\n");
			exit(99);
		}

		digits[setidx]=w;
		if (w != 0) {
			highestusedidx=setidx;
			notzero=1;
		}
		setidx++;
		p=pl-1;
	} // p

	//printf("2) (as|%s|) s0=%i (vorz%i)\n",as, s0,vorz);

	if (notzero == 0) {
		vorz=0;
		digits[0]=0;
	}

	/*
	printf("\n  set(|%s|)\nBigInt:setstr result: ",as);
	ausgabe(stdout);
	printf("\n");
	*/

	return 0;

}

int32_t BigInt::shiftLeftBase(const int32_t toleft) {
	int32_t error=0;

	if ( (highestusedidx + toleft) >= (MAXBIGINTDIGITS-4) ) {
		setToZero();
		return -1; // error

		LOGMSG("\nError. Overflow bigInt shift\n");
		exit(99);
	}

	for(int32_t i=highestusedidx;i>=0;i--) {
		digits[i+toleft]=digits[i];
	} // i
	for(int32_t i=0;i<toleft;i++) {
		digits[i]=0;
	}
	highestusedidx += toleft;

	if (error != 0) {
		setToZero();
	}

	return error;
}

int32_t BigInt::shiftLeft10(const int32_t toleft) {
	int32_t error=0;

	// in 10-base, not BIGINTBASE
	if (toleft == 0) return 0;
	if (toleft < 0) {
		setToZero();

		return -1; // error
		LOGMSG("\nError. Implementation. ::shiftLeft woth negative argument\n");
		exit(99);
	}

	// now if toleft >= 9
	// shift in digits[] array positions to the left
	// and fill with 0
	int32_t shift=toleft;
	int32_t how9=0;
	if (shift == 9) { how9=1; shift=0; }
	else {
		how9=shift / 9; // downward rounding
		shift -= how9*9;
	}
	//printf("\ntoleft %i how9 %i shift %i\n",toleft,how9,shift);

	if ( (shift < 0) || (shift > 8) ) {
		// error
		setToZero();

		return -1;

		LOGMSG("\nError. Implementation. bigint::ShiftLeft div 9\n");
		exit(99);
	}

	if ( (highestusedidx + how9) >= (MAXBIGINTDIGITS-2) ) {
		setToZero();

		return -1;

		LOGMSG("\nError. bigInt::ShiftLeft overflow.\n");
		exit(99);
	}

	if (how9 > 0) {
		for(int32_t i=(highestusedidx+how9);i>=how9;i--) {
			digits[i]=digits[i-how9];
		}
		for(int32_t i=(how9-1);i>=0;i--) {
			digits[i]=0;
		}

		highestusedidx += how9;
	} // how9

	if (shift > 0) {
		// then if sth is left
		// do a multiplication by 10^appropriate
		BigInt tmp;
		error += bigintMul_digit_TDA(tmp,tenpower[shift],*this);
		copyFrom(tmp);
	}

	if (error != 0) {
		setToZero();
	}

	return error;
}

void BigInt::load(FILE* f) {
	if (!f) {
		LOGMSG("\nError. Implementation. BigInt::load file\n");
		exit(99);
	}

	int32_t r=0;
	r += fread(&vorz,1,sizeof(vorz),f);
	r += fread(&highestusedidx,1,sizeof(highestusedidx),f);
	if (r != (sizeof(vorz)+sizeof(highestusedidx)) ) {
		LOGMSG("\nError. Reading bigint. Probably invalid matrix file. Deleting recommended.\n");
		exit(99);
	}

	if (highestusedidx >= MAXBIGINTDIGITS) {
		LOGMSG("\nerror. not enough precision to read BigInt.\n");
		LOGMSG3("present=%i, needed=%i\n",MAXBIGINTDIGITS,highestusedidx);
		exit(99);
	}

	for(int32_t i=0;i<=highestusedidx;i++) {
		if (fread(&digits[i],1,sizeof(digits[i]),f) != sizeof(digits[i])) {
			LOGMSG("\nError. Reading bigint. Probably invalid matrix file. Deleting recommended.\n");
			exit(99);
		}
	} // i
}

void BigInt::save(FILE* f) {
	if (!f) {
		LOGMSG("\nError. Implementation. BigInt::save file\n");
		exit(99);
	}

	int32_t r=0;
	r += fwrite(&vorz,1,sizeof(vorz),f);
	r += fwrite(&highestusedidx,1,sizeof(highestusedidx),f);
	if (r != (sizeof(vorz)+sizeof(highestusedidx)) ) {
		LOGMSG("\nError. Writing bigint. Probably invalid matrix file. Deleting recommended.\n");
		exit(99);
	}
	for(int32_t i=0;i<=highestusedidx;i++) {
		if (fwrite(&digits[i],1,sizeof(digits[i]),f) != sizeof(digits[i])) {
			LOGMSG("\nError. Writing bigint. Probably invalid matrix file. Deleting recommended.\n");
			exit(99);
		}
	} // i
}

void BigInt::setTo10power(const int32_t pow) {

	if (BIGINTBASE != 1000000000) {
		LOGMSG("\nError. Implementation. BigInt::setTo10power needs base to be 10^9\n");
		exit(99);
	}

	highestusedidx=-1;
	digits[0]=0;
	vorz=1;

	int32_t p=pow;
	while (p >= 9) {
		highestusedidx++;
		if (highestusedidx > gmbiginthighestusedidx) gmbiginthighestusedidx=highestusedidx;
		digits[highestusedidx]=0;
		p -= 9;
	} // while

	highestusedidx++;
	if (highestusedidx > gmbiginthighestusedidx) gmbiginthighestusedidx=highestusedidx;
	digits[highestusedidx]=tenpower[p];

}

int32_t BigInt::dec(void) {
	int32_t error=0;

	if (vorz <= 0) {
		LOGMSG("\nError. BigInt. Decrement only works on positive values.\n");
		exit(99);
	}

	for(int32_t i=0;i<=highestusedidx;i++) {
		if (digits[i] > 0) {
			digits[i]--;
			break;
		}
		digits[i]=BIGINTBASE - 1;
	}

	while (highestusedidx > 0) {
		if (digits[highestusedidx] == 0) highestusedidx--;
		else break;
	}

	// is value zero
	vorz=0;
	for(int32_t i=0;i<=highestusedidx;i++) {
		if (digits[i] != 0) {
			vorz=1;
			break;
		}
	} // check

	if (error != 0) {
		setToZero();
	}

	return error;

}

int32_t BigInt::inc(void) {
	int32_t error=0;

	if (vorz < 0) {
		LOGMSG("\nError. BigInt. Increment only works on positive values.\n");
		exit(99);
	}

	vorz=1; // could have been 0 from the start

	int64_t carryover=1; // increment

	for(int32_t i=0;i<=highestusedidx;i++) {
		int64_t sum=carryover + (int64_t)digits[i];
		if (sum >= BIGINTBASE) {
			sum -= BIGINTBASE;
			if (sum >= BIGINTBASE) {
				LOGMSG("\nError. BigInt.inc\n");
				exit(99);
			}
			digits[i]=sum;
			carryover=1;
		} else {
			digits[i]=sum;
			if (error != 0) {
				setToZero();
			}

			return error;
		}
	} // i

	if (carryover > 0) {
		if ( (highestusedidx+1) >= MAXBIGINTDIGITS) {
			setToZero();
			return -1; // error,

			LOGMSG("\nError. BigInt. Increment. Overflow.\n");
			exit(99);
		}
		highestusedidx++;
		if (highestusedidx > gmbiginthighestusedidx) gmbiginthighestusedidx=highestusedidx;
		digits[highestusedidx]=1;
	}

	if (error != 0) {
		setToZero();
	}

	return error;

}

void BigInt::setToZero(void) {
	vorz=0;
	highestusedidx=0;
	digits[0]=0;
}

void BigInt::getStr(DynSlowString& erg) {
	erg.setEmpty();
	if (vorz == 0) {
		erg.add("+0");
		return;
	}

	char tmp[256];
	if (vorz < 0) erg.add("-");

	for(int32_t i=highestusedidx;i>=0;i--) {
		if (i != highestusedidx) {
			if (BIGINTBASE == 1000000000) {
				sprintf(tmp,"%09i",digits[i]);
			} else {
				sprintf(tmp,"/%i/",digits[i]);
			}
		}
		else sprintf(tmp,"%i",digits[i]);
		erg.add(tmp);
	}
}

void BigInt::getStrLZ(DynSlowString& erg) {
	erg.setEmpty();
	if (vorz == 0) {
		erg.add("+0");
		return;
	}

	char tmp[256];
	if (vorz < 0) erg.add("-");

	for(int32_t i=highestusedidx;i>=0;i--) {
		if (BIGINTBASE == 1000000000) {
			sprintf(tmp,"%09i",digits[i]);
		} else {
			sprintf(tmp,"/%i/",digits[i]);
		}

		erg.add(tmp);
	}
}

int32_t BigInt::subTo(BigInt& av) {
	int32_t error=0;

	BigInt res;
	error += bigintSub_TAB(res,*this,av);
	copyFrom(res);

	if (error != 0) {
		res.setToZero();
	}

	return error;
}

void BigInt::copyFrom(BigInt& av) {
	vorz=av.vorz;
	highestusedidx=av.highestusedidx;
	for(int32_t i=0;i<=highestusedidx;i++) {
		digits[i]=av.digits[i];
	}

}

void BigInt::copyFrom(const BigInt& av) {
	vorz=av.vorz;
	highestusedidx=av.highestusedidx;
	for(int32_t i=0;i<=highestusedidx;i++) {
		digits[i]=av.digits[i];
	}

}

int32_t BigInt::addTo(BigInt& av) {
	int32_t error=0;

	BigInt res;
	error += bigintAdd_TAB(res,*this,av);
	copyFrom(res);

	if (error != 0) {
		res.setToZero();
	}

	return error;
}

BigInt::BigInt() {
	// do nothing as fastest construction
}

BigInt::BigInt(const BigInt& b) {
	vorz=b.vorz;
	highestusedidx=b.highestusedidx;
	for(int32_t i=0;i<=highestusedidx;i++) {
		digits[i]=b.digits[i];
	}
}

void BigInt::set_int64(const int64_t aw) {
	// as one negates aw if negative it can happen that it
	// will overflow/underflow as positive and negative inetegrs
	// in int64_t are not identical in range
	// so test beforehand by higher precision number type

	int64_t w2=(int64_t)BIGINTBASE * BIGINTBASE;
	if ( (aw >= w2) || (aw <= (-w2) ) ) {
		LOGMSG2("\nError. Not possible to set value ||>= %I64d directly\n",w2);
		exit(99);
	}

	if (aw == 0) {
		setToZero();
	} else {
		int64_t setv;
		if (aw > 0) {
			vorz=1;
			setv=aw;
		} else {
			vorz=-1;
			// ATTN: can overflow if aw is too negative
			// (the most negative int64_t there is)
			// then the positive one does not exist and
			// it stays negative, hence the pre-test
			// at the beginning with float128
			setv=-aw;
		}

		if (setv <= 0) {
			LOGMSG("\nImplementation error. setv non-positve.\n");
			exit(99);
		}

		digits[0]=setv % BIGINTBASE;
		digits[1]=(setv-digits[0]) / BIGINTBASE;

		if (digits[1] >= BIGINTBASE) {
			LOGMSG("\nError. Setting BigInt too large\n");
			exit(99);
		}

		if ( (digits[0] < 0) || (digits[1] < 0) ) {
			LOGMSG("\nImplementation error. Setting BigInt negative indices\n");
			LOGMSG2("\n%I64d\n",aw);
			exit(99);
		}

		if (digits[1] != 0) highestusedidx=1;
		else highestusedidx=0;
	}
}

int32_t bigintMul_digit_TDA(
	BigInt& res,
	int32_t& D,
	BigInt& A
) {
	res.setToZero();

	int32_t error=0;

	// D is a BIGINTBASE-digit
	if ( (D < 0) || (D >= BIGINTBASE) ) {
		LOGMSG2("\nError. Implementation. bigintmul_digit outside range 0..BIGINTBASE[ : %i\n",D);
		exit(99);
	}

	if (D == 0) {
		res.setToZero();
		return 0;
	}

	int32_t signD;
	if (D < 0) signD=-1;
	else if (D == 0) signD=0;
	else if (D > 0) signD=1;

	// conservatively judged
	// current highstusedidx is smaller, so adding 1 does not overflow
	if ( (A.highestusedidx+1) >= MAXBIGINTDIGITS) {
		// treat this as an erro and return -1
		res.setToZero();
		return -1; // error

		LOGMSG("\nError. Overflow. bigintMul\n");
		exit(99);
	}

	for(int32_t i=0;i<=(A.highestusedidx+1);i++) {
		res.digits[i]=0; // has to be set to zero as continues addition occurs later
	}

	for(int32_t idxa=0;idxa<=A.highestusedidx;idxa++) {
		int64_t w=(int64_t)A.digits[idxa] * D;
		// add w % BASE to res[i]
		// add w / BASE to res[i+1] and shift carry-overs

		// adds value 0 < DIG < BIGINTBASE at res[i] and
		// considers carry-over
		#define ADDAT(IDX,DIG) \
		{\
			int64_t carryover=0;\
			for(int32_t lokali=(IDX);lokali<=(A.highestusedidx+1);lokali++) {\
				int64_t sum=carryover + (int64_t)res.digits[lokali];\
				if (lokali == (IDX)) sum += (DIG);\
				\
				if (sum < BIGINTBASE) {\
					carryover=0;\
					res.digits[lokali]=sum;\
					break;\
				} else {\
					carryover=1;\
					sum -= (int64_t)BIGINTBASE;\
					if (sum >= BIGINTBASE) {\
						LOGMSG("\nError. Implementation. Overflow. AddAT/2\n");\
						exit(99);\
					}\
					res.digits[lokali]=sum;\
				}\
				\
				if (carryover == 0) break;\
			}\
			\
			if (carryover != 0) {\
				LOGMSG("\nError. Implementation. ADDAT. Overflow.\n");\
				exit(99);\
			}\
		}

		if (w < BIGINTBASE) {
			if (w != 0) {
				ADDAT(idxa,w)
			}
		} else {
			int64_t w1=w % BIGINTBASE;
			int64_t w2=(w-w1) / BIGINTBASE;
			if (w1 != 0) {
				ADDAT(idxa,w1)
			}
			if (w2 != 0) {
				ADDAT(idxa+1,w2)
			}
		}
	} // idxa

	// res sign: sign(T)*sign(A)

	for(int32_t i=(A.highestusedidx+1);i>=0;i--) {
		if (res.digits[i] > 0) {
			res.vorz=A.vorz * signD;
			res.highestusedidx=i;
			if (res.highestusedidx > gmbiginthighestusedidx) gmbiginthighestusedidx=res.highestusedidx;

			if (error != 0) {
				res.setToZero();
			}

			return error;
		}
	}

	// all digits 0
	res.vorz=0;
	res.highestusedidx=0;

	if (error != 0) {
		res.setToZero();
	}

	return error;

}

void getDigits_abs_TANN(
	BigInt& erg,
	BigInt& src,
	const int32_t i0,
	const int32_t i1
) {
	if (src.vorz == 0) {
		erg.set_int64(0);
		return;
	}

	if (
		(src.highestusedidx < i1) ||
		(i1 >= (MAXBIGINTDIGITS-4))
	) {
		LOGMSG("\nError. BigInt abs TANN. Wrong index\n");
		exit(99);
	}

	erg.vorz=0;
	erg.digits[0]=0;
	erg.highestusedidx=0;

	//printf("[%i..%i]\n",i0,i1);

	for(int32_t i=i0;i<=i1;i++) {
		erg.digits[i-i0]=src.digits[i];
		if (src.digits[i] != 0) {
			if (i > erg.highestusedidx) {
				erg.highestusedidx=i;
				erg.vorz=1;
			}
		}
	} // i

	erg.highestusedidx -= i0;

	// leading zeros are automatically trimmed
}

int32_t bigintMul_TAB(
	BigInt& res,
	BigInt& A,
	BigInt& B
) {
	res.setToZero();

	// digit-wise
	int32_t error=0;

	if (
		(A.vorz == 0) ||
		(B.vorz == 0)
	) {
		res.setToZero();
		return 0;
	}

	int32_t w=A.highestusedidx+B.highestusedidx;
	if (w >= MAXBIGINTDIGITS) {
		res.setToZero();
		return -1; // error

		// might work for some combinations of numbers at those
		// used indices, but it is easier just to increase the MAXDIGIT
		LOGMSG("\nError. bigintMul. Overflow.\n");
		exit(99);
	}

	for(int32_t i=0;i<=w;i++) res.digits[i]=0;
	res.highestusedidx=w;
	res.vorz=1;

	for(int32_t idxa=0;idxa<=A.highestusedidx;idxa++) {
		BigInt tmp;
		if (A.digits[idxa] == 0) continue;
		if (A.digits[idxa] == 1) {
			tmp.copyFrom(B);
		} else {
			error += bigintMul_digit_TDA(tmp,A.digits[idxa],B);
		}

		// shift result tmp i digits to the right
		if (idxa>0) {
			for(int32_t k=(tmp.highestusedidx+idxa);k>=idxa;k--) {
				tmp.digits[k]=tmp.digits[k-idxa];
			}
			for(int32_t k=(idxa-1);k>=0;k--) tmp.digits[k]=0;
		}
		tmp.highestusedidx += idxa;

		// add to res
		BigInt tmp2;
		error += bigintAdd_abs_TAB(tmp2,res,tmp);
		res.copyFrom(tmp2);

	} // idxa

	res.vorz=A.vorz*B.vorz;
	if (res.highestusedidx > gmbiginthighestusedidx) gmbiginthighestusedidx=res.highestusedidx;

	if (error != 0) {
		res.setToZero();
	}

	return error;

}

int32_t bigintAdd_TAB(
	BigInt& res,
	BigInt& A,
	BigInt& B
) {
	res.setToZero();
	int32_t error=0;

	if (A.vorz==0) {
		if (B.vorz==0) {
			res.setToZero();
			return 0;
		} else {
			res.copyFrom(B);
			return 0;
		}
	} else
	if (B.vorz==0) {
		res.copyFrom(A);
		// is not zero
		return 0;
	}

	// check sign
	if (A.vorz == B.vorz) {
		error += bigintAdd_abs_TAB(res,A,B);
		res.vorz=A.vorz;
		if (error != 0) {
			res.setToZero();
		}

		return error;
	} else {
		if ( (A.vorz > 0) && (B.vorz < 0) ) {
			// erg=a-|b|
			BigInt B2(B);
			B2.vorz=1;
			error += bigintSub_TAB(res,A,B2);
			if (error != 0) {
				res.setToZero();
			}

			return error;
		} else {
			// a < 0, b > 0
			// res = -|a|+b = b-|a|
			BigInt A2(A);
			A2.vorz=1;
			error += bigintSub_TAB(res,B,A2);
			if (error != 0) {
				res.setToZero();
			}

			return error;
		}
	}

	LOGMSG("\nError. Implementation. bigint::add.\n");
	exit(99);

	if (error != 0) {
		res.setToZero();
	}

	return error;

}

int32_t bigintSub_TAB(
	BigInt& res,
	BigInt& A,
	BigInt& B
) {
	res.setToZero();
	int32_t error=0;

	if (A.vorz == 0) {
		// a=0 => a-b=-b
		if (B.vorz == 0) {
			// -b=0
			res.setToZero();
			return 0;
		} else {
			// -b
			res.copyFrom(B);
			res.vorz*=-1;
			return 0;
		}
	} else if (B.vorz == 0) {
		// a-b=a-0
		res.copyFrom(A);
		// A != 0
		return 0;
	}

	if (A.vorz>0) {
		// a>0: a-b=|a|-b
		if (B.vorz>0) {
			// a-b=|a|-|b|
			error += bigintSub_abs_TAB(res,A,B);
			if (error != 0) {
				res.setToZero();
			}

			return error;
		} else {
			// a>0, b<0: a-b=a+|b|
			// a-b=|a|-(-|b|)=|a|+|b|
			error += bigintAdd_abs_TAB(res,A,B);
			if (error != 0) {
				res.setToZero();
			}

			return error;
		}
	} else {
		// a<0, a-b=-|a|-b
		if (B.vorz>0) {
			// a-b=-|a|-|b|=-(|a|+|b|)
			error += bigintAdd_abs_TAB(res,A,B);
			res.vorz*=-1;
			if (error != 0) {
				res.setToZero();
			}

			return error;
		} else {
			// a<0,b<0: a-b=-|a|-b=-|a|-(-|b|)=|b|-|a|
			error += bigintSub_abs_TAB(res,B,A);
			if (error != 0) {
				res.setToZero();
			}
			return error;
		}
	}

	if (error != 0) {
		res.setToZero();
	}

	return error;

	LOGMSG("\nError. Implementation. bigintSub\n");
	exit(99);
}

int8_t bigintVgl_AB(BigInt& A,BigInt& B) {
	if (A.vorz == 0) {
		if (B.vorz == 0) return 0;
		if (B.vorz < 0) return +1;
		if (B.vorz > 0) return -1;
	} else if (A.vorz > 0) {
		if (B.vorz <= 0) return +1;

		// both positive
		if (A.highestusedidx > B.highestusedidx) return +1;
		if (A.highestusedidx < B.highestusedidx) return -1;

		for(int32_t i=A.highestusedidx;i>=0;i--) {
			if (A.digits[i] > B.digits[i]) return +1;
			if (A.digits[i] < B.digits[i]) return -1;
		}

		return 0;
	} else if (A.vorz < 0) {
		if (B.vorz >= 0) return -1;

		// both negative
		if (A.highestusedidx > B.highestusedidx) return -1;
		if (A.highestusedidx < B.highestusedidx) return +1;

		for(int32_t i=A.highestusedidx;i>=0;i--) {
			if (A.digits[i] > B.digits[i]) return -1;
			if (A.digits[i] < B.digits[i]) return +1;
		}

		return 0;
	}

	LOGMSG("\nError. Implementation. bigIntVgl\n");
	exit(99);
}

int8_t bigintVgl_abs_AB(BigInt& A,BigInt& B) {
	// zero is considered, but sign not relevant otherwise
	if (A.vorz == 0) {
		if (B.vorz == 0) return 0;
		return -1; // 0 < |B|
	}

	if (B.vorz == 0) return +1; // |A|>0

	if (A.highestusedidx > B.highestusedidx) return +1;
	if (A.highestusedidx < B.highestusedidx) return -1;

	for(int32_t i=A.highestusedidx;i>=0;i--) {
		if (A.digits[i] > B.digits[i]) return +1;
		if (A.digits[i] < B.digits[i]) return -1;
	}

	return 0;
}

int32_t bigintAdd_abs_TAB(
	BigInt& res,
	BigInt& A,
	BigInt&B
) {
	res.setToZero();

	// the sign of A,B is NOT considered
	// only the value itself is used

	int32_t error=0;

	if (A.vorz==0) {
		if (B.vorz==0) res.setToZero();
		else res.copyFrom(B);
		return 0;
	} else if (B.vorz==0) {
		res.copyFrom(A);
		return 0;
	}

	// neither A nor B are zero
	// just add digits fromlowest to highest with carry-over
	int64_t carryover=0;

	int32_t m=A.highestusedidx;
	if (B.highestusedidx > m) m=B.highestusedidx;

	res.vorz=1; // will be positive as 0 can not happen
	res.highestusedidx=0;
	res.digits[0]=0;

	for(int32_t i=0 /* lowest order */;i<=m;i++) {
		// as digits are UINT32, and adding three numbers can
		// onlyhave at most a carry over of 2,the sum will fit
		// into int64
		int64_t sum=carryover;
		if (i <= A.highestusedidx) sum += (int64_t)A.digits[i];
		if (i <= B.highestusedidx) sum += (int64_t)B.digits[i];

		if (sum >= BIGINTBASE) {
			sum -= BIGINTBASE;
			if (sum >= BIGINTBASE) {
				LOGMSG("\nError. Implementation/1. bigint::add\n");
				exit(99);
			}
			carryover=1;
			res.digits[i]=sum;
		} else {
			carryover=0;
			res.digits[i]=sum;
		}
	} // i

	res.highestusedidx=m;

	if (carryover > 0) {
		res.highestusedidx=m+1;

		if (res.highestusedidx >= MAXBIGINTDIGITS) {
			res.setToZero();

			return -1;

			LOGMSG("\nError. Overflow bigint::Add\n");
			exit(99);
		}
		res.digits[res.highestusedidx]=carryover;
	}

	if (res.highestusedidx > gmbiginthighestusedidx) gmbiginthighestusedidx=res.highestusedidx;

	if (error != 0) {
		res.setToZero();
	}

	return error;

}

int32_t bigintSub_abs_ovgl_TAB(
	BigInt& res,
	BigInt& A,
	BigInt& B
) {
	res.setToZero();
	// computes |A|-|B| and it is assumed (and checked outside)
	// that |A| >  |B|
	// school method

	int32_t error=0;

	if (A.vorz==0) {
		if (B.vorz==0) {
			res.setToZero();
			return 0;
		}
		else {
			// Fehler
			LOGMSG("\nError. Implementation/1. sub_abs_vgl out of range\n");
			exit(99);
		}
	} else if (B.vorz==0) {
		res.copyFrom(A);
		return 0;
	}

	// signed here necessary as variable w below can get negative
	int64_t carryover=0;
	int32_t m=A.highestusedidx;
	if (B.highestusedidx > m) m=B.highestusedidx;

	for(int32_t i=0;i<=m;i++) {
		int64_t vala,valb;
		if (i <= A.highestusedidx) vala=(int64_t)A.digits[i]; else vala=0;
		if (i <= B.highestusedidx) valb=(int64_t)B.digits[i]; else valb=0;
		// as only 32 bits in digits and carry-over are used at max
		// there is no overflow or underflow by adding or subtracting
		// up to 3 numbers int int64_t
		int64_t w=vala - (valb + carryover);
		if (w < 0) {
			w += (int64_t)BIGINTBASE;
			if (w < 0) {
				LOGMSG("\nError. Implementation. subabs_ovgl/2\n");
				printf("\nvala %I64d\nvalb %I64d\ncarry %I64d\n",
					vala,valb,carryover);
				printf("\nw + BIGINTBASE(%I64d): %I64d\n",
					(int64_t)BIGINTBASE,w);
                printf("\nA-B (without size check)\n");
                printf("A(vz%i,hui%i)=",A.vorz,A.highestusedidx);
                A.ausgabe(stdout);
                printf("\nB(vz%i,hui%i)=",B.vorz,B.highestusedidx);
                B.ausgabe(stdout);
				exit(99);
			}
			// carry over
			carryover=1;
			res.digits[i]=w % BIGINTBASE;
		} else {
			if (w >= BIGINTBASE) {
				LOGMSG("\nError. Implementation. bigintsubb_ovgl_abs/2\n");
				exit(99);
			}

			res.digits[i]=w;
			carryover=0;
		}
	}
	res.highestusedidx=m;
	// leading zeros => trim
	while (res.highestusedidx > 0) {
		if (res.digits[res.highestusedidx] == 0) res.highestusedidx--;
		else break;
	}

	if (carryover != 0) {
		res.setToZero();
		return -1;

		LOGMSG2("\nError. Overflow bigint_sub_abs_ovgl. highestusedidx=%i but still carry-over\n",
			res.highestusedidx);
		LOGMSG("\n|A| (ignore sign): "); A.ausgabe(stdout);
		A.ausgabe(flog);
		LOGMSG("\n minus |B| (ignore sign): "); B.ausgabe(stdout);
		B.ausgabe(flog);
		exit(99);
	}

	// as |A| > |B|, |A|-|B| > 0
	res.vorz=1;

	if (error != 0) {
		res.setToZero();
	}

	return error;

}

int32_t bigintSub_abs_TAB(
	BigInt& res,
	BigInt& A,
	BigInt& B
) {
	// only absolute value used, sign is NOT considered
	res.setToZero();
	int32_t error=0;

	int32_t vgl=bigintVgl_abs_AB(A,B);

	if (vgl == 0) {
		res.setToZero();
		return 0;
	}

	if (vgl > 0) {
		// |a| > |b|
		error += bigintSub_abs_ovgl_TAB(res,A,B);
		res.vorz=1;
		if (error != 0) {
			res.setToZero();
		}
		return error;
	} else {
		// |a| < |b|
		error += bigintSub_abs_ovgl_TAB(res,B,A);
		res.vorz=-1;
		if (error != 0) {
			res.setToZero();
		}
		return error;
	}

	if (error != 0) {
		res.setToZero();
	}

	return error;

	LOGMSG("\nError. Implementation. bigint::subAbs\n");
	exit(99);
}

BigInt& BigInt::operator=(const BigInt avalue) {
	if (this != &avalue) {
		vorz=avalue.vorz;
		highestusedidx=avalue.highestusedidx;
		for(int32_t i=0;i<=highestusedidx;i++) {
			digits[i]=avalue.digits[i];
		}
	}

	return *this;
}

void BigInt::ausgabe(FILE *f) {
	if (vorz==0) {
		fprintf(f,"+0");
		return;
	}

	if (vorz<0) fprintf(f,"-"); else fprintf(f,"+");

	for(int32_t i=highestusedidx;i>=0;i--) {
		if (i != highestusedidx) {
			if (BIGINTBASE == 1000000000) {
				fprintf(f,"%09i",digits[i]);
			} else {
				fprintf(f,"/%i/",digits[i]);
			}
		} else fprintf(f,"%i",digits[i]);
	}
}

#define CHECKBIGINTDIV \
{\
	BigInt tmp1,tmp2;\
	error += bigintMul_TAB(tmp1,dividing,B);\
	error += bigintAdd_TAB(tmp2,tmp1,remainder);\
	if (bigintVgl_AB(tmp2,A) != 0) {\
		LOGMSG("\nError. BigInt Div incorrect test\n");\
		exit(99);\
	}\
	\
	error += bigintAdd_abs_TAB(tmp2,tmp1,remainder);\
	if (bigintVgl_abs_AB(tmp2,A) > 0) {\
		LOGMSG("\nError. BigInt Div incorrect test/2\n");\
		exit(99);\
	}\
}

int32_t bigintDiv9_abs_TRAB(
	BigInt& dividing,
	BigInt& remainder,
	BigInt& A,
	BigInt& B
) {
	// compute |A|/|B|
	// it holds: 0 <= |A|/|B| < 10

	int32_t error=0;

	if (
		(A.vorz < 0) ||
		(B.vorz < 0)
	) {
		LOGMSG("\nError. bigintDiv9_abs_TRAB negative sign\n");
		exit(99);
	}

	if (A.vorz == 0) {
		dividing.setToZero();
		remainder.setToZero();
		return 0;
	}

	if (B.vorz == 0) {
		dividing.setToZero();
		remainder.setToZero();
		return -1; // error

		LOGMSG("\nError. bigintDiv9_abs_TRAB division by zero");
		exit(99);
	}

	// A,B are positive

	// naive implementation
	BigInt tmp;
	for(int32_t i=9;i>=0;i--) {
		error += bigintMul_digit_TDA(tmp,i,B);
		int32_t vgl=bigintVgl_abs_AB(tmp,A);
		if (vgl == 0) {
			dividing.set_int64(i);
			remainder.setToZero();
			#ifdef _INVOKECLAIMVERIFICATIONS
			CHECKBIGINTDIV
			#endif
			if (error != 0) {
				dividing.setToZero();
				remainder.setToZero();
			}
			return error;
		} else if (vgl < 0) {
			// as descending loop => i times
			dividing.set_int64(i);
			error += bigintSub_abs_TAB(remainder,A,tmp);
			#ifdef _INVOKECLAIMVERIFICATIONS
			CHECKBIGINTDIV
			#endif
			if (error != 0) {
				dividing.setToZero();
				remainder.setToZero();
			}
			return error;
		}
	} // i

	if (error != 0) {
		dividing.setToZero();
		remainder.setToZero();
	}

	return error;

	LOGMSG("\nError. Implementation. bigintDiv9_abs_TRAB. End of while loop\n");
	exit(99);
}


int32_t bigintDiv_TRAB(
	BigInt& dividing,
	BigInt& remainder,
	BigInt& A,
	BigInt& B
) {
	dividing.setToZero();
	remainder.setToZero();

	// division is performed in |A| / |B|
	// and sign is adjusted to fit the following

	int32_t error=0;

	// at return, it holds:
	//		dividing * B + remainder = A
	// and
	//		|dividing| * |B| + |remainder| <= |A|

	if (A.vorz == 0) {
		// it holds: 0
		dividing.set_int64(0);
		remainder.setToZero();
		#ifdef _INVOKECLAIMVERIFICATIONS
		CHECKBIGINTDIV
		#endif
		if (error != 0) {
			dividing.setToZero();
			remainder.setToZero();
		}
		return error;
	}

	if (B.vorz == 0) {
		LOGMSG("\nError. division by zero. BigInt.\n");
		exit(99);
	}

	if (
		(B.highestusedidx == 0) &&
		(B.digits[0] == 1)
	) {
		// it holds: dividing * 1 + remainder = A
		dividing.copyFrom(A);
		dividing.vorz *= B.vorz; // A / -1 = -A
		remainder.setToZero();
		#ifdef _INVOKECLAIMVERIFICATIONS
		CHECKBIGINTDIV
		#endif
		if (error != 0) {
			dividing.setToZero();
			remainder.setToZero();
		}
		return error;
	}

	int32_t vgl=bigintVgl_abs_AB(A,B);

	if (vgl == 0) {
		if (A.vorz == B.vorz) {
			// it holds: 1 * B + 0 = A
			remainder.setToZero();
			dividing.set_int64(1);
		} else if (A.vorz < 0) {
			// div * |B| + rem = -|A|
			remainder.setToZero();
			dividing.set_int64(-1);
		} else { // A > 0, B < 0
			// div * -|B| + rem = |A|
			remainder.setToZero();
			dividing.set_int64(-1);
		}

		#ifdef _INVOKECLAIMVERIFICATIONS
		CHECKBIGINTDIV
		#endif
		if (error != 0) {
			dividing.setToZero();
			remainder.setToZero();
		}

		return error;

	} // |A|=|B|

	if (vgl < 0) {
		// |A| < |B|
		// it holds: 0 * B + A = A no matter the sign
		remainder.copyFrom(A);
		dividing.set_int64(0);
		#ifdef _INVOKECLAIMVERIFICATIONS
		CHECKBIGINTDIV
		#endif
		if (error != 0) {
			dividing.setToZero();
			remainder.setToZero();
		}
		return error;
	}

	// |A| > |B|. Sign is not used during division

	BigInt rest;
	BigInt absA(A);
	BigInt absB(B);
	if (absA.vorz != 0) absA.vorz=1;
	if (absB.vorz != 0) absB.vorz=1;

	rest.copyFrom(absA);
	dividing.setToZero();

	int32_t bten=B.tendigitcount();

	while (1) {

		if (rest.vorz == 0) {
			// dividing is already correct
			remainder.setToZero();
			if (A.vorz > 0) {
				if (B.vorz > 0) {
					// all correct
				} else {
					// -|B|*dividing + 0 = |A|
					dividing.vorz *= -1;
				}
			} else {
				if (B.vorz > 0) {
					// |B|*dividing + 0 = -|A|
					dividing.vorz *= -1;
				} else {
					// -|B|*dividing + 0 = -|A|
				}
			} // A.vorz <= 0

			#ifdef _INVOKECLAIMVERIFICATIONS
			CHECKBIGINTDIV
			#endif
			if (error != 0) {
				dividing.setToZero();
				remainder.setToZero();
			}
			return error;
		}

		int32_t vgl=bigintVgl_abs_AB(rest,B);
		if (vgl < 0) {
			remainder.copyFrom(rest);

			if (A.vorz > 0) {
				if (B.vorz > 0) {
					// |B|*dividing + rem = |A|
					// all positive
				} else {
					// -|B|*dividing + rem = |A|
					dividing.vorz *= -1;
				}
			} else {
				if (B.vorz > 0) {
					// |B|*dividing + rem = -|A|
					dividing.vorz *= -1;
					remainder.vorz *= -1;
				} else {
					// -|B|*dividing + rem = -|A|
					remainder.vorz *= -1;
				}
			}

			#ifdef _INVOKECLAIMVERIFICATIONS
			CHECKBIGINTDIV
			#endif
			if (error != 0) {
				dividing.setToZero();
				remainder.setToZero();
			}
			return error;

		}

		// use a shifted version of |B|
		// compute |rest| / |shifted B|
		// is 0 <= .. < 10

		int32_t shiftdigits=rest.tendigitcount() - bten;
		BigInt Bshifted(absB);
		error += Bshifted.shiftLeft10(shiftdigits);

		// now Bshifted could be larger than rest
		// then shift one to the right
		if (bigintVgl_abs_AB(Bshifted,rest) > 0) {
			shiftdigits--;
			Bshifted.copyFrom(absB);
			error += Bshifted.shiftLeft10(shiftdigits);
		}

		if (shiftdigits < 0) {
			LOGMSG("\nError. bigintDiv/5 shifted too large\n");
			exit(99);
		}

		#ifdef _INVOKECLAIMVERIFICATIONS
		if (bigintVgl_abs_AB(Bshifted,A) > 0) {
			LOGMSG("\nError. bigintDiv/4 shifted too large\n");
			exit(99);
		}
		#endif

		BigInt div9,rem9;
		error += bigintDiv9_abs_TRAB(div9,rem9,rest,Bshifted);

		if (error != 0) return error;

		if (div9.vorz == 0) {
			// error
			LOGMSG("\nError. bigintDiv. B shifted divides 0 times\n");
			exit(99);
		}

		div9.shiftLeft10(shiftdigits);
		dividing.addTo(div9);
		rest.copyFrom(rem9);

	} // while

	if (error != 0) {
		dividing.setToZero();
		remainder.setToZero();
	}

	return error;

}


int32_t BigInt::tendigitcount(void) {
	// how many digits in base 10 would the current number
	// need (value relevant for division)

	if (BIGINTBASE != 1000000000) {
		LOGMSG("\nError. Implementation. Bigint::tendigitcount needs base to be 10^9\n");
		exit(99);
	}

	int32_t ct=highestusedidx * 9;

	for(int32_t t=8;t>=0;t--) {
		if (digits[highestusedidx] >= tenpower[t]) {
			ct = sum_int32t(ct,t+1);
			break;
		}
	}

	return ct;

}

int8_t BigInt::isZero(void) {
	if (vorz == 0) return 1;
	return 0;
}

void BigInt::invertSign(void) {
	vorz *= -1;
}

int32_t bigintPow_TAE(
	BigInt& res,
	BigInt& value,
	const int32_t power
) {
	int32_t error=0;
	res.setToZero();

	if (power < 0) {
		LOGMSG("\nError. pow bigInt not defined for negative exponents\n");
		exit(99);
	}

	res.set_int64(1);

	for(int32_t i=1;i<=power;i++) {
		BigInt w;
		error += bigintMul_TAB(w,res,value);
		res.copyFrom(w);
	} // i

	if (error != 0) {
		res.setToZero();
	}

	return error;

}

int8_t BigInt::isPositiveOne(void) {
	if (
		(vorz > 0) &&
		(highestusedidx == 0) &&
		(digits[0] == 1)
	) return 1;

	return 0;

}

int32_t bigintGcd_abs_TAB(
	BigInt& gcd,
	BigInt& A,
	BigInt& B
) {
    // ret != 0 => error, else success
	// sign not relevant
	// |A| > |B|, if not => return gcd(B,A) switch places

	int8_t vgl=bigintVgl_abs_AB(A,B);
	if (vgl == 0) {
		gcd.copyFrom(A);
		gcd.vorz=1;
		return 0; // success
	}

	if (vgl < 0) {
        return bigintGcd_abs_TAB(gcd,B,A);
	}

	BigInt zaehler,nenner;
	gcd.set_int64(0);
	zaehler.copyFrom(A);
	zaehler.vorz=1;
	nenner.copyFrom(B);
	nenner.vorz=1;
	BigInt lastrem;
	int8_t first=1;

	while (1) {
		BigInt div,rem;
		bigintDiv_TRAB(div,rem,zaehler,nenner);

		if (rem.vorz == 0) {
			if (first > 0) {
				gcd.copyFrom(B);
			} else {
				gcd.copyFrom(lastrem);
			}
			return 0;
		}

		first=0;

		lastrem.copyFrom(rem);
		zaehler.copyFrom(nenner);
		nenner.copyFrom(rem);

	} // while

	return 0;

}

#endif
