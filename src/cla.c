#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

//EXAMPLE DATA STRUCTURE DESIGN AND LAYOUT FOR CLA
#define input_size 1024
#define block_size 8
//Do not touch these defines
#define digits (input_size+1)
#define bits digits*4
#define ngroups bits/block_size
#define nsections ngroups/block_size
#define nsupersections nsections/block_size

//Global definitions of the various arrays used in steps for easy access
int gi[bits] = {0};
int pi[bits] = {0};
int ci[bits] = {0};
int ggj[ngroups] = {0};
int gpj[ngroups] = {0};
int gcj[ngroups] = {0};
int sgk[nsections] = {0};
int spk[nsections] = {0};
int sck[nsections] = {0};
int ssgl[nsupersections] = {0} ;
int sspl[nsupersections] = {0} ;
int sscl[nsupersections] = {0} ;
int sumi[bits] = {0};
//Integer array of inputs in binary form
unsigned int bin1[bits];
unsigned int bin2[bits];
//Character array of inputs in hex form
char hex1[digits];
char hex2[digits];
//character array to store the final hex answer
char hexAns[digits];

/**
 * simple hex char to binary conversion
 * @param hex: a single hex value in character form (0-F)
 * @param bin: an int array pointer at which to store the 4 calculated binary values
 */
void hexToBin(char hex, unsigned int* bin) {
	int decimal = (9*!!(hex&64))+(hex&15);
	*bin = decimal&1;
	*(bin+1) = decimal>>1&1;
	*(bin+2) = decimal>>2&1;
	*(bin+3) = decimal>>3&1;
}

/**
 * simple binary to hex char conversion
 * @param bin: an int array pointer at which the 4 calculated binary values are stored
 */
char binToHex(int* bin) {
	int decimalVal = *bin + 2**(bin+1) + 4**(bin+2) + 8**(bin+3);
	return decimalVal < 10 ? decimalVal + '0' : 55+decimalVal;
}

/**
 * step 1: calculate gi and pi, using bin1 and bin2
 * gi is the generate function. It tells us if we generated a carry in the ith stage (occurs when a[i] and b[i] are both 1) for the current bit.
 * pi is the propagate function. It tells us if we propagated a carry in the ith stage for the current bit.
 */
void calc_gi_pi() {
	for (int i = 0; i < bits; ++i) {
		gi[i] = bin1[i] & bin2[i];
		pi[i] = (bin1[i] | bin2[i]);
	}
}

/**
 * step 2: calculate ggj and gpj, using gi and pi
 * ggj is the generate function for the current 8-bit group.
 * gpj is the propagate function for the current 8-bit group.
 */
void calc_ggj_gpj() {
	for (int i = 0; i < ngroups; ++i) {
		int iblock = block_size*i;
		ggj[i] = gi[iblock+block_size-1];
		for (int j = iblock+block_size-2; j >= iblock; --j) {
			int curLine = pi[iblock+block_size-1] & gi[j];
			for (int k = iblock+block_size-2; k > j; curLine &= pi[k--]);
			ggj[i] |= curLine;
		}

		gpj[i] = pi[iblock+block_size-1];
		for (int j = block_size-2; j >= 0; gpj[i] &= pi[iblock+j--]);
	}
}

/**
 * step 3: calculate sgk and spk, using ggj and gpj
 * sgk is the generate function for the current 64-bit section.
 * spk is the propagate function for the current 64-bit section.
 */
void calc_sgk_spk() {
	for (int i = 0; i < nsections; ++i) {
		int iblock = block_size*i;
		sgk[i] = ggj[iblock+block_size-1];
		for (int j = iblock+block_size-2; j >= iblock; --j) {
			int curLine = gpj[iblock+block_size-1] & ggj[j];
			for (int k = iblock+block_size-2; k > j; curLine &= gpj[k--]);
			sgk[i] |= curLine;
		}

		spk[i] = gpj[iblock+block_size-1];
		for (int j = block_size-2; j >= 0; spk[i] &= gpj[iblock+j--]);
	}
}

/**
 * step 4: calculate ssgl and sspl, using sgk and spk
 * ssgl is the generate function for the current 512-bit super section.
 * sspl is the propagate function for the current 512-bit super section.
 */
void calc_ssgl_sspl() {
	for (int i = 0; i < nsupersections; ++i) {
		int iblock = block_size*i;
		ssgl[i] = sgk[iblock+block_size-1];
		for (int j = iblock+block_size-2; j >= iblock; --j) {
			int curLine = spk[iblock+block_size-1] & sgk[j];
			for (int k = iblock+block_size-2; k > j; curLine &= spk[k--]);
			ssgl[i] |= curLine;
		}

		sspl[i] = spk[iblock+block_size-1];
		for (int j = block_size-2; j >= 0; sspl[i] &= spk[iblock+j--]);
	}
}

/**
 * step 5: calculate sscl, using ssgl and sspl, and 0 for ssc-1
 * sscl stores whether or not there is a carry bit for the current 512-bit super section.
 */
void calc_sscl() {
	sscl[0] = ssgl[0];
	for (int i = 1; i < nsupersections; ++i) {
		sscl[i] = (ssgl[i] | (sspl[i]&sscl[i-1]));
	}
}
/**
 * step 6: calculate sck, using sgk and spk, as well as correct sscl
 * sck stores whether or not there is a carry bit for the current 64-bit section
 */
void calc_sck() {
	sck[0] = sgk[0];
	for (int i = 1; i < nsections; ++i) {
		sck[i] = (sgk[i] | (spk[i]&(i%block_size==0 ? sscl[i/block_size-1] : sck[i-1])));
	}
}
/**
 * step 7: calculate gcj, using ggj and gpj, as well as correct sck
 * gcj stores whether or not there is a carry bit for the current 8-bit group
 */
void calc_gcj() {
	gcj[0] = ggj[0];
	for (int i = 1; i < ngroups; ++i) {
		gcj[i] = (ggj[i] | (gpj[i]&(i%block_size==0 ? sck[i/block_size-1] : gcj[i-1])));
	}
}
/**
 * step 8: calculate ci, using gi and pi, as well as correct gcj
 * ci stores whether or not there is a carry bit for the current bit
 */
void calc_ci() {
	ci[0] = gi[0];
	for (int i = 1; i < bits; ++i) {
		ci[i] = (gi[i] | (pi[i]&(i%block_size==0 ? gcj[i/block_size-1] : ci[i-1])));
	}
}

/**
 * step 9: calculate sumi, using bin1 and bin2, as well as ci
 * sumi stores the final sum for the current bit
 */
void calc_sumi() {
	sumi[0] = bin1[0] ^ bin2[0];
	for (int i = 1; i < bits; ++i) {
		sumi[i] = bin1[i] ^ bin2[i] ^ ci[i-1];
	}
}

/**
 * main cla routine; calls all of the different sub-steps in order
 */
void cla() {
	calc_gi_pi();
	calc_ggj_gpj();
	calc_sgk_spk();
	calc_ssgl_sspl();
	calc_sscl();
	calc_sck();
	calc_gcj();
	calc_ci();
	calc_sumi();
}

/**
 * convert the calculated sumi back into hex for final solution
 */
void convertAnswerToHex() {
	for (int i = 0; i < input_size; hexAns[i] = binToHex(sumi+4*(input_size-i-1)), ++i);
	printf("0%s\n",hexAns);
}

/**
 * parse the piped input, reading in two 1024 hex strings and converting them to 4096 byte arrays
 */
void parseInput() {
	//first hex string
	scanf("%s", hex1);
	for (int i = 0; i < input_size; hexToBin(hex1[i],bin1+4*(input_size-i-1)), ++i);
	//second hex string
	scanf("%s", hex2);
	for (int i = 0; i < input_size; hexToBin(hex2[i],bin2+4*(input_size-i-1)), ++i);
}

int main(void) {
	parseInput();
	cla();
	convertAnswerToHex();
}
