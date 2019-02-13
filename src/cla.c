#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <mpi.h>
#include <stdbool.h>

//EXAMPLE DATA STRUCTURE DESIGN AND LAYOUT FOR CLA
#define input_size 262144
#define block_size 32
//Do not touch these defines
#define digits (input_size+1)
#define bits input_size*4
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
int *subSumi = NULL;
//Integer array of inputs in binary form
unsigned int bin1[bits];
unsigned int bin2[bits];
//Character array of inputs in hex form
char hex1[digits];
char hex2[digits];
//character array to store the final hex answer
char hexAns[digits];
//MPI data
int numRanks = -1;
int rank = -1;
float rankFactor = -1;
int elementsPerProc = -1;
bool usingBarrier = true;
//new vars for storing data subsets
int *subBin1 = NULL;
int *subBin2 = NULL;
int *fullSumi = NULL;

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
	for (int i = 0; i < (int)(bits*rankFactor); ++i) {
		//NOTE: switched bin1 and bin2 to subBin1 and subBin2 so that each rank uses its scattered bits, rather than an empty array
		gi[i] = subBin1[i] & subBin2[i];
		pi[i] = (subBin1[i] | subBin2[i]);
	}
}

/**
 * step 2: calculate ggj and gpj, using gi and pi
 * ggj is the generate function for the current 8-bit group.
 * gpj is the propagate function for the current 8-bit group.
 */
void calc_ggj_gpj() {
	for (int i = 0; i < (int)(ngroups*rankFactor); ++i) {
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
	for (int i = 0; i < (int)(nsections*rankFactor); ++i) {
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
	for (int i = 0; i < (int)(nsupersections*rankFactor); ++i) {
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
	for (int i = 1; i < (int)(nsupersections*rankFactor); ++i) {
		sscl[i] = (ssgl[i] | (sspl[i]&sscl[i-1]));
	}
}
/**
 * step 6: calculate sck, using sgk and spk, as well as correct sscl
 * sck stores whether or not there is a carry bit for the current 64-bit section
 */
void calc_sck() {
	sck[0] = sgk[0];
	for (int i = 1; i < (int)(nsections*rankFactor); ++i) {
		sck[i] = (sgk[i] | (spk[i]&(i%block_size==0 ? sscl[i/block_size-1] : sck[i-1])));
	}
}
/**
 * step 7: calculate gcj, using ggj and gpj, as well as correct sck
 * gcj stores whether or not there is a carry bit for the current 8-bit group
 */
void calc_gcj() {
	gcj[0] = ggj[0];
	for (int i = 1; i < (int)(ngroups*rankFactor); ++i) {
		gcj[i] = (ggj[i] | (gpj[i]&(i%block_size==0 ? sck[i/block_size-1] : gcj[i-1])));
	}
}
/**
 * step 8: calculate ci, using gi and pi, as well as correct gcj
 * ci stores whether or not there is a carry bit for the current bit
 */
void calc_ci() {
	ci[0] = gi[0];
	for (int i = 1; i < (int)(bits*rankFactor); ++i) {
		ci[i] = (gi[i] | (pi[i]&(i%block_size==0 ? gcj[i/block_size-1] : ci[i-1])));
	}
}

/**
 * step 9: calculate sumi, using bin1 and bin2, as well as ci
 * sumi stores the final sum for the current bit
 */
void calc_sumi() {
	//NOTE: changed bin1 and bin2 to subBin1 and subBin2 here as well
	subSumi[0] = subBin1[0] ^ subBin2[0];
	for (int i = 1; i < (int)(bits*rankFactor); ++i) {
		subSumi[i] = subBin1[i] ^ subBin2[i] ^ ci[i-1];
	}
}

/**
 * scatter the input binary so that everyone gets their own chunk
 */
void scatterData() {
	subBin1 = malloc(sizeof(int) * elementsPerProc);
	subBin2 = malloc(sizeof(int) * elementsPerProc);

	MPI_Scatter(bin1, elementsPerProc, MPI_INT, subBin1, elementsPerProc, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(bin2, elementsPerProc, MPI_INT, subBin2, elementsPerProc, MPI_INT, 0, MPI_COMM_WORLD);
}

/**
 * master collects the calculated output binary
 */
void gatherData() {
	if (rank == 0) {
		fullSumi = malloc(sizeof(int) * bits);
	}
	MPI_Gather(subSumi, 1, MPI_INT, fullSumi, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

/**
 * main cla routine; calls all of the different sub-steps in order, with barriers between each step if enabled
 */
void cla() {
	scatterData();
	if (usingBarrier) MPI_Barrier(MPI_COMM_WORLD);
	calc_gi_pi();
	if (usingBarrier) MPI_Barrier(MPI_COMM_WORLD);
	calc_ggj_gpj();
	if (usingBarrier) MPI_Barrier(MPI_COMM_WORLD);
	calc_sgk_spk();
	if (usingBarrier) MPI_Barrier(MPI_COMM_WORLD);
	calc_ssgl_sspl();
	if (usingBarrier) MPI_Barrier(MPI_COMM_WORLD);
	calc_sscl();
	if (usingBarrier) MPI_Barrier(MPI_COMM_WORLD);
	calc_sck();
	if (usingBarrier) MPI_Barrier(MPI_COMM_WORLD);
	calc_gcj();
	if (usingBarrier) MPI_Barrier(MPI_COMM_WORLD);
	calc_ci();
	if (usingBarrier) MPI_Barrier(MPI_COMM_WORLD);
	calc_sumi();
	if (usingBarrier) MPI_Barrier(MPI_COMM_WORLD);
	gatherData();
}

/**
 * convert the calculated sumi back into hex for final solution
 */
void convertAnswerToHex() {
	for (int i = 0; i < input_size; hexAns[i] = binToHex(fullSumi+4*(input_size-i-1)), ++i);
	printf("%s\n",hexAns);
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

int main(int argc, char* argv[]) {
	//usage check
	if (argc != 3) {
		fprintf(stderr,"Error: %d input argument[s] supplied, but 2 were expected. Usage: mpirun -np X ./cla test1.txt test1-output.txt",argc-1);
		exit(1);
	}

	//init MPI and grab size and rank
	MPI_Init( &argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	rankFactor = 1/(float)numRanks;
	elementsPerProc = bits*rankFactor;
	subSumi = malloc(sizeof(int) * elementsPerProc);
	printf("my rank is %d | total size is %d | rank factor is 1\\%d = %f | elementsPerProc is %d\n",rank,numRanks,numRanks,rankFactor,elementsPerProc);

	//treat test1.txt as stdin, and test1-output.txt as stdout
	if (rank == 0) {
		freopen(argv[1], "r", stdin);
		freopen(argv[2], "w", stdout);
		//master reads and formats the input
		parseInput();
	}

	//everybody runs the main routine
	cla();

	//master handles the final output
	if (rank == 0) {
		convertAnswerToHex();
		free(fullSumi);
	}
	free(subBin1);
	free(subBin2);
	free(subSumi);
	MPI_Finalize();
}
