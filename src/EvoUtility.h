///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#define NLAYER3 3 //number of output nodes

#define MLAYER1 75 //the second network
#define MLAYER2 40
#define MLAYER3 3

#define ALPHA 0.001 //learning rate
#define LAMDA 0.9 //momentum: keep last weight
#define MAX_LEN 1000//4000
#define LenScale 1000

#define	GAP_OPEN_PENALTY -11.0 // gap open penalty
#define GAP_XTN_PENALTY -11.0 // -1.0 // gap extension penalty (not added to gapo)

namespace utility {

class Blosum62 {
	public:
		/* Convert AA letter to numeric code (0-22) */
                int aanum(int ch);
        	/*  BLOSUM 62 */
                short aamat(int i,int j);
	};

}

typedef struct _SequenceData {
		int len;
		char *seq;
		char *ss1,*ss2,*sa1,*sa2;
}SequenceData;


namespace Text {
  int ncharLine(FILE* fp);
  long int nnewline(FILE* fp);
}


#endif

