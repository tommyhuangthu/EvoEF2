/*******************************************************************************************************************************
This file is a part of the EvoDesign physical Energy Function (EvoEF)

Copyright (c) 2019 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/

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

