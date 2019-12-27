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

/******************************************************************************/
/*                                                                            */
/*                                 SeqAlign.h                                 */
/*                                                                            */
/*    A class implementing the alignment between two character sequences.     */
/*                                                                            */
/******************************************************************************/

#ifndef SEQ_ALIGN_H
#define SEQ_ALIGN_H

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "EvoAminoName.h"
#include "EvoUtility.h"


namespace CharSeq {

	class SeqAlign {

		int len;
    //SequenceData seqData;
		float** prf;
		float weight_ss, weight_sa, weight_angle;
		char seq[MAX_LEN],secstr1[MAX_LEN],secstr2[MAX_LEN],solacc1[MAX_LEN],solacc2[MAX_LEN];
		float *phiT,*psiT,*phiQ,*psiQ;

		enum Denum {D1E, D2E, D3E}; // identifiers for the alignment matrices

		// type for any item in an alignment matrix
		struct MtxItem {
			float score;  // (optimal) alignment score
			bool pvmtx[3]; // tells which matrices provided the optimal previous-step sub-score

			MtxItem() {};
			MtxItem(const float s, const bool pv1, const bool pv2, const bool pv3);
			MtxItem(const float max, const float d1, const float d2, const float d3);
		};

    // DX[i][j].score = optimal align. score of s1[0],...,s1[i] and s2[0],...,s2[j]
		MtxItem** D1; // D1[i][j].score = score of the alignment ending with (s1[i], -)
		MtxItem** D2; // D2[i][j].score = score of the alignment ending with (-, s2[j])
		MtxItem** D3; // D3[i][j].score = score of the alignment ending with (s1[i], s2[j])

		MtxItem** getMtx() const;

		void updateD1(const int i, const int j);
		void updateD2(const int i, const int j);
		void updateD3(const int i, const int j);

		void readMkPrf(char prffile[500]);
		void readPhiPsi(char phipsipath[500]);
		float Score(const int i1, const int i2);

	public:
		SeqAlign(SequenceData dsInfo,char prffile[500],char path[500],float wss, float wsa, float wang);
		~SeqAlign();
		float printMtxInfo() const;
  };
}

#endif
