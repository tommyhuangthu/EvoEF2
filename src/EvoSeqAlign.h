///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

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
