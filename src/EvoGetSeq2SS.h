///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#ifndef _getSeq2SS_h_
#define _getSeq2SS_h_

#include "EvoUtility.h"

namespace SSPrediction {

  class getSeq2SS {

  /* few constant values */
  /*  #define NLAYER3 3 //number of output nodes

  #define MLAYER1 75 //the second network
  #define MLAYER2 40
  #define MLAYER3 3
  
  #define ALPHA 0.001 //learning rate
  #define LAMDA 0.9 //momentum: keep last weight
  #define MAX_LEN 4000
  #define LenScale 1000 */

  /* variable declaration */
  int Len,SEG,NLAYER1,NLAYER2;
  char StrA1[MAX_LEN],StrA2[MAX_LEN];

  float Score1[MAX_LEN][21],Score[MAX_LEN][4];
  float nout[MAX_LEN][NLAYER3],tgtout[NLAYER3],mout[MLAYER3];

  float im1[MLAYER1],im2[MLAYER2], mout2[MLAYER2],im3[MLAYER3];
  float in1[525],in2[150],out2[150],in3[NLAYER3],out[NLAYER3];

  float v1[MLAYER1][MLAYER2],Dv1[MLAYER1][MLAYER2];
  float v2[MLAYER2][MLAYER3],Dv2[MLAYER2][MLAYER3];
  float w1[525][150],Dw1[525][150];
  float w2[150][NLAYER3],Dw2[150][NLAYER3];

  char wgtfile[500],propfile[500];

/* function declaration */
  public:
    getSeq2SS(char prdA2[MAX_LEN],char fasta[MAX_LEN],char path[500]);
    int Read_wgt();
    int Read_score1(char fasta[MAX_LEN],float divisor); 
    int Read_score();
    void FeedNet(int I);
    void Netout();
    void FeedNet1(int I);
    void Netout1();
  };
}

#endif
