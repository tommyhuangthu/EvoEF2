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
