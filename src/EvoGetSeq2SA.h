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

#ifndef _getSeq2SA_h_
#define _getSeq2SA_h_

#include "EvoUtility.h"


namespace SAPrediction {
  class getSeq2SA {

  /* few constant values */
  //  #define SEG 12 //window size.
  
  //  #define NLAYER1 636 //number of input nodes=SEG*53
  //  #define NLAYER2 120  //number of  hidden nodes

  /* variable declaration */
  int Len,SEG,NLAYER1,NLAYER2;

  char StrA1[MAX_LEN],StrA2[MAX_LEN];
  float Score1[MAX_LEN][21],Score2[MAX_LEN][23],Score[MAX_LEN][4],Score3[MAX_LEN][4];
  float nout[MAX_LEN][NLAYER3],pout[MAX_LEN][NLAYER3];
  float w1[636][120],Dw1[636][120];
  float in1[636];
  float w2[120][NLAYER3],Dw2[120][NLAYER3];
  float in2[120],out2[120],in3[NLAYER3],out[NLAYER3];
  float tgtout[NLAYER3];
  float v1[MLAYER1][MLAYER2],Dv1[MLAYER1][MLAYER2];
  float v2[MLAYER2][MLAYER3],Dv2[MLAYER2][MLAYER3];
  float im1[MLAYER1],im2[MLAYER2], mout2[MLAYER2],im3[MLAYER3],mout[MLAYER3];

  float dw2[120][NLAYER3];
  float dv1[MLAYER1][MLAYER2],dv2[MLAYER2][MLAYER3]; // Derivatives

  char wgtfile[500],propfile[500];
 
  public:
  /* function declaration */
    getSeq2SA(char prdA2[MAX_LEN],char fasta[MAX_LEN],char sstruct[MAX_LEN],char path[500]);
    void Read_wgt();
    int Read_score1(char fasta[MAX_LEN]); 
    int Read_score(char sstruct[MAX_LEN]); 
    void FeedNet(int I);
    void Netout();
    void FeedNet1(int I);
    void Netout1();
    float SANDER(int i);
    int getID(char c);
  };
}

/*
  #define SEG 12 //window size.

  #define NLAYER1 636 //number of input nodes=SEG*53
  #define NLAYER2 120  //number of  hidden nodes
  #define NLAYER3 3 //number of output nodes

  #define MLAYER1 75 //the second network
  #define MLAYER2 40
  #define MLAYER3 3

  #define ALPHA 0.001 //learning rate
  #define LAMDA 0.9 //momentum: keep last weight
  #define REPEAT 200 //200
  #define MAX_LEN 4000
  #define LenScale 1000

  const int SEG,NLAYER1,NLAYER2,NLAYER3,MLAYER1,MLAYER2,MLAYER3,MAX_LEN,LenScale;
  const float ALPHA,LAMDA;

  char StrA1[MAX_LEN],StrA2[MAX_LEN];
  float Score1[MAX_LEN][21],Score2[MAX_LEN][23],Score[MAX_LEN][4],Score3[MAX_LEN][4];
  float nout[MAX_LEN][NLAYER3],pout[MAX_LEN][NLAYER3];
  int Len;
  float w1[NLAYER1][NLAYER2],Dw1[NLAYER1][NLAYER2];
  float in1[NLAYER1];
  float w2[NLAYER2][NLAYER3],Dw2[NLAYER2][NLAYER3];
  float in2[NLAYER2],out2[NLAYER2],in3[NLAYER3],out[NLAYER3];
  float tgtout[NLAYER3];
  float v1[MLAYER1][MLAYER2],Dv1[MLAYER1][MLAYER2];
  float v2[MLAYER2][MLAYER3],Dv2[MLAYER2][MLAYER3];
  float im1[MLAYER1],im2[MLAYER2], mout2[MLAYER2],im3[MLAYER3],mout[MLAYER3];

  float dw2[NLAYER2][NLAYER3];
  float dv1[MLAYER1][MLAYER2],dv2[MLAYER2][MLAYER3]; // Derivatives
*/

#endif
