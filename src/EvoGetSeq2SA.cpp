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

#pragma warning(disable:4996)
#pragma warning(disable:4305)
#pragma warning(disable:4244)
#include "EvoGetSeq2SA.h"


SAPrediction::getSeq2SA::getSeq2SA(char prdA2[MAX_LEN],char fasta[MAX_LEN],char sstruct[MAX_LEN],char evolutionpath[500]){
  int i,j;
    
  SEG=12;
  NLAYER1= SEG*53; //number of input nodes=SEG*53
  NLAYER2= 120;    //number of  hidden nodes

  sprintf(wgtfile,"%s/param/wgtSA",evolutionpath);
  sprintf(propfile,"%s/param/SApropensity.txt",evolutionpath);

  Read_wgt();
  Read_score1(fasta);    // Score1 is from Blosum62 matrix, normalized by 2.0 and taken as algebric sigmoid function
  Read_score(sstruct);   // Score is directly taken from propensity file. Score2 is composition mutiplied by Sander A-X-A. Score 3 is SS info

  for(i=1;i<=Len;i++){
    FeedNet(i);
    Netout(); 
    for(j=0;j<NLAYER3;j++){
      nout[i][j]=out[j];
    }
  }
    
  for(i=1;i<=Len;i++){
    FeedNet1(i);
    Netout1();
    if(mout[0]>=mout[1] && mout[0]>=mout[2]){
      prdA2[i]='1';          // prdA2[i]='B'; 
    }
    else if(mout[1]>=mout[0] && mout[1]>=mout[2]){
      prdA2[i]='2';          // prdA2[i]='I';  
    }
    else if(mout[2]>=mout[0] && mout[2]>=mout[1]) {
      prdA2[i]='3';          // prdA2[i]='E'; 
    }
  } 

  for(i=0;i<Len;i++){
    prdA2[i]=prdA2[i+1];
  }
}


//ABCDEFGHIKLMNPQRSTVWXYZ
inline float SAPrediction::getSeq2SA::SANDER(int i){
  float sANDER[23]={106,160,135,163,194,197,84,184,169,205,164,188,157,136,198,248,130,142,142,227,180,222,196};
  return sANDER[i];
}


void SAPrediction::getSeq2SA::Read_wgt(){
  FILE *fpout;
  int i,j;

  fpout=fopen(wgtfile,"r");     //weight file
  for(i=0;i<NLAYER1;i++){       //330
    for(j=0;j<NLAYER2;j++){     //150
      fscanf(fpout,"%f\n",&w1[i][j]);
      Dw1[i][j]=0.0;
    }
  }
  for(i=0;i<NLAYER2;i++){       //150
    for(j=0;j<NLAYER3;j++){     //3
      fscanf(fpout,"%f\n",&w2[i][j]);
      Dw2[i][j]=0.0;
    }
  }
  for(i=0;i<MLAYER1;i++){       //75
    for(j=0;j<MLAYER2;j++){     //40
      fscanf(fpout,"%f\n",&v1[i][j]);
      Dv1[i][j]=0.0;
    }
  }
  for(i=0;i<MLAYER2;i++){       //40
    for(j=0;j<MLAYER3;j++){     //3
      fscanf(fpout,"%f\n",&v2[i][j]);
      Dv2[i][j]=0.0;
    }
  }
  fclose(fpout);
}

/****************************************************************/
int SAPrediction::getSeq2SA::Read_score(char sstruct[MAX_LEN]){
  int i,j,k;
  FILE *fp=NULL;
  char line[5000]="XXX",amn='X',amino[MAX_LEN]="XXX";
  float tmp1=0,tmp2=0,tmp3=0,aa[24]={0};
  float propensity[MAX_LEN][4]={0};
  
  fp=fopen(propfile,"r");
  fgets(line,5000,fp);
  i=0;
  while(!feof(fp)){
    fgets(line,5000,fp);
    if(feof(fp)) break;
    sscanf(line,"%c %f %f %f",&amn,&tmp1,&tmp2,&tmp3);
    amino[i]=amn;
    propensity[i][1]=tmp1;    // B
    propensity[i][2]=tmp2;    // I 
    propensity[i][3]=tmp3;    // E
    i++;
  }
  fclose(fp);
  for(i=1;i<=Len;i++){
    for(j=0;j<20;j++) if(StrA1[i]==amino[j]) break;
    for(k=0;k<3;k++){
      tmp1=propensity[j][k+1];
      Score[i][k]=tmp1/sqrt(1.0+tmp1*tmp1);
    }
  }

  for(i=1;i<=Len;i++){
    for(j=0;j<23;j++) aa[j]=0.0;
    for(j=(i-SEG/2);j<(i-SEG/2)+SEG;j++) {
      if(j>0 && j<Len)
        aa[getID(StrA1[j])]+=SANDER(getID(StrA1[j]));
    }

    for(k=0;k<23;k++){
      tmp1=aa[k]/150.0;
      Score2[i][k]=tmp1/sqrt(1.0+tmp1*tmp1);
    }

    Score3[i][0]=Score3[i][1]=Score3[i][2]=0.0;
    if(sstruct[i-1]=='1')  Score3[i][0]=1.0;
    if(sstruct[i-1]=='2')  Score3[i][1]=1.0;
    if(sstruct[i-1]=='3')  Score3[i][2]=1.0;
  }  
  
  return(0);
}

/****************************************************************/
int SAPrediction::getSeq2SA::Read_score1(char fasta[MAX_LEN]){
  int i,j,aa_num=0;
  char *ncbicodes = "XAXCDEFGHIKLMNPQRSTVWXYXXX";
  float tmp3=0;
  utility::Blosum62 bsm;

  // read Score1[][]---------------->
  Len=strlen(fasta);
  for(i=1;i<=Len;i++){
    StrA1[i]=fasta[i-1];
  }
  for(i=1;i<=Len;i++){
    aa_num=bsm.aanum(StrA1[i]);
    Score1[i][0]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[1]));
    Score1[i][4]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[3]));
    Score1[i][3]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[4]));
    Score1[i][6]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[5]));
    Score1[i][13]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[6]));
    Score1[i][7]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[7]));
    Score1[i][8]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[8]));
    Score1[i][9]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[9]));
    Score1[i][11]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[10]));
    Score1[i][10]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[11]));
    Score1[i][12]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[12]));
    Score1[i][2]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[13]));
    Score1[i][14]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[14]));
    Score1[i][5]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[15]));
    Score1[i][1]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[16]));
    Score1[i][15]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[17]));
    Score1[i][16]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[18]));
    Score1[i][19]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[19]));
    Score1[i][17]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[20]));
    Score1[i][18]=bsm.aamat(aa_num,bsm.aanum(ncbicodes[22]));

    for(j=0;j<20;j++) {
      tmp3=Score1[i][j];
      Score1[i][j]=tmp3/sqrt(1.0+tmp3*tmp3);
    }
    Score1[i][20]=0;
  }

  return Len;
}


int SAPrediction::getSeq2SA::getID(char c){
  switch (c){
    case 'A': return 0;
    case 'B': return 1;
    case 'C': return 2;
    case 'D': return 3;
    case 'E': return 4;
    case 'F': return 5;
    case 'G': return 6;
    case 'H': return 7;
    case 'I': return 8;
    case 'K': return 9;
    case 'L': return 10;
    case 'M': return 11;
    case 'N': return 12;
    case 'P': return 13;
    case 'Q': return 14;
    case 'R': return 15;
    case 'S': return 16;
    case 'T': return 17;
    case 'V': return 18;
    case 'W': return 19;
    case 'X': return 20;
    case 'Y': return 21;
    case 'Z': return 22;
  };
  return -1;
}


/********************************************************************/
void SAPrediction::getSeq2SA::FeedNet(int I){
  int i,j,k;
  
  k=0;
  for(i=(I-SEG/2);i<(I-SEG/2)+SEG;i++){
    //---Score1-------->
    if(i<1 || i>Len){
      for(j=0;j<20;j++){
        in1[k]=0.0;
        k++;
      }
      in1[k]=(float)Len/LenScale;
      k++;
    }
    else{
      for(j=0;j<20;j++){
        in1[k]=Score1[i][j];
        k++;
      }
      in1[k]=(float)Len/LenScale;
      k++;
    }
    //---Score-------->
    if(i<1 || i>Len){
      for(j=0;j<3;j++){
        in1[k]=0.0;
        k++;
      }
      in1[k]=(float)Len/LenScale;
      k++;
    }
    else{
      for(j=0;j<3;j++){
        in1[k]=Score[i][j];
        k++;
      }
      in1[k]=(float)Len/LenScale;
      k++;
    }
    //---Score2-------->
    if(i<1 || i>Len){
      for(j=0;j<23;j++){
        in1[k]=0.0;
        k++;
      }
      in1[k]=(float)Len/LenScale;
      k++;
    }
    else{
      for(j=0;j<23;j++) {
        in1[k]=Score2[i][j];
        k++;
      }
      in1[k]=(float)Len/LenScale;
      k++;
    }
    //---Score3-------->
    if(i<1 || i>Len){
      for(j=0;j<3;j++){
        in1[k]=0.0;
        k++;
      }
      in1[k]=(float)Len/LenScale;
      k++;
    }
    else{
      for(j=0;j<3;j++){
        in1[k]=Score3[i][j];
        k++;
      }
      in1[k]=(float)Len/LenScale;
      k++;
    }
  }
  
  for(i=k;i<=NLAYER1;i++)
    in1[i]=0;
  
  for(i=0;i<NLAYER3;i++) 
    tgtout[i]=0.0;

  tgtout[(StrA2[I]-48)]=1.0; //StrA2[I]=48,49,50, Buried, Intermediate, Exposed (BIE)
}

/**************************************************************************/
void SAPrediction::getSeq2SA::Netout(){
  int i,j;
  
  for(i=0;i<NLAYER2;i++){
    in2[i]=0.0;
    for(j=0;j<NLAYER1;j++) 
      in2[i]+=w1[j][i]*in1[j];
    out2[i]=1.0/(1.0+exp(-1.0*in2[i]));
  }
  
  for(i=0;i<NLAYER3;i++){
    in3[i]=0.0;
    for(j=0;j<NLAYER2;j++)
      in3[i]+=w2[j][i]*out2[j];
    out[i]=1.0/(1.0+exp(-1.0*in3[i]));
  }
}

/***************************************************************/
void SAPrediction::getSeq2SA::FeedNet1(int I){
  int i,j,k;
  k=0;
  for(i=(I-SEG/2);i<(I-SEG/2)+SEG;i++){
    if(i<1 || i>Len){
      for(j=0;j<NLAYER3;j++){
        im1[k]=0.0;
        k++;
      }
      im1[k]=1.0;
      k++;
      im1[k]=(float)Len/LenScale;
      k++;
    }
    else{
      for(j=0;j<NLAYER3;j++){
        im1[k]=nout[i][j];
        k++;
      }
      im1[k]=0.0;
      k++;
      im1[k]=(float)Len/LenScale;
      k++;
    }
  } 
  for(i=0;i<MLAYER3;i++){
    tgtout[i]=0.0;
  }
  tgtout[(StrA2[I]-48)]=1.0; 
}

/**************************************************************************/
void SAPrediction::getSeq2SA::Netout1(){
  int i,j;
  for(i=0;i<MLAYER2;i++){
    im2[i]=0.0;
    for(j=0;j<MLAYER1;j++) 
      im2[i]+=v1[j][i]*im1[j];
    mout2[i]=1.0/(1.0+exp(-1.0*im2[i]));
  }

  for(i=0;i<MLAYER3;i++){
    im3[i]=0.0;
    for(j=0;j<MLAYER2;j++) 
      im3[i]+=v2[j][i]*mout2[j];
    mout[i]=1.0/(1.0+exp(-1.0*im3[i]));
  }
}

