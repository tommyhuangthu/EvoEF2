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
#include "EvoGetSeq2SS.h"

SSPrediction::getSeq2SS::getSeq2SS(char prdA2[MAX_LEN],char fasta[MAX_LEN],char evolutionpath[500]){
  /* variable declaration */
  FILE *fp=NULL;
  int i,j;
  float pout[MAX_LEN][NLAYER3]={0};
  char filename[500]="";
  
  SEG=21;
  NLAYER1=SEG*25;
  NLAYER2=150;

  Len=Read_score1(fasta,2.0);      // Score1 is from Blosum62 matrix, normalized by 2.0 and taken as algebric sigmoid function
  
  sprintf(wgtfile,"%s/param/wgtSS1",evolutionpath);
  //printf("XX%s\n", wgtfile);

  Read_wgt(); 
  sprintf(propfile,"%s/param/SSpropensity.txt",evolutionpath);
  Read_score();          // Score is directly taken from propensity file. It is used by next run also.

  for(i=1;i<=Len;i++) {
    FeedNet(i);
    Netout(); 
    for(j=0;j<NLAYER3;j++) {
      nout[i][j]=out[j];
    }
  }
    
  for(i=1;i<=Len;i++) {
    FeedNet1(i);
    Netout1();
    
    pout[i][0]=mout[0];          // 0 > H ; 1 > E ; 2 > C
    if(mout[1]>mout[0] && mout[1]>mout[2]) {
      pout[i][1]=mout[1]*1.3;    // 0 > H ; 1 > E ; 2 > C
    }
    if(mout[2]>mout[0] && mout[2]>mout[1]) {
      pout[i][2]=mout[2]*1.15;  // 0 > H ; 1 > E ; 2 > C
    }
  } 

  SEG=17;
  NLAYER1=SEG*25;
  Len=Read_score1(fasta,1.0);      // Score1 is from Blosum62 matrix, normalized by 1.0 and taken as algebric sigmoid function
  sprintf(wgtfile,"%s/param/wgtSS2",evolutionpath);
  Read_wgt(); 

  for(i=1;i<=Len;i++) {
    FeedNet(i);
    Netout(); 
    for(j=0;j<NLAYER3;j++) {
      nout[i][j]=out[j];
    }
  }
    
  for(i=1;i<=Len;i++) {
    FeedNet1(i);
    Netout1();

    pout[i][0]+=mout[0];        // 0 > H ; 1 > E ; 2 > C
    if(mout[1]>mout[0] && mout[1]>mout[2]) {
      pout[i][1]+=mout[1]*1.3;    // 0 > H ; 1 > E ; 2 > C
    }
    if(mout[2]>mout[0] && mout[2]>mout[1]) {
      pout[i][2]+=mout[2]*1.15;    // 0 > H ; 1 > E ; 2 > C
    }
  } 

  for(i=1;i<=Len;i++) {
    if(pout[i][0]>=pout[i][1] && pout[i][0]>=pout[i][2]) {
      prdA2[i]='1';            // H
    }
    else if(pout[i][1]>=pout[i][0] && pout[i][1]>=pout[i][2]) {
      prdA2[i]='2';            // E
    }
    else if(pout[i][2]>=pout[i][0] && pout[i][2]>=pout[i][1]) {
      prdA2[i]='3';            // C
    }
  }

  //smooth the predictions 
  for(i=1;i<=Len;i++) {  // EEHEE => EEEEE  otherwise   --H-- => CCCCC
    if((prdA2[i-2]=='2') && (prdA2[i-1]=='2') && (prdA2[i]=='1') && (prdA2[i+1]=='2') && (prdA2[i+2]=='2')){
      prdA2[i]='2';
    }
    else if((prdA2[i-2]!='1') && (prdA2[i-1]!='1') && (prdA2[i]=='1') && (prdA2[i+1]!='1') && (prdA2[i+2]!='1')){
      prdA2[i]='3';
    }
  }
  for(i=1;i<=Len;i++){  // EEHHEE => EEEEEE       otherwise        --HH-- => CCCCCC
    if((prdA2[i-2]=='2') && (prdA2[i-1]=='2') && (prdA2[i]=='1') && (prdA2[i+1]=='1') && (prdA2[i+2]=='2') && (prdA2[i+3]=='2')){
      prdA2[i]=prdA2[i+1]='2';
    }
    else if((prdA2[i-2]!='1') && (prdA2[i-1]!='1') && (prdA2[i]=='1') && (prdA2[i+1]=='1') && (prdA2[i+2]!='1') && (prdA2[i+3]!='1')){
      prdA2[i]=prdA2[i+1]='3';
    }
  } 

  //sprintf(filename,"%s/protein.ss2",path);
  sprintf(filename,"./protein.ss2");
  fp=fopen(filename,"w");
  fprintf(fp,"# SS prediction from single sequence developed by Pralay Mitra\n\n");
  for(i=1;i<=Len;i++) {
    fprintf(fp,"%4d %c ",i,fasta[i-1]);  
    if(prdA2[i]=='1')   fprintf(fp,"H ");  
    else if(prdA2[i]=='2')  fprintf(fp,"E ");
    else       fprintf(fp,"C ");
    fprintf(fp,"%7.3f%7.3f%7.3f\n",pout[i][0]/2.0,pout[i][1]/2.0,pout[i][2]/2.0);  
  }
  fclose(fp);

  for(i=0;i<Len;i++)   {
    prdA2[i]=prdA2[i+1];
  }
}

int SSPrediction::getSeq2SS::Read_score(){
  int i,j,k;
  FILE *fp;
  char line[5000]="XXX",aa='X';
  char amino[22]="XXX";
  float propensity[21][4]={0};
  float tmp1=0,tmp2=0,tmp3=0;

  if((fp=fopen(propfile,"r"))==NULL) {
    printf("cannot open propensity.txt\n");
    return(-1);
  }

  fgets(line,5000,fp);
  i=0;
  while(!feof(fp)) {
    fgets(line,5000,fp);
    if(feof(fp)) break;
    sscanf(line,"%c %f %f %f",&aa,&tmp1,&tmp2,&tmp3);
    amino[i]=aa;
    propensity[i][1]=tmp1;
    propensity[i][2]=tmp2;
    propensity[i][3]=tmp3;
    i++;
  }
  fclose(fp);

  for(i=1;i<=Len;i++) {
    for(j=0;j<20;j++) if(StrA1[i]==amino[j]) break;
    for(k=0;k<3;k++) {
      Score[i][k]=propensity[j][k+1];
    }
  }

  return(0);
}


int SSPrediction::getSeq2SS::Read_wgt(){
  int i,j;
  FILE *fpout=NULL;
  
  fpout=fopen(wgtfile,"r");       //weight file
  for(i=0;i<NLAYER1;i++) {       //330
    for(j=0;j<NLAYER2;j++) {     //150
      fscanf(fpout,"%f\n",&w1[i][j]);
      Dw1[i][j]=0.0;
    }
  }
  for(i=0;i<NLAYER2;i++) {       //150
    for(j=0;j<NLAYER3;j++){     //3
      fscanf(fpout,"%f\n",&w2[i][j]);
      Dw2[i][j]=0.0;
    }
  }
  for(i=0;i<MLAYER1;i++) {       //75
    for(j=0;j<MLAYER2;j++) {     //40
      fscanf(fpout,"%f\n",&v1[i][j]);
      Dv1[i][j]=0.0;
    }
  }
  for(i=0;i<MLAYER2;i++) {       //40
    for(j=0;j<MLAYER3;j++) {     //3
      fscanf(fpout,"%f\n",&v2[i][j]);
      Dv2[i][j]=0.0;
    }
  }
  fclose(fpout);

  return 0;
}

/****************************************************************/

int SSPrediction::getSeq2SS::Read_score1(char fasta[MAX_LEN],float divisor){
  int i,j,aa_num=0,Len;
  char *ncbicodes = "XAXCDEFGHIKLMNPQRSTVWXYXXX";
  float tmp3=0;

  utility::Blosum62 bsm;
  
  // read Score1[][]---------------->
  Len=strlen(fasta); 
  for(i=1;i<=Len;i++) {
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
      tmp3=Score1[i][j]/divisor;
      Score1[i][j]=tmp3/sqrt(1.0+tmp3*tmp3);
    }
    Score1[i][20]=0;
  }
  
  return Len;
}

/********************************************************************/
void SSPrediction::getSeq2SS::FeedNet(int I){
  int i,j,k;
  
  k=0;
  for(i=(I-SEG/2);i<(I-SEG/2)+SEG;i++){
    //---Score1-------->
    if(i<1 || i>Len){
      for(j=0;j<20;j++) {
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
  }
  
  for(i=k;i<=NLAYER1;i++)
    in1[i]=0;
  
  for(i=0;i<NLAYER3;i++)
    tgtout[i]=0.0;

  tgtout[(StrA2[I]-48)]=1.0; //StrA2[I]=48,49,50, native SS
}

/**************************************************************************/
void SSPrediction::getSeq2SS::Netout(){
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
void SSPrediction::getSeq2SS::FeedNet1(int I){
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
    else {
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
void SSPrediction::getSeq2SS::Netout1(){
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


