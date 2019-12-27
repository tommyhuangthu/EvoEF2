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
/*                                xtd-seqali                                  */
/*                                                                            */
/* Computes the optimum alignment between two amino acid sequences (S1 & S2). */
/* The score of each pair (i,j) of positions takes into account not only      */
/* the similarity between the amino acids but also other features that        */
/* characterize those positions.                                              */
/*                                                                            */
/* In the current implementation, the score of each pair (i,j) of positions   */
/* is given by the secondary structures at positions i and j in S1 and S2,    */
/* respectively, and by the value of S1[i] at position j in the PSSM of S2.   */
/*                                                                            */
/* ARGUMENTS:                                                                 */
/* ---> S1:  path to the file containing amino acid sequence S1.              */
/*           The sequence is all on one line, starts at the start of the      */
/*           file, and is immediately followed by a line terminator.          */
/* ---> PF2: path to the file containing the PSSM of S2.                      */
/*           The PSSM can be output either by program makemat from the        */
/*           Impala suite, or by program mkprf (see notes).                   */
/* ---> SS1: path to the file containing the secondary structure of S1.       */
/*           The secondary structure is a horizontal sequence, with one       */
/*           character per residue position. It starts at the start of the    */
/*           file and is immediately followed by a line terminator.           */
/*           The i-th character represents the secondary structure at the     */
/*           i-th residue position (i=0,...,N-1, where N is the number of     */
/*           residues).            .                                          */
/* ---> SS2: path to the file containing the secondary structure of S2.       */
/*           The file format is the same as for SS1.                          */
/* ---> WSS: weight of the secondary structure in the alignment.              */
/*                                                                            */
/* NOTES:                                                                     */
/* - the PSSM output by program mkprf consists of a header line, N data lines,*/
/*   and a footer line.                                                       */
/*     The header line contains 21 columns, with one amino acid type per      */
/*   column, expressed in one-letter code. The last column, labeled '-'       */
/*   refers to the gap. The columns are separated by a certain number of      */
/*   blanks and the last column is immediately followd by a line terminator.  */
/*     The i-th data line contains the profile data for the i-th residue      */
/*   position, for i=0,...,N-1, where N is the number of residue positions.   */
/*   The j-th column is the score of the j-th amino acid type at that         */
/*   position, for j=0,...,19. The 20-th column is the score of the gap at    */
/*   that position. The 20-th column is immediately followed by a line        */
/*   terminator.                                                              */
/*     The footer line is identical to the header line. After its line        */
/*   terminator the file terminates.                                          */
/*                                                                            */
/******************************************************************************/
#pragma warning(disable:4996)
#pragma warning(disable:4305)
#pragma warning(disable:4244)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "EvoSeqAlign.h"
#include "EvoGetSeq2SA.h"
#include "EvoGetSeq2SS.h"
#include "EvoGetPhiPsi.h"
#include "EvoUtility.h"
#include "ErrorTracker.h"
#include "EvoAminoName.h"

float wss = 1.58;
float wsa = 2.45;
float wang= 1.00;

using namespace CharSeq;


namespace Text {
  int ncharLine(FILE* fp);
}

float EvolutionScore(char* programpath, char* prffile, char* ssfile, char* safile, char* phipsifile, char* fasseqfile){
  char evolutiondir[500]="";
  sprintf(evolutiondir,"%s/evolution",programpath);
  printf("evolution parameter file path is: %s\n", evolutiondir);

  FILE* fp = fopen(fasseqfile, "r");
  if(fp==NULL){
    printf("In %s %s %d, cannot open file for reading\n", __FILE__, __FUNCTION__, __LINE__);
    return IOError;
  }
  SequenceData dsInfo;
  dsInfo.len = Text::ncharLine(fp);
  dsInfo.seq = new char[dsInfo.len+1];
  fseek(fp, 0, SEEK_SET);
  fscanf(fp, "%s", dsInfo.seq);
  fclose(fp);


  dsInfo.ss1 = new char[dsInfo.len+1];
  SSPrediction::getSeq2SS SS(dsInfo.ss1,dsInfo.seq, evolutiondir);
  dsInfo.sa1 = new char[dsInfo.len+1];
  SAPrediction::getSeq2SA SA(dsInfo.sa1,dsInfo.seq,dsInfo.ss1, evolutiondir);

  //phi-psi prediction
  PhiPsiPrediction::getPhiPsi(dsInfo, evolutiondir);

  dsInfo.ss2 = new char[dsInfo.len+1];
  fp = fopen(ssfile, "r");
  fscanf(fp, "%s", dsInfo.ss2);
  fclose(fp);
  dsInfo.sa2 = new char[dsInfo.len+1];
  fp = fopen(safile, "r");
  fscanf(fp, "%s", dsInfo.sa2);
  fclose(fp);


  wsa /= (dsInfo.len*1.0);    // weight for solvent accessibility
  wss /= (dsInfo.len*1.0);    // weight for secondary structure
  wang/= (dsInfo.len*180.0);  // weight for phi-psi prediction

  CharSeq::SeqAlign sa(dsInfo, prffile, phipsifile, wss, wsa, wang);

  float max = sa.printMtxInfo();
  printf("%f\n",max);

  delete[] dsInfo.ss1;  delete[] dsInfo.ss2;
  delete[] dsInfo.sa1;  delete[] dsInfo.sa2;
  delete[] dsInfo.seq;

  return max;
  //return 1.0;
}


float EvolutionScore2(char* programpath, char* prffile, char* ssfile, char* safile, char* phipsifile, char* fasseq){
  char evolutiondir[500]="";
  sprintf(evolutiondir,"%s/evolution",programpath);
  //printf("evolution parameter file path is: %s\n", evolutiondir);

  SequenceData dsInfo;
  char seq[MAX_LEN];
  strcpy(seq,fasseq);
  dsInfo.len=strlen(seq);
  dsInfo.seq = new char[dsInfo.len+1];
  strcpy(dsInfo.seq,seq);

  dsInfo.ss1 = new char[dsInfo.len+1];
  SSPrediction::getSeq2SS SS(dsInfo.ss1,dsInfo.seq, evolutiondir);
  dsInfo.sa1 = new char[dsInfo.len+1];
  SAPrediction::getSeq2SA SA(dsInfo.sa1,dsInfo.seq,dsInfo.ss1, evolutiondir);

  //phi-psi prediction
  PhiPsiPrediction::getPhiPsi(dsInfo, evolutiondir);

  dsInfo.ss2 = new char[dsInfo.len+1];
  FILE* fp = fopen(ssfile, "r");
  fscanf(fp, "%s", dsInfo.ss2);
  fclose(fp);
  dsInfo.sa2 = new char[dsInfo.len+1];
  fp = fopen(safile, "r");
  fscanf(fp, "%s", dsInfo.sa2);
  fclose(fp);


  wsa /= (dsInfo.len*1.0);    // weight for solvent accessibility
  wss /= (dsInfo.len*1.0);    // weight for secondary structure
  wang/= (dsInfo.len*180.0);  // weight for phi-psi prediction

  CharSeq::SeqAlign sa(dsInfo, prffile, phipsifile, wss, wsa, wang);

  float max = sa.printMtxInfo();
  //printf("%f\n",max);

  delete[] dsInfo.ss1;  delete[] dsInfo.ss2;
  delete[] dsInfo.sa1;  delete[] dsInfo.sa2;
  delete[] dsInfo.seq;

  return -1.0*max;
}


float EvolutionScore3(char* programpath, char* prffile, char* fasseq){
  char evolutiondir[500]="";
  sprintf(evolutiondir,"%s/evolution",programpath);
  //printf("evolution parameter file path is: %s\n", evolutiondir);

  char seq[MAX_LEN];
  strcpy(seq,fasseq);
  int len=strlen(seq);

  //read profile
  float** prf;
  FILE *fp=NULL;
  fp=fopen(prffile,"r");
  prf = new float*[len];
  for(int i=0; i<len; ++i)
    prf[i] = new float[AMINOS];
  AminoName map[AMINOS];
  char ami[2];
  for(int i=0; i<AMINOS; ++i) {
    fscanf(fp, "%s", ami);
    map[i] = charToAmino(ami[0]);
  }
  fscanf(fp, "%*s");
  for(int i=0; i<len; ++i) {
    fscanf(fp, "%*s");
    for(int j=0; j<AMINOS; ++j) {
      fscanf(fp, "%f", &prf[i][map[j]]);
    }
    fscanf(fp, "%*s");
  }
  fclose(fp);

  //calculate profile score without alignment
  float score=0;
  for(int i=0; i<len; i++){
    AminoName map=charToAmino(seq[i]);
    score += prf[i][map];
  }

  for(int i=0; i<len; i++) delete [] prf[i];
  delete [] prf;

  return -1.0*score;
}
