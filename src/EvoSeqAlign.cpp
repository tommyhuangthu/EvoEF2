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
/*                                SeqAlign.cpp                                */
/*                                                                            */
/******************************************************************************/
#pragma warning(disable:4996)
#pragma warning(disable:4244)
#include "EvoSeqAlign.h"
#include "EvoAminoName.h"
#include <string.h>


namespace Text {
  long int nnewline(FILE* fp);
}


// CONSTRUCTOR: CharSeq::SeqAlign::SeqAlign()
CharSeq::SeqAlign::SeqAlign(SequenceData dsInfo,char prffile[500],char phipsifile[500],float wss,float wsa,float wang) 
{
	int i,j;

	len=dsInfo.len;
	weight_ss=wss;
	weight_sa=wsa;
	weight_angle=wang;

	strcpy(seq,dsInfo.seq);
	strcpy(secstr1,dsInfo.ss1);
	strcpy(secstr2,dsInfo.ss2);
	strcpy(solacc1,dsInfo.sa1);
	strcpy(solacc2,dsInfo.sa2);

	readMkPrf(prffile);			// read the profile
	readPhiPsi(phipsifile);	// read target and designed phi-psi

	  // initialize alignment matrices
	D1 = getMtx();
	for(j=0; j<len; ++j) {
		MtxItem tmp(2*GAP_OPEN_PENALTY + j*GAP_XTN_PENALTY, false, false, false);
		D1[0][j] = tmp;
	}

	D2 = getMtx();
	for(i=0; i<len; ++i) {
		MtxItem tmp(2*GAP_OPEN_PENALTY + i*GAP_XTN_PENALTY, false, false, false);
		D2[i][0] = tmp;
	}

	D3 = getMtx();
	MtxItem tmp(Score(0,0), false, false, false);
	D3[0][0] = tmp;
	for(j=1; j<len; ++j) {
		MtxItem tmp(GAP_OPEN_PENALTY + (j-1)*GAP_XTN_PENALTY + Score(0,j), false, false, false);
		D3[0][j] = tmp;
	}
	for(i=1; i<len; ++i) {
		MtxItem tmp(GAP_OPEN_PENALTY + (i-1)*GAP_XTN_PENALTY + Score(i,0), false, false, false);
		D3[i][0] = tmp;
	}

	for(i=1; i<len; ++i) {
		updateD1(i, 0);
	}

	for(j=1; j<len; ++j) {
		updateD2(0, j);
	}

	// compute optimal alignment
	for(i=1; i<len; ++i) {
		for(j=1; j<len; ++j) {
			updateD1(i,j);
			updateD2(i,j);
			updateD3(i,j);
		}
	}
}


void CharSeq::SeqAlign::readPhiPsi(char phipsifile[500])
{
	FILE *fp;
	int i;
	float ang,ang1;

	phiT=new float[len];
	psiT=new float[len];
	phiQ=new float[len];
	psiQ=new float[len];

	fp=fopen(phipsifile,"r");
	if(fp==NULL) {
		for(i=0;i<len;i++) phiT[i]=psiT[i]=0.0;
	}
  else{
		i=0;
		while(fscanf(fp,"%f %f",&ang1,&ang)!=EOF){
      //fscanf(fp,"%*d %f %f",&ang1,&ang);
			//fscanf(fp,"%f %f",&ang1,&ang);
			phiT[i]=ang1;
			psiT[i]=ang;
			i++;
		}
    if(i!=len) printf("error in reading phi-psi file %d of phi-psi length but %d of sequence length\n",i,len);
    fclose(fp);
	}

	fp=fopen("./qphi.txt","r");
	if(fp==NULL) {
		for(i=0;i<len;i++) phiQ[i]=0.0;
	}
  else{
		i=0;
		while(fscanf(fp,"%*d %f",&ang)!=EOF){
			phiQ[i++]=ang;
		}
    if(i!=len) printf("error in reading qphi file %d of qphi length but %d of sequence length\n",i,len);
    fclose(fp);
	}

	fp=fopen("./qpsi.txt","r");
	if(fp==NULL){
		for(i=0;i<len;i++) psiQ[i]=0.0;
	}
  else {
		i=0;
		while(fscanf(fp,"%*d %f",&ang)!=EOF){
			psiQ[i++]=ang;
		}
    if(i!=len) printf("error in reading qpsi file %d of qpsi length but %d of sequence length\n",i,len);
    fclose(fp);
	}
}


// FUNCTION: Protein::Profile::readMkPrf()
// PURPOSE: reads a profile from the output file of the mkprf program.
// ARGUMENT:	---> pointer to the start of the file containing the profile.
// NOTES:	- the information about gaps is ignored.
void CharSeq::SeqAlign::readMkPrf(char* prffile) {
	FILE *fp=NULL;
	fp=fopen(prffile,"r");
	// count the number of residues
	
//	len = Text::nnewline(fp) - 2;
//	fseek(fp, 0, SEEK_SET);

	// allocate the profile
	prf = new float*[len];
	for(int i=0; i<len; ++i)
	    prf[i] = new float[AMINOS];

	// map amino acid indexes in input profile to AminoName enumerators
	AminoName map[AMINOS];
	char ami[2];
	for(int i=0; i<AMINOS; ++i) {
		fscanf(fp, "%s", ami);
		map[i] = charToAmino(ami[0]);
	}
	fscanf(fp, "%*s");

	// read the profile
	for(int i=0; i<len; ++i) {
		fscanf(fp, "%*s");
		for(int j=0; j<AMINOS; ++j) {
			fscanf(fp, "%f", &prf[i][map[j]]);
    }
		fscanf(fp, "%*s");
  }
	fclose(fp);
}

// CONSTRUCTOR: CharSeq::SeqAlign::MtxItem::MtxItem()
CharSeq::SeqAlign::MtxItem::MtxItem(const float s, const bool pv1, const bool pv2, const bool pv3){
	score = s;
	pvmtx[D1E] = pv1;
	pvmtx[D2E] = pv2;
	pvmtx[D3E] = pv3;
}

// CONSTRUCTOR: CharSeq::SeqAlign::MtxItem::MtxItem()
// ARGUMENTS:	--> max: maximum value among d1, d2, and d3 ---> d1, d2, d3: scores from D1, D2, and D3, respectively.
// PURPOSE:	assigns max to score and sets to true the flags of those matrices whose score equals max.
CharSeq::SeqAlign::MtxItem::MtxItem(const float max, const float d1, const float d2, const float d3){
	score=max;
	pvmtx[D1E] = (max == d1) ? true : false;
	pvmtx[D2E] = (max == d2) ? true : false;
	pvmtx[D3E] = (max == d3) ? true : false;
}


//	i1 = Query (Template)
//	i2 = Target (Decoy)
float CharSeq::SeqAlign::Score(const int i1, const int i2){
	float scosp,scoss,scosa,scoangle,phi,psi;
	AminoName aa = charToAmino(seq[i1]);	
	
	scosp = prf[i2][aa];
  //scoss = (secstr1[i1] == secstr2[i2]) ? 1.0 : 0.0;
	if(secstr1[i1] == secstr2[i2]) {
		scoss=1.0;
	}
  else if(secstr1[i1] == '3' || secstr2[i2] =='3') {
		scoss=0;
	}
  else {
		scoss=-1.0;
	}	
	if(solacc1[i1] == solacc2[i2]) {
		scosa=1.0;
	}
  else if(solacc1[i1] == '2' || solacc2[i2] == '2') {
		scosa=0.0;
	}
  else {
		scosa=-1.0;
	}
	
  psi=fabs(psiT[i2]-psiQ[i1]);
  if(psi >= 180.0 && psi <= 360.0) psi = 360.0 - fabs(psiT[i2]-psiQ[i1]);
  phi=fabs(phiT[i2]-phiQ[i1]);
  if(phi >= 180.0 && phi <= 360.0) phi = 360.0 - fabs(phiT[i2]-phiQ[i1]);

	if(i1==0 || i2==0)
    scoangle=psi; // N-Terminal
	else if(i1==len-1 || i2==len-1)
    scoangle=phi; // C-Terminal
	else
    scoangle=phi+psi;
	
	return scosp + weight_ss*scoss + weight_sa*scosa - weight_angle*scoangle;
}


// FUNCTION: CharSeq::SeqAlign::getMtx()		// PURPOSE: returns an empty alignment matrix
CharSeq::SeqAlign::MtxItem** CharSeq::SeqAlign::getMtx() const {
  MtxItem** mtx = new MtxItem*[len];
  for(int i=0; i<len; ++i)
    mtx[i] = new MtxItem[len];

  return mtx;
}


// FUNCTION: CharSeq::SeqAlign::updateD1()		// PURPOSE: fills in an item in D1
// ARGUMENTS:	 ---> i: item's row	 ---> j: item's column
void CharSeq::SeqAlign::updateD1(const int i, const int j) {
  const float d1 = D1[i-1][j].score + GAP_XTN_PENALTY;
  const float d2 = D2[i-1][j].score + GAP_OPEN_PENALTY;
  const float d3 = D3[i-1][j].score + GAP_OPEN_PENALTY;

	float max;

  if(d1>d2)
    max=d1;
  else
    max=d2;

  if(d3>max)
    max=d3;

  D1[i][j] = MtxItem(max, d1, d2, d3);
}

// FUNCTION: CharSeq::SeqAlign::updateD2()		// PURPOSE: fills in an item in D2
// ARGUMENTS:	 ---> i: item's row	---> j: item's column
void CharSeq::SeqAlign::updateD2(const int i, const int j) {
  const float d1 = D1[i][j-1].score + GAP_OPEN_PENALTY;
  const float d2 = D2[i][j-1].score + GAP_XTN_PENALTY;
  const float d3 = D3[i][j-1].score + GAP_OPEN_PENALTY;

	float max;

  if(d1>d2)
    max=d1;
  else
    max=d2;

  if(d3>max)
    max=d3;

  D2[i][j] = MtxItem(max, d1, d2, d3);
}

// FUNCTION: CharSeq::SeqAlign::updateD3()		// PURPOSE: fills in an item in D3
// ARGUMENTS:	 ---> i: item's row			---> j: item's column
void CharSeq::SeqAlign::updateD3(const int i, const int j) {
  
  const float d1 = D1[i-1][j-1].score + Score(i,j);
  const float d2 = D2[i-1][j-1].score + Score(i,j);
  const float d3 = D3[i-1][j-1].score + Score(i,j);

  float max;

	if(d1>d2)
    max=d1;
	else
    max=d2;

  if(d3>max)
    max=d3;

  D3[i][j] = MtxItem(max, d1, d2, d3);
}

// FUNCTION: CharSeq::SeqAlign::printMtxInfo() 		// PURPOSE: prints information from the D matrices, such as the scores and the	"previous" fields.
float CharSeq::SeqAlign::printMtxInfo() const {

	int a,b;
	float max,d1,d2,d3;

	a=len-1;
	b=len-1;
	d1=D1[a][b].score;
	d2=D2[a][b].score;
	d3=D3[a][b].score;

	if(d1>d2)
    max = d1;
	else
    max = d2;

  if(d3>max)
    max = d3;

	return max;
//	printf("%f\n", max);
}

// DESTRUCTOR: CharSeq::SeqAlign::~SeqAlign() 
CharSeq::SeqAlign::~SeqAlign(){
	for(int i=0; i<len; ++i)	delete[] D3[i];
	delete[] D3; 
	for(int i=0; i<len; ++i)	delete[] D2[i];
	delete[] D2; 
	for(int i=0; i<len; ++i)	delete[] D1[i];
	delete[] D1; 
	for(int i=0; i<len; ++i)
	delete[] prf[i];
	delete[] prf;
	
	delete[] phiT; 
	delete[] psiT; 
	delete[] phiQ; 
	delete[] psiQ; 
}

