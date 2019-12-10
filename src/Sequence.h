///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "Structure.h"
#include "EnergyMatrix.h"
#include "EnergyFunction.h"

#define  MAX_SEQ_LEN  2000

typedef struct _Sequence{
	double        etot;
  double        ephy;
  double        eevo;
  double        ebpf;
	int           designSiteCount;
	IntArray      rotamerIndices;
} Sequence;

int SequenceCreate(Sequence* pThis);
int SequenceDestroy(Sequence* pThis);
int SequenceCopy(Sequence* pThis, Sequence* pOther);
int SequenceWriteDesignRotamer(Sequence* pThis,Structure*pStructure,int index,FILE* pFile);
int SequenceWriteDesignFasta(Sequence* pThis,Structure*pStructure,char*deschns,int index,FILE* pFile2);


int StructureShowMinEnergyStructure(Structure* pStructure, Sequence* bestSeq, char* pdbfile);
int StructureGetWholeSequence(Structure* pStructure,Sequence* sequence,char*desChnName,char* fastaseq);
int StructureShowDesignedResiduesOnly(Structure* pStructure, Sequence* sequence, char* pdbfile);
int SequenceWriteDesignFastaForSCP(Sequence* pThis,Structure*pStructure,char*deschns,int index,FILE* pFile2);
int SequenceWriteSCPtorsionAndRmsd(Sequence* pThis,Structure*pStructure,char*deschns,FILE* pTorsion,FILE* pRmsd);

#endif // SEQUENCE_H
