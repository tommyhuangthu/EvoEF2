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

#ifndef CHAIN_H
#define CHAIN_H

#include "Residue.h"

typedef enum _Type_Chain{
  Type_Chain_Protein, 
  Type_Chain_SmallMol, 
  Type_Chain_Nucleotide, 
  Type_Chain_MetalIon, 
  Type_Chain_Water
}Type_Chain;

int ChainTypeConvertFromString(char* typeName, Type_Chain* type);
Type_Chain ChainTypeIdentifiedFromResidueName(char *resiName);

typedef struct _Chain{
  Residue* residues;
  char name[MAX_LENGTH_CHAIN_NAME+1];
  int residueNum;
  Type_Chain type;
} Chain;

int ChainCreate(Chain* pThis);
int ChainDestroy(Chain* pThis);
int ChainCopy(Chain* pThis, Chain* pOther);
char* ChainGetName(Chain* pThis);
int ChainSetName(Chain* pThis, char* newName);
Type_Chain ChainGetType(Chain* pThis);
int ChainSetType(Chain* pThis, Type_Chain newType);
int ChainGetResidueCount(Chain* pThis);
Residue* ChainGetResidue(Chain* pThis, int index);
int ChainInsertResidue(Chain* pThis, int index, Residue* pNewResi);
int ChainRemoveResidue(Chain* pThis, int index);
int ChainAppendResidue(Chain* pThis, Residue* pNewResi);
int ChainReadCoordinate(Chain* pThis, char* coordinateFile);
int ChainCalcAllAtomXYZ(Chain* pThis, ResiTopoSet* topos);
int ChainShowInPDBFormat(Chain* pThis, int resiIndex, int atomIndex, BOOL showHydrogen, FILE* pFile);
int ChainFindResidueByPosInChain(Chain* pThis, int posInchain, int *index);
int ChainShowAtomParameter(Chain* pThis);
int ChainShowBondInformation(Chain* pThis);
#endif