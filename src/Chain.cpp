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

#include "Chain.h"
#include <string.h>

int ChainTypeConvertFromString(char* typeName, Type_Chain* type){
  if(strcmp(typeName, "PROTEIN")==0){
    *type = Type_Chain_Protein;
  }
  else if(strcmp(typeName, "SMALLMOL")==0){
    *type = Type_Chain_SmallMol;
  }
  else if(strcmp(typeName, "METALION")==0){
    *type = Type_Chain_MetalIon;
  }
  else if(strcmp(typeName, "NUCLEOTIDE")==0){
    *type = Type_Chain_Nucleotide;
  }
  else if(strcmp(typeName, "WATER")==0){
    *type = Type_Chain_Water;
  }
  else{
    return ValueError;
  }
  return Success;
}

Type_Chain ChainTypeIdentifiedFromResidueName(char *resiName){
  if(strcmp(resiName, "ALA") == 0
    || strcmp(resiName, "ARG") == 0
    || strcmp(resiName, "ASN") == 0
    || strcmp(resiName, "ASP") == 0
    || strcmp(resiName, "CYS") == 0
    || strcmp(resiName, "GLN") == 0
    || strcmp(resiName, "GLU") == 0
    || strcmp(resiName, "GLY") == 0
    || strcmp(resiName, "HIS") == 0
    || strcmp(resiName, "HSD") == 0
    || strcmp(resiName, "HSE") == 0
    || strcmp(resiName, "HSP") == 0
    || strcmp(resiName, "ILE") == 0
    || strcmp(resiName, "LEU") == 0
    || strcmp(resiName, "LYS") == 0
    || strcmp(resiName, "MET") == 0
    || strcmp(resiName, "PHE") == 0
    || strcmp(resiName, "PRO") == 0
    || strcmp(resiName, "SER") == 0
    || strcmp(resiName, "THR") == 0
    || strcmp(resiName, "TRP") == 0
    || strcmp(resiName, "TYR") == 0
    || strcmp(resiName, "VAL") == 0){
    return Type_Chain_Protein;
  }
  else if(strcmp(resiName, "GUA") == 0
    || strcmp(resiName, "ADE") == 0
    || strcmp(resiName, "CYT") == 0
    || strcmp(resiName, "THY") == 0
    || strcmp(resiName, "URA") == 0){
    return Type_Chain_Nucleotide;
  }
  else if(strcmp(resiName, "HOH") == 0
    || strcmp(resiName, "H2O") == 0
    || strcmp(resiName, "WAT") == 0){
    return Type_Chain_Water;
  }
  else{
    return Type_Chain_SmallMol;
  }
}

int ChainCreate(Chain* pThis){
  ChainSetName(pThis, "");
  pThis->residueNum = 0;
  pThis->residues = NULL;
  return Success;
}

int ChainDestroy(Chain* pThis){
  for(int i=0;i<pThis->residueNum;i++){
    ResidueDestroy(&pThis->residues[i]);
  }
  free(pThis->residues);
  pThis->residues = NULL;
  pThis->residueNum = 0;
  return Success;
}

int ChainCopy(Chain* pThis, Chain* pOther){
  ChainDestroy(pThis);
  ChainSetName(pThis, pOther->name);
  ChainSetType(pThis, pOther->type);

  pThis->residueNum = pOther->residueNum;
  pThis->residues = (Residue*)malloc(sizeof(Residue)*pThis->residueNum);
  for(int i=0;i<pThis->residueNum;i++){
    int result;
    ResidueCreate(&pThis->residues[i]);
    result = ResidueCopy(&pThis->residues[i], &pOther->residues[i]);
    if(FAILED(result))
      return result;
  }
  return Success;
}

char*  ChainGetName(Chain* pThis){
  return pThis->name;
}

int ChainSetName(Chain* pThis, char* newName){
  if(strlen(newName)>MAX_LENGTH_CHAIN_NAME){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s function %s() line %d, name is too long", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->name, newName);
  return Success;
}

Type_Chain ChainGetType(Chain* pThis){
  return pThis->type;
}

int ChainSetType(Chain* pThis, Type_Chain newType){
  pThis->type = newType;
  return Success;
}

int ChainGetResidueCount(Chain* pThis){
  return pThis->residueNum;
}

Residue* ChainGetResidue(Chain* pThis, int index){
  if(index<0 || index>=ChainGetResidueCount(pThis)){
    return NULL;
  }
  else{
    return pThis->residues + index;
  }

}

int ChainInsertResidue(Chain* pThis, int index, Residue* pNewResi){
  if(index<0 || index>pThis->residueNum)
    return IndexError;

  int newCount = pThis->residueNum+1;
  pThis->residues = (Residue*)realloc(pThis->residues, sizeof(Residue)*newCount);

  ResidueCreate(&pThis->residues[newCount-1]);
  for(int i=newCount-1;i>index;i--){
    ResidueCopy(&pThis->residues[i], &pThis->residues[i-1]);
    ResidueSetPosInChain(&pThis->residues[i], pThis->residues[i-1].posInChain);
  }

  ResidueCopy(&pThis->residues[index], pNewResi);
  ResidueSetChainName(&pThis->residues[index], pThis->name);
  ResidueSetPosInChain(&pThis->residues[index], pNewResi->posInChain);
  pThis->residueNum = newCount;

  return Success;
}

int ChainRemoveResidue(Chain* pThis, int index){
  if(index<0 || index>=pThis->residueNum)
    return IndexError;

  for(int i=index;i<pThis->residueNum-1;i++){
    ResidueCopy(&pThis->residues[i], &pThis->residues[i+1]);
    ResidueSetPosInChain(&pThis->residues[i], ResidueGetPosInChain(&pThis->residues[i+1]));
  }
  ResidueDestroy(&pThis->residues[pThis->residueNum-1]);
  (pThis->residueNum)--;
  return Success;
}

int ChainAppendResidue(Chain* pThis, Residue* pNewResi){
  return ChainInsertResidue(pThis, pThis->residueNum, pNewResi);
}


int ChainCalcAllAtomXYZ(Chain* pThis, ResiTopoSet* topos){
  for(int i=0;i<ChainGetResidueCount(pThis);i++){
    Residue* prevResi = ChainGetResidue(pThis, i-1);
    Residue* nextResi = ChainGetResidue(pThis, i+1);
    Residue* curResi = ChainGetResidue(pThis, i);

    int result = ResidueCalcAllAtomXYZ(curResi, topos, prevResi, nextResi);
    if(FAILED(result)){
      char errMsg[MAX_LENGTH_ERR_MSG+1];
      sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
      TraceError(errMsg, result);
      return result;
    }
  }
  return Success;
}

int ChainShowInPDBFormat(Chain* pThis, int resiIndex, int atomIndex, BOOL showHydrogen, FILE* pFile){
  char header[7];
  char chainName[2];
  if(pThis->type==Type_Chain_Protein){
    strcpy(header, "ATOM  ");
    strcpy(chainName, " ");
    chainName[0] = pThis->name[strlen(pThis->name)-1];
  }
  else{
    strcpy(header, "HETATM");
    strcpy(chainName, " ");
  }
  int atomCounter = atomIndex;
  for(int i=0;i<pThis->residueNum;i++){
    ResidueShowInPDBFormat(&pThis->residues[i], header, chainName, atomCounter, ResidueGetPosInChain(&pThis->residues[i])+resiIndex, showHydrogen, pFile);
    atomCounter += ResidueGetAtomCount(&pThis->residues[i]);
  }
  return Success;
}

int ChainFindResidueByPosInChain(Chain* pThis, int posInchain, int *index){
  for(int i=0; i<ChainGetResidueCount(pThis); ++i){
    Residue* pResi = ChainGetResidue(pThis, i);
    if(ResidueGetPosInChain(pResi)==posInchain){
      *index=i;
      return Success;
    }
  }
  return DataNotExistError;
}

int ChainShowAtomParameter(Chain* pThis){
  for(int i=0; i<ChainGetResidueCount(pThis); i++){
    ResidueShowAtomParameter(ChainGetResidue(pThis,i));
  }
  return Success;
}


int ChainShowBondInformation(Chain* pThis){
  printf("------------------Bonds information for Chain %s--------------\n", ChainGetName(pThis));
  for(int i=0; i<ChainGetResidueCount(pThis); i++){
    Residue* pResidue=ChainGetResidue(pThis,i);
    printf("Residue %s PosInChain %d:\n", ResidueGetName(pResidue),ResidueGetPosInChain(pResidue));
    ResidueShowBondInformation(pResidue);
  }
  return Success;
}