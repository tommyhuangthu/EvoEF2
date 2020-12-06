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

#include "Residue.h"
#include <string.h>
#include <ctype.h>


int ResidueCreate(Residue* pThis){
  strcpy(pThis->name, "");
  strcpy(pThis->chainName, "");
  pThis->posInChain = -1;
  pThis->designSiteType = Type_ResidueDesignType_Fixed;
  AtomArrayCreate(&pThis->atoms);
  StringArrayCreate(&pThis->patches);
  BondSetCreate(&pThis->bonds);
  pThis->resiTerm = Type_ResidueIsNotTerminal;
  pThis->nCbIn8A = 0;
  pThis->internalEnergy = 0;
  pThis->backboneEnergy = 0;
  pThis->phipsi[0]=-120;
  pThis->phipsi[1]=0;
  DoubleArrayCreate(&pThis->xtorsions,0);
  pThis->dunbrack=0;
  return Success;
}

int ResidueDestroy(Residue* pThis){
  BondSetDestroy(&pThis->bonds);
  AtomArrayDestroy(&pThis->atoms);
  StringArrayDestroy(&pThis->patches);
  DoubleArrayDestroy(&pThis->xtorsions);
  return Success;
}

int ResidueCopy(Residue* pThis, Residue* pOther){
  ResidueDestroy(pThis);
  strcpy(pThis->name, pOther->name);
  strcpy(pThis->chainName, pOther->chainName);
  pThis->posInChain = pOther->posInChain;
  pThis->designSiteType = pOther->designSiteType;
  pThis->nCbIn8A = pOther->nCbIn8A;
  AtomArrayCreate(&pThis->atoms);
  AtomArrayCopy(&pThis->atoms, &pOther->atoms);
  StringArrayCreate(&pThis->patches);
  StringArrayCopy(&pThis->patches, &pOther->patches);
  BondSetCopy(&pThis->bonds, &pOther->bonds);
  pThis->resiTerm = pOther->resiTerm;
  pThis->phipsi[0]=pOther->phipsi[0];
  pThis->phipsi[1]=pOther->phipsi[1];
  DoubleArrayCopy(&pThis->xtorsions,&pOther->xtorsions);
  pThis->dunbrack = pOther->dunbrack;
  return Success;
}

char* ResidueGetName(Residue* pThis){
  return pThis->name;
}

int ResidueSetName(Residue* pThis, char* newName){
  if(strlen(newName)>MAX_LENGTH_RESIDUE_NAME){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->name, newName);
  return Success;
}

char* ResidueGetChainName(Residue* pThis){
  return pThis->chainName;
}

int ResidueSetChainName(Residue* pThis, char* newChainName){
  int i;
  if(strlen(newChainName)>MAX_LENGTH_CHAIN_NAME)
  {
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->chainName, newChainName);
  for(i=0;i<ResidueGetAtomCount(pThis);i++){
    Atom* pCurAtom = ResidueGetAtom(pThis, i);
    AtomSetChainName(pCurAtom, newChainName);
  }
  return Success;
}

int ResidueGetPosInChain(Residue* pThis){
  return pThis->posInChain;
}

int ResidueSetPosInChain(Residue* pThis, int newPosInChain){
  int i;
  //if(newPosInChain < 0){
  //  char errMsg[MAX_LENGTH_ERR_MSG+1];
  //  int errorCode = ValueError;
  //  sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
  //  TraceError(errMsg, errorCode);
  //  return errorCode;
  //}
  pThis->posInChain = newPosInChain;
  for(i=0;i<ResidueGetAtomCount(pThis);i++){
    Atom* pCurAtom = ResidueGetAtom(pThis, i);
    AtomSetPosInChain(pCurAtom, newPosInChain);
  }

  return Success;
}

int ResidueSetDesignSiteFlag(Residue* pThis, Type_ResidueDesignType newFlag){
  pThis->designSiteType = newFlag;
  return Success;
}

double ResidueGetCharge(Residue* pThis){
  return AtomArrayCalcTotalCharge(&pThis->atoms);
}

int ResidueGetPolarity(Residue* pThis, Type_ResiduePolarity* pPolarity){
  // if the total charge is not zero, the type is Type_ResiduePolarity_Charged;
  // if the total charge is zero, but any non-root atom has non-zero charge, the type is Type_ResiduePolarity_Polar;
  // if the total charge is zero and all non-root atoms have zero charge, the type is Type_ResiduePolarity_NonPolar;
  BOOL nonRootAtomCharged = FALSE;
  for(int i=0;i<ResidueGetAtomCount(pThis);i++){
    Atom* pAtom;
    BOOL nonRoot;
    BOOL charged;
    pAtom = ResidueGetAtom(pThis, i);
    nonRoot = !(pAtom->isBBAtom);
    charged = fabs(pAtom->charge)>MIN_ZERO_TOLERANCE;
    if(nonRoot && charged)
      nonRootAtomCharged = TRUE;
  }
  BOOL residueCharged = fabs(ResidueGetCharge(pThis))>MIN_ZERO_TOLERANCE;
  if(residueCharged){
    *pPolarity = Type_ResiduePolarity_Charged;
  }
  else if(nonRootAtomCharged){
    *pPolarity = Type_ResiduePolarity_Polar;
  }
  else{
    *pPolarity = Type_ResiduePolarity_NonPolar;
  }
  return Success;
}

int ResidueGetAtomCount(Residue* pThis){
  return AtomArrayGetCount(&pThis->atoms);
}

Atom* ResidueGetAtom(Residue* pThis, int index){
  return AtomArrayGet(&pThis->atoms, index);
}

Atom* ResidueGetAtomByName(Residue* pThis, char* atomName){
  int index;
  int result = ResidueFindAtom(pThis, atomName, &index);
  if(FAILED(result)){
    return NULL;
  }
  else{
    return ResidueGetAtom(pThis, index);
  }
}

int ResidueFindAtom(Residue* pThis, char* atomName, int* pIndex){
  return AtomArrayFind(&pThis->atoms, atomName, pIndex);
}

int ResidueGetAtomXYZ(Residue* pThis, char* atomName, XYZ* pXYZ){
  Atom* pAtom = ResidueGetAtomByName(pThis, atomName);
  if( pAtom == NULL || pAtom->isXyzValid==FALSE){
    return DataNotExistError;
  }
  *pXYZ = pAtom->xyz;
  return Success;
}

AtomArray* ResidueGetAllAtoms(Residue* pThis){
  return &pThis->atoms;
}

int ResidueInsertAtom(Residue* pThis, int newIndex, Atom* pNewAtom){
  int index;
  if(newIndex < 0 || newIndex > AtomArrayGetCount(&pThis->atoms)){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = IndexError;
    sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  else if(FAILED(ResidueFindAtom(pThis, AtomGetName(pNewAtom), &index))){
    return AtomArrayInsert(&pThis->atoms, newIndex, pNewAtom);
  }
  else{
    return AtomCopy(ResidueGetAtom(pThis, index), pNewAtom);
  }
}

int ResidueAddAtom(Residue* pThis, Atom* pNewAtom){
  int index;
  if(FAILED(ResidueFindAtom(pThis, AtomGetName(pNewAtom), &index)))
  {
    return AtomArrayAppend(&pThis->atoms, pNewAtom);
  }else
  {
    return AtomCopy(ResidueGetAtom(pThis, index), pNewAtom);
  }
}

int ResidueAddAtomsFromAtomParams(Residue* pThis, AtomParamsSet* pAtomParams){
  int count;
  int result = AtomParamsSetGetAtomCount(pAtomParams, ResidueGetName(pThis), &count);
  if(FAILED(result)){
    return result;
  }

  for(int i=0;i<count;i++){
    BOOL atomExist = FALSE;
    Atom newAtom;
    AtomCreate(&newAtom);
    AtomParamsSetGetAtomParam(pAtomParams, ResidueGetName(pThis), i, &newAtom);
    for(int j = 0; j < ResidueGetAtomCount(pThis); j++){
      Atom *pAtom1 = ResidueGetAtom(pThis, j);
      if(strcmp(AtomGetName(pAtom1), AtomGetName(&newAtom)) == 0){
        Atom tempAtom;
        AtomCreate(&tempAtom);
        AtomCopy(&tempAtom, pAtom1);
        AtomCopy(pAtom1, &newAtom);
        AtomSetChainName(pAtom1,AtomGetChainName(&tempAtom));
        AtomSetPosInChain(pAtom1,AtomGetPosInChain(&tempAtom));
        pAtom1->xyz = tempAtom.xyz;
        pAtom1->isXyzValid = tempAtom.isXyzValid;
        atomExist = TRUE;
        AtomDestroy(&tempAtom);
        break;
      }
    }
    if(atomExist == FALSE){
      AtomSetChainName(&newAtom,ResidueGetChainName(pThis));
      AtomSetPosInChain(&newAtom, ResidueGetPosInChain(pThis));
      ResidueAddAtom(pThis, &newAtom);
    }
    AtomDestroy(&newAtom);
  }
  return Success;
}

int ResidueDeleteAtom(Residue* pThis, char* atomName){
  int index;
  if(FAILED(ResidueFindAtom(pThis, atomName, &index))){
    return DataNotExistError;
  }
  else{
    return AtomArrayRemove(&pThis->atoms, index);
  }
}


int ResidueReadXYZFromPDB(Residue* pThis, FileReader* pPDBFileReader){
  // the following three lines show the format of a pdb file
  // type, serial, name, altLoc, resName, chainID, resPos, ..., X,   Y,   Z
  // 0,    6,      12,   16,     17,      21,      22,     ..., 30,  38,  46
  // 6,    5,       4,    1,      4,       1,       5,   ...,    8,   8,  8
  BOOL firstLine = TRUE;
  char iniSeqPos[MAX_LENGTH_ONE_LINE_IN_FILE+1] = "UNKNOWN";
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(!FAILED(FileReaderGetNextLine(pPDBFileReader, line))){
    char strType[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    char strAtomName[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    char strResName[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    char strResPos[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    char strX[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    char strY[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    char strZ[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    ExtractTargetStringFromSourceString(strType, line, 0, 4);
    
    //for smallmol residue, read its vdw energies if possible
    //if(strcmp(strType, "ENER") == 0){
    //  StringArray strings;
    //  StringArrayCreate(&strings);
    //  StringArraySplitString(&strings, line, ' ');
    //  pThis->internalEnergy = atof(StringArrayGet(&strings, 2));
    //  pThis->backboneEnergy = atof(StringArrayGet(&strings, 4));
    //  StringArrayDestroy(&strings);
    //}


    if(strcmp(strType, "ATOM")!=0){
      continue;
    }
    ExtractTargetStringFromSourceString(strAtomName, line, 12, 4);
    ExtractTargetStringFromSourceString(strResName, line, 17, 4);
    ExtractTargetStringFromSourceString(strResPos, line, 22, 5);
    ExtractTargetStringFromSourceString(strX, line, 30, 8);
    ExtractTargetStringFromSourceString(strY, line, 38, 8);
    ExtractTargetStringFromSourceString(strZ, line, 46, 8);
    // if the residue is HIS, use HSD or HSE
    if(strcmp(strResName, "HIS")==0  &&  strcmp(pThis->name, "HSD")==0){
      strcpy(strResName, "HSD");
    }
    else if(strcmp(strResName, "HIS")==0  &&  strcmp(pThis->name, "HSE")==0){
      strcpy(strResName, "HSE");
    }
    // the CD1 atom of ILE in pdb is altered into CD
    if(strcmp(strResName, "ILE")==0 && strcmp(strAtomName, "CD1")==0){
      strcpy(strAtomName, "CD");
    }
    // read the hydrogen atoms, sometimes the hydrogen atom can have a name like '1HG1' in residue VAL
    if(strAtomName[0]=='H'){
      continue;
    }
    else if(isdigit(strAtomName[0]) && strAtomName[1]=='H' && (int)strlen(strAtomName)==4){
      char tempName[MAX_LENGTH_ATOM_NAME+1];
      tempName[0]=strResName[1]; tempName[1]=strResName[2]; tempName[2]=strResName[3]; tempName[3]=strResName[0]; tempName[4]='\0';
      strcpy(strAtomName,tempName);
      continue;
    }

    if(firstLine){
      strcpy(iniSeqPos, strResPos);
      firstLine = FALSE;
    }
    else{
      // read the line of a new residue
      if(strcmp(iniSeqPos, strResPos) != 0){
        FileReaderSetCurrentPos(pPDBFileReader, FileReaderGetCurrentPos(pPDBFileReader)-1);
        return Success;
      }
    }
    // read OXT coordinate from PDB file instead of recalculation
    if(strcmp(strAtomName, "OXT") == 0){
      Atom atomOXT;
      AtomCreate(&atomOXT);
      AtomCopy(&atomOXT, ResidueGetAtomByName(pThis, "O"));
      atomOXT.xyz.X = atof(strX); atomOXT.xyz.Y = atof(strY); atomOXT.xyz.Z = atof(strZ);
      atomOXT.isXyzValid = TRUE;
      atomOXT.charge = -0.55;
      strcpy(atomOXT.name, "OXT");
      ResidueGetAtomByName(pThis, "O")->charge = -0.55;
      ResidueGetAtomByName(pThis, "C")->charge = 0.1;
      AtomArrayAppend(&pThis->atoms, &atomOXT);
      AtomDestroy(&atomOXT);
    }
    else{
      Atom* pAtom = ResidueGetAtomByName(pThis, strAtomName);
      if(pAtom==NULL || pAtom->isXyzValid==TRUE) continue;
      //if(pAtom == NULL){
      //  int result = FormatError;
      //  char errMsg[MAX_LENGTH_ERR_MSG+1];
      //  sprintf(errMsg, "in file %s function %s() line %d, residue %s do not have atom %s when reading line\n%s", __FILE__, __FUNCTION__, __LINE__, strResName, strAtomName, line);
      //  TraceError(errMsg, result);
      //  continue;
      //}
      pAtom->xyz.X = atof(strX); pAtom->xyz.Y = atof(strY); pAtom->xyz.Z = atof(strZ);
      pAtom->isXyzValid = TRUE;
    }
  }
  return Success;
}

BondSet* ResidueGetBonds(Residue* pThis){
  return &pThis->bonds;
}
int ResidueAddBondsFromResiTopos(Residue* pThis, ResiTopoSet* pResiTopoCollection){
  ResidueTopology resiTopo;
  ResidueTopologyCreate(&resiTopo);
  if(FAILED(ResiTopoSetGet(pResiTopoCollection, pThis->name, &resiTopo))){
    return DataNotExistError;
  }
  BondSetCopy(&pThis->bonds, ResidueTopologyGetBonds(&resiTopo));
  ResidueTopologyDestroy(&resiTopo);
  return Success;
}

int ResidueShowInPDBFormat(Residue* pThis, char* header, char* chainName, int atomIndex, int resiIndex, BOOL showHydrogen, FILE* pFile){
  AtomArrayShowInPDBFormat(&pThis->atoms, header, ResidueGetName(pThis), chainName, atomIndex, resiIndex, showHydrogen, pFile);
  return Success;
}


int ResiduePatch(Residue* pThis, char* patchName, AtomParamsSet* pAtomParam, ResiTopoSet* pTopos){
  // delete old atoms
  ResidueTopology patchResiTopo;
  ResidueTopologyCreate(&patchResiTopo);
  int result = ResiTopoSetGet(pTopos, patchName, &patchResiTopo);
  if(FAILED(result)) return result;
  StringArray* deleteAtoms = ResidueTopologyGetDeletes(&patchResiTopo);
  for(int i=0;i<StringArrayGetCount(deleteAtoms);i++){
    ResidueDeleteAtom(pThis, StringArrayGet(deleteAtoms, i));
  }

  // Temporarily reset its name to that of the patch residue, to find topology in the topology collection
  char resiName[MAX_LENGTH_RESIDUE_NAME];
  strcpy(resiName, ResidueGetName(pThis));
  ResidueSetName(pThis, patchName);
  result = ResidueAddAtomsFromAtomParams(pThis, pAtomParam);
  if(FAILED(result)) return result;
  ResidueSetName(pThis, resiName);

  // new patches are stored at the head
  StringArrayInsert(&pThis->patches, 0, patchName);

  // deal with the bonds
  BondSet* pPatchBonds = ResidueTopologyGetBonds(&patchResiTopo);
  for(int i=0;i<BondSetGetCount(pPatchBonds);i++){
    Bond* pCurBond = BondSetGet(pPatchBonds, i);
    BondSetAdd(&pThis->bonds, BondGetFromName(pCurBond), BondGetToName(pCurBond), BondGetType(pCurBond));
  }
  BondSet newBonds;
  BondSetCreate(&newBonds);
  for(int i=0;i<BondSetGetCount(&pThis->bonds);i++){
    Bond* pCurBond = BondSetGet(&pThis->bonds, i);
    char* atomName = BondGetFromName(pCurBond);
    if( atomName[0]!='+' && atomName[0]!='-' && ResidueGetAtomByName(pThis, atomName)==NULL){
      continue;
    }
    atomName = BondGetToName(pCurBond);
    if( atomName[0]!='+' && atomName[0]!='-' && ResidueGetAtomByName(pThis, atomName)==NULL){
      continue;
    }
    BondSetAdd(&newBonds, BondGetFromName(pCurBond), BondGetToName(pCurBond), BondGetType(pCurBond));
  }
  BondSetCopy(&pThis->bonds, &newBonds);
  BondSetDestroy(&newBonds);
  ResidueTopologyDestroy(&patchResiTopo);
  return Success;
}

// there is a bug in previous ResiduePatch for CTER patching, because it will recalculate the coordinate of atom O and OXT
// but most time the coordinates of O and OXT will be determined, so it should be directly read from file instead of recalculating
int ResiduePatchCTER(Residue* pThis, char* patchName, AtomParamsSet* pAtomParam,ResiTopoSet* pTopos){
  // delete old atoms
  ResidueTopology patchResiTopo;
  ResidueTopologyCreate(&patchResiTopo);
  int result = ResiTopoSetGet(pTopos, "CTER", &patchResiTopo);
  if(FAILED(result)) return result;
  // do not delete atom O;
  //StringArray* deleteAtoms = ResidueTopologyGetDeletes(&patchResiTopo);
  //for(int i=0;i<StringArrayGetCount(deleteAtoms);i++){
  //  ResidueDeleteAtom(pThis, StringArrayGet(deleteAtoms, i));
  //}

  // Temporarily reset its name to that of the patch residue, to find topology in the topology collection
  char resiName[MAX_LENGTH_RESIDUE_NAME];
  strcpy(resiName, ResidueGetName(pThis));
  ResidueSetName(pThis, patchName);
  result = ResidueAddAtomsFromAtomParams(pThis, pAtomParam);
  // after patching, please check the chain name of the atom
  //for(int i=1; i<ResidueGetAtomCount(pThis); i++){
  //  Atom* pAtom = ResidueGetAtom(pThis,i);
  //  if(strcmp(AtomGetChainName(pAtom),"")==0){
  //    AtomSetChainName(pAtom,AtomGetChainName(ResidueGetAtom(pThis,0)));
  //  }
  //}
  if(FAILED(result)) return result;
  ResidueSetName(pThis, resiName);

  // new patches are stored at the head
  StringArrayInsert(&pThis->patches, 0, patchName);
  // deal with the bonds
  BondSet* pPatchBonds = ResidueTopologyGetBonds(&patchResiTopo);
  for(int i=0;i<BondSetGetCount(pPatchBonds);i++){
    Bond* pCurBond = BondSetGet(pPatchBonds, i);
    BondSetAdd(&pThis->bonds, BondGetFromName(pCurBond), BondGetToName(pCurBond), BondGetType(pCurBond));
  }
  BondSet newBonds;
  BondSetCreate(&newBonds);
  for(int i=0;i<BondSetGetCount(&pThis->bonds);i++){
    Bond* pCurBond = BondSetGet(&pThis->bonds, i);
    char* atomName = BondGetFromName(pCurBond);
    if( atomName[0]!='+' && atomName[0]!='-' && ResidueGetAtomByName(pThis, atomName)==NULL){
      continue;
    }
    atomName = BondGetToName(pCurBond);
    if( atomName[0]!='+' && atomName[0]!='-' && ResidueGetAtomByName(pThis, atomName)==NULL){
      continue;
    }
    BondSetAdd(&newBonds, BondGetFromName(pCurBond), BondGetToName(pCurBond), BondGetType(pCurBond));
  }
  BondSetCopy(&pThis->bonds, &newBonds);
  BondSetDestroy(&newBonds);
  ResidueTopologyDestroy(&patchResiTopo);

  // remove the useless bonds
  char deleteBondWithPreviousOrNextResidue = '+';
  for(int i=0;i<BondSetGetCount(&pThis->bonds);i++){
    Bond* pCurBond = BondSetGet(&pThis->bonds, i);
    char* fromAtomName = BondGetFromName(pCurBond);
    char* toAtomName = BondGetToName(pCurBond);
    if( fromAtomName[0] == deleteBondWithPreviousOrNextResidue || toAtomName[0] == deleteBondWithPreviousOrNextResidue){
      BondSetRemove(&pThis->bonds, fromAtomName, toAtomName);
      break;
    }
  }
  return Success;
}

int ResiduePatchNTERorCTER(Residue* pThis, char* NTERorCTER,AtomParamsSet* pAtomParam,ResiTopoSet* pTopos){
  int result;
  char deleteBondWithPreviousOrNextResidue = ' ';
  if(strcmp(NTERorCTER, "NTER")==0){
    if(strcmp(ResidueGetName(pThis), "GLY")==0){
      result = ResiduePatch(pThis, "GLYP", pAtomParam, pTopos);
    }
    else if(strcmp(ResidueGetName(pThis), "PRO")==0){
      result = ResiduePatch(pThis, "PROP", pAtomParam, pTopos);
    }
    else{
      result = ResiduePatch(pThis, "NTER", pAtomParam, pTopos);
    }
    deleteBondWithPreviousOrNextResidue = '-';
  }
  else if(strcmp(NTERorCTER, "CTER")==0){
    result = ResiduePatch(pThis, "CTER", pAtomParam, pTopos);
    deleteBondWithPreviousOrNextResidue = '+';
  }
  else{
    result = ValueError;
  }

  if(FAILED(result)){
    return result;
  }
  else{
    Bond* pCurBond;
    char* fromAtomName;
    char* toAtomName;
    int i;
    for(i=0;i<BondSetGetCount(&pThis->bonds);i++){
      pCurBond = BondSetGet(&pThis->bonds, i);
      fromAtomName = BondGetFromName(pCurBond);
      toAtomName = BondGetToName(pCurBond);
      if( fromAtomName[0] == deleteBondWithPreviousOrNextResidue || toAtomName[0]   == deleteBondWithPreviousOrNextResidue){
        BondSetRemove(&pThis->bonds, fromAtomName, toAtomName);
        break;
      }
    }
  }
  return Success;
}

StringArray* ResidueGetPatchingHistory(Residue* pThis){
  return &pThis->patches;
}


int ResidueCalcAtomXYZ(Residue* pThis, ResiTopoSet* pResiTopos,Residue* pPrevResi, Residue* pNextResi, char* atomName, XYZ* pDestXYZ){
  // find IC from patch topology orderly;
  BOOL icFound = FALSE;
  CharmmIC ic;
  ResidueTopology topo;
  CharmmICCreate(&ic);
  ResidueTopologyCreate(&topo);
  for(int i=0;i<StringArrayGetCount(&pThis->patches);i++){
    ResiTopoSetGet(pResiTopos, StringArrayGet(&pThis->patches, i), &topo);
    if(!FAILED(ResidueTopologyFindCharmmIC(&topo, atomName, &ic))){
      icFound = TRUE;
      break;
    }          
  }
  // find IC in the residue topology;
  if(icFound == FALSE){
    ResiTopoSetGet(pResiTopos, ResidueGetName(pThis), &topo);
    if(FAILED(ResidueTopologyFindCharmmIC(&topo, atomName, &ic))){
      // if the IC is still not found in the residue topology, return an error;
      int result = DataNotExistError;
      char errMsg[MAX_LENGTH_ERR_MSG+1];
      sprintf(errMsg, "in file %s function %s() line %d, cannot find charmm IC of atom '%s' in residue %s %d %s", __FILE__, __FUNCTION__, __LINE__,
        atomName, ResidueGetName(pThis), ResidueGetPosInChain(pThis), ResidueGetChainName(pThis));
      TraceError(errMsg, result);
      return result;
    }
  }

  char* namesOfAtomABC[3];
  XYZ xyzsOfAtomABC[3];
  namesOfAtomABC[0] = CharmmICGetAtomA(&ic);
  namesOfAtomABC[1] = CharmmICGetAtomB(&ic);
  namesOfAtomABC[2] = CharmmICGetAtomC(&ic);
  // find atom coordinate from the current residue, preceding residue and next residue;
  for(int i=0;i<3;i++){
    if( namesOfAtomABC[i][0] == '-' && pPrevResi!=NULL && !FAILED(ResidueGetAtomXYZ(pPrevResi, namesOfAtomABC[i]+1, &xyzsOfAtomABC[i]))){
      continue;
    }
    else if( namesOfAtomABC[i][0] == '+' && pNextResi!=NULL && !FAILED(ResidueGetAtomXYZ(pNextResi, namesOfAtomABC[i]+1, &xyzsOfAtomABC[i]))){
      continue;
    }
    else if( namesOfAtomABC[i][0] != '+' && namesOfAtomABC[i][0] != '-' && !FAILED(ResidueGetAtomXYZ(pThis, namesOfAtomABC[i], &xyzsOfAtomABC[i]))){
      continue;
    }
    else{
      ResidueTopologyDestroy(&topo);
      CharmmICDestroy(&ic);
      return DataNotExistError;
    }
  }

  GetFourthAtom(&xyzsOfAtomABC[0], &xyzsOfAtomABC[1], &xyzsOfAtomABC[2], CharmmICGetICParams(&ic), pDestXYZ);
  ResidueTopologyDestroy(&topo);
  CharmmICDestroy(&ic);
  return Success;
}


//int ResidueCalcAtomXYZ(Residue* pThis, ResiTopoSet* pResiTopos,Residue* pPrevResi, Residue* pNextResi, char* atomName, XYZ* pDestXYZ)
//{
//  int i;
//  char* namesOfAtomABC[3];
//  XYZ xyzsOfAtomABC[3];
//  BOOL icFound;
//
//  //CharmmIC ic;
//  //ResidueTopology topo;
//  //CharmmICCreate(&ic);
//  //ResidueTopologyCreate(&topo);
//
//  icFound = FALSE;
//  int targetTopIndex = -1;
//  int targetIcIndex = -1;
//  // find IC from patch topology orderly;
//  for(i=0;i<StringArrayGetCount(&pThis->patches);i++){
//    ResiTopoCollectionGetIndex(pResiTopos, StringArrayGet(&pThis->patches, i), &targetTopIndex);
//    //ResiTopoCollectionGet(pResiTopos, StringArrayGet(&pThis->patches, i), &topo);
//    if(!FAILED(ResidueTopologyFindCharmmICIndex(&pResiTopos->topos[targetTopIndex], atomName, &targetIcIndex))){
//      icFound = TRUE;
//      break;
//    }          
//  }
//  // find IC in the residue topology;
//  if(icFound == FALSE){
//    ResiTopoCollectionGetIndex(pResiTopos, ResidueGetName(pThis), &targetTopIndex);
//    //ResiTopoCollectionGet(pResiTopos, ResidueGetName(pThis), &topo);
//    if(FAILED(ResidueTopologyFindCharmmICIndex(&pResiTopos->topos[targetTopIndex], atomName, &targetIcIndex))){
//      // if the IC is still not found in the residue topology, return an error;
//      int result = DataNotExistError;
//      char errMsg[MAX_LENGTH_ERR_MSG+1];
//      sprintf(errMsg, "in file %s function %s() line %d, cannot find charmm IC of atom '%s' in residue %s %d %s", __FILE__, __FUNCTION__, __LINE__,
//        atomName, ResidueGetName(pThis), ResidueGetPosInChain(pThis), ResidueGetChainName(pThis));
//      TraceError(errMsg, result);
//      return result;
//    }
//  }
//
//  namesOfAtomABC[0] = CharmmICGetAtomA(&pResiTopos->topos[targetTopIndex].ics[targetIcIndex]);
//  namesOfAtomABC[1] = CharmmICGetAtomB(&pResiTopos->topos[targetTopIndex].ics[targetIcIndex]);
//  namesOfAtomABC[2] = CharmmICGetAtomC(&pResiTopos->topos[targetTopIndex].ics[targetIcIndex]);
//
//  // find atom coordinate from the current residue, preceding residue and next residue;
//  for(i=0;i<3;i++){
//    if( namesOfAtomABC[i][0] == '-' && pPrevResi!=NULL && !FAILED(ResidueGetAtomXYZ(pPrevResi, namesOfAtomABC[i]+1, &xyzsOfAtomABC[i]))){
//      continue;
//    }
//    else if( namesOfAtomABC[i][0] == '+' && pNextResi!=NULL && !FAILED(ResidueGetAtomXYZ(pNextResi, namesOfAtomABC[i]+1, &xyzsOfAtomABC[i]))){
//      continue;
//    }
//    else if( namesOfAtomABC[i][0] != '+' && namesOfAtomABC[i][0] != '-' && !FAILED(ResidueGetAtomXYZ(pThis, namesOfAtomABC[i], &xyzsOfAtomABC[i]))){
//      continue;
//    }
//    else{
//      //ResidueTopologyDestroy(&topo);
//      //CharmmICDestroy(&ic);
//      return DataNotExistError;
//    }
//  }
//
//  GetFourthAtom(&xyzsOfAtomABC[0], &xyzsOfAtomABC[1], &xyzsOfAtomABC[2], CharmmICGetICParams(&pResiTopos->topos[targetTopIndex].ics[targetIcIndex]), pDestXYZ);
//  //ResidueTopologyDestroy(&topo);
//  //CharmmICDestroy(&ic);
//  return Success;
//}


int ResidueCalcAllAtomXYZ(Residue* pThis, ResiTopoSet* pResiTopos, Residue* pPrevResi, Residue* pNextResi){
  BOOL allAtomsXYZAreValid = FALSE;
  BOOL done = FALSE;
  while(!done){
    done = TRUE;
    allAtomsXYZAreValid = TRUE;
    for(int i=0;i<ResidueGetAtomCount(pThis);i++){
      XYZ newXYZ;
      Atom* pCurAtom = ResidueGetAtom(pThis, i);
      if(pCurAtom->isXyzValid) continue;
      int result = ResidueCalcAtomXYZ(pThis, pResiTopos, pPrevResi, pNextResi, AtomGetName(pCurAtom), &newXYZ);
      if(FAILED(result)){
        allAtomsXYZAreValid = FALSE;
        continue;
      }
      else{
        pCurAtom->xyz = newXYZ;
        pCurAtom->isXyzValid = TRUE;
        done = FALSE; // New atom XYZ has been calculated in this 'while' loop, go on and try to find more
      }
    }
  }

  for(int i=0;i<ResidueGetAtomCount(pThis);i++){
    Atom* pAtom=ResidueGetAtom(pThis,i);
    if(pAtom->isXyzValid==FALSE){
      printf("coordinate of atom %s on residue %s %s %d is invalid, please check\n",AtomGetName(pAtom),ResidueGetChainName(pThis),ResidueGetName(pThis),ResidueGetPosInChain(pThis));
    }
  }

  if(!allAtomsXYZAreValid){
    //int result = Warning;
    int result=ValueError;
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    sprintf(errMsg, "in file %s function %s() line %d, not all atoms' XYZ can be calculated for residue %s in chain %s %d", __FILE__, __FUNCTION__, __LINE__, ResidueGetName(pThis), ResidueGetChainName(pThis), ResidueGetPosInChain(pThis));
    TraceError(errMsg, result);
    return result;
  }

  return Success;
}

int ResidueShowAtomParameter(Residue* pThis){
  for(int i=0; i<ResidueGetAtomCount(pThis); ++i){
    Atom* pAtom=ResidueGetAtom(pThis,i);
    AtomShowAtomParameter(pAtom);
  }
  return Success;
}

int ResidueShowBondInformation(Residue* pThis){
  BondSetShow(ResidueGetBonds(pThis));
  return Success;
}

double ResidueAndResidueSidechainRMSD(Residue* pThis, Residue* pOther){
  double rmsd=0.0;
  int count=0;
  if(strcmp(ResidueGetName(pThis),ResidueGetName(pOther))==0 ||
    (strcmp(ResidueGetName(pThis),"HSD")==0 && strcmp(ResidueGetName(pOther),"HSE")==0) ||
    (strcmp(ResidueGetName(pThis),"HSE")==0 && strcmp(ResidueGetName(pOther),"HSD")==0)){
    for(int i=0; i<ResidueGetAtomCount(pThis); ++i){
      Atom* pAtom1=ResidueGetAtom(pThis,i);
      if(AtomIsHydrogen(pAtom1)||pAtom1->isBBAtom||strcmp(pAtom1->name,"CB")==0) continue;
      Atom* pAtom2=ResidueGetAtomByName(pOther,AtomGetName(pAtom1));
      rmsd+=(pAtom1->xyz.X-pAtom2->xyz.X)*(pAtom1->xyz.X-pAtom2->xyz.X)+(pAtom1->xyz.Y-pAtom2->xyz.Y)*(pAtom1->xyz.Y-pAtom2->xyz.Y)+(pAtom1->xyz.Z-pAtom2->xyz.Z)*(pAtom1->xyz.Z-pAtom2->xyz.Z);
      count++;
    }
    if(count>0) return sqrt(rmsd/count);
    else return 0.0;
  }
  //when residues are different, return a big value
  else{
    return 1e8;
  }
}

BOOL LigandResidueNameConflictWithAminoAcid(char* ligname){
  BOOL result=FALSE;
  char resiname[][MAX_LENGTH_RESIDUE_NAME+1]={
    "ALA","CYS","ASP","GLU","PHE","GLY","HIS","HSE","HSD","ILE","LYS","LEU",
    "MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"};
  for(int i=0; i<22; i++){
    if(strcmp(ligname,resiname[i])==0){
      result=TRUE;
      break;
    }
  }
  return result;
}


int ResidueCheckAtomCoordinateValidity(Residue* pThis){
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  for(int i=0;i<ResidueGetAtomCount(pThis);i++){
    Atom* pAtom=ResidueGetAtom(pThis,i);
    if(AtomIsHydrogen(pAtom)==TRUE) continue;
    if(pAtom->isBBAtom==TRUE && pAtom->isXyzValid==FALSE){
      sprintf(errMsg,"in file %s function %s line %d, the coordinate of backbone atom %s on residue %s%d%s is invalid",
        __FILE__,__FUNCTION__,__LINE__,AtomGetName(pAtom),ResidueGetChainName(pThis),ResidueGetPosInChain(pThis),ResidueGetName(pThis));
      TraceError(errMsg,ValueError);
    }
    else if(pAtom->isBBAtom==FALSE && pAtom->isXyzValid==FALSE){
      sprintf(errMsg,"in file %s function %s line %d, the coordinate of sidechain atom %s on residue %s%d%s is invalid",
        __FILE__,__FUNCTION__,__LINE__,AtomGetName(pAtom),ResidueGetChainName(pThis),ResidueGetPosInChain(pThis),ResidueGetName(pThis));
      TraceError(errMsg,Warning);
    }
  }
  return Success;
}