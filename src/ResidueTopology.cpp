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

#include "ResidueTopology.h"
#include <string.h>

int BondCreate(Bond* pThis){
  strcpy(pThis->atomFromName, "");
  strcpy(pThis->atomToName, "");
  pThis->type = Type_Bond_None;
  return Success;
}

int BondDestroy(Bond* pThis){
  strcpy(pThis->atomFromName, "");
  strcpy(pThis->atomToName, "");
  pThis->type = Type_Bond_None;
  return Success;
}

int BondCopy(Bond* pThis, Bond* pOther){
  *pThis = *pOther;
  return Success;
}

char* BondGetFromName(Bond* pThis){
  return pThis->atomFromName;
}

int BondSetFromName(Bond* pThis, char* from){
  if(strlen(from)>MAX_LENGTH_ATOM_NAME){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s function %s() line %d",__FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->atomFromName, from);
  return Success;
}

char* BondGetToName(Bond* pThis){
  return pThis->atomToName;
}

int BondSetToName(Bond* pThis, char* to){
  if(strlen(to)>MAX_LENGTH_ATOM_NAME){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->atomToName, to);
  return Success;
}

Type_Bond BondGetType(Bond* pThis){
  return pThis->type;
}

int BondSetType(Bond* pThis, Type_Bond newType){
  if(newType<Type_Bond_Single || newType>Type_Bond_None){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s function %s() line %d",__FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  pThis->type = newType;
  return Success;
}

int BondShow(Bond* pThis){
  char bondType[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  switch(BondGetType(pThis)){
  case Type_Bond_None:
    strcpy(bondType, "NONE  ");break;
  case Type_Bond_Single:
    strcpy(bondType, "BOND  ");break;
  case Type_Bond_double:
    strcpy(bondType, "DOUBLE");break;
  case Type_Bond_Triple:
    strcpy(bondType, "TRIPLE");break;
  default:
    strcpy(bondType, "");break;
  }
  printf("%s  %s  %s", bondType, BondGetFromName(pThis),BondGetToName(pThis));
  return Success;
}





int BondSetCreate(BondSet* pThis){
  pThis->count = 0;
  pThis->bonds = NULL;
  return Success;
}

int BondSetDestroy(BondSet* pThis){
  for(int i=0;i<pThis->count;i++){
    BondDestroy(&pThis->bonds[i]);
  }
  free(pThis->bonds);
  pThis->bonds = NULL;
  pThis->count = 0;
  return Success;
}

int BondSetCopy(BondSet* pThis, BondSet* pOther){
  BondSetDestroy(pThis);
  BondSetCreate(pThis);
  pThis->count = pOther->count;
  pThis->bonds = (Bond*)malloc(sizeof(Bond)*pOther->count);
  for(int i=0;i<pThis->count;i++){
    BondCreate(&pThis->bonds[i]);
    BondCopy(&pThis->bonds[i], &pOther->bonds[i]);
  }
  return Success;
}

int BondSetGetCount(BondSet* pThis){
  return pThis->count;
}

int BondSetAdd(BondSet* pThis, char* atom1, char* atom2, Type_Bond bondType){
  if(BondSetFind(pThis, atom1, atom2) == Type_Bond_None){
    Bond* pNewlyAddedBond;
    (pThis->count)++;
    pThis->bonds = (Bond*)realloc(pThis->bonds, sizeof(Bond)*pThis->count);
    pNewlyAddedBond = &pThis->bonds[pThis->count-1];
    BondCreate(pNewlyAddedBond);
    if(FAILED(BondSetFromName(pNewlyAddedBond, atom1)) || FAILED(BondSetToName(pNewlyAddedBond, atom2)) || FAILED(BondSetType(pNewlyAddedBond, bondType))){
      return ValueError;
    }
  }
  return Success;
}

int BondSetRemove(BondSet* pThis, char* atom1, char* atom2){
  for(int i=0;i<pThis->count;i++){
    Bond* pCurBond = &pThis->bonds[i];
    if( (strcmp(BondGetFromName(pCurBond), atom1)==0 && strcmp(BondGetToName(pCurBond), atom2)==0) || (strcmp(BondGetFromName(pCurBond), atom2)==0 && strcmp(BondGetToName(pCurBond), atom1)==0) ){
      // Remove this bond, swap this bond with the last one in the set, and decrease the bond counter
      // Because these is at least one bond when reaching here, the operations are safe
      BondCopy(pCurBond, &pThis->bonds[pThis->count-1]); 
      BondDestroy(&pThis->bonds[pThis->count-1]);
      (pThis->count)--;
      return Success;
    }
  }
  return DataNotExistError;
}

Type_Bond BondSetFind(BondSet* pThis, char* atom1, char* atom2){
  for(int i=0;i<pThis->count;i++){
    Bond* pCurBond = &pThis->bonds[i];
    if( (strcmp(BondGetFromName(pCurBond), atom1)==0 && strcmp(BondGetToName(pCurBond), atom2)==0) || (strcmp(BondGetFromName(pCurBond), atom2)==0 && strcmp(BondGetToName(pCurBond), atom1)==0) ){
      return BondGetType(pCurBond);
    }
  }
  return Type_Bond_None;
}

Bond* BondSetGet(BondSet* pThis, int index){
  if(index<0 || index>=pThis->count){
    return NULL;
  }
  return &pThis->bonds[index];
}

int BondSetShow(BondSet* pThis){
  for(int i=0;i<pThis->count;i++){
    BondShow(&pThis->bonds[i]);
    printf("\n");
  }
  return Success;
}


//------------Internal Coordinates-----------------------

int CharmmICCreate(CharmmIC* pThis){
  for(int i=0;i<4;i++)
    strcpy(pThis->atomNames[i], "");
  return Success;
}

int CharmmICCreateFromStringArray(CharmmIC* pThis, StringArray* params){
  CharmmICCreate(pThis);
  if(StringArrayGetCount(params) != 10){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s function %s() line %d",  __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }

  for(int i=0;i<4;i++){
    if(strlen(StringArrayGet(params, i+1))>MAX_LENGTH_ATOM_NAME){
      char errMsg[MAX_LENGTH_ERR_MSG+1];
      int errorCode = ValueError;
      sprintf(errMsg, "in file %s function %s() line %d, atom name is too long", __FILE__, __FUNCTION__, __LINE__);
      TraceError(errMsg, errorCode);
      return errorCode;
    }
    strcpy(pThis->atomNames[i], StringArrayGet(params, i+1));
  }

  for(int i=0;i<5;i++){
    pThis->icParam[i] = atof(StringArrayGet(params, i+5));
    if(i>=1 && i<=3){
      // covert degree to arc
      pThis->icParam[i] = DegToRad(pThis->icParam[i]);
    }  
  }

  char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  strcpy(buffer, CharmmICGetAtomC(pThis));
  // deal with improper torsion angle
  if(buffer[0]=='*'){
    pThis->torsionProperFlag = FALSE;
    strcpy(pThis->atomNames[2], buffer+1);
  }
  else{
    pThis->torsionProperFlag = TRUE;
  }
  return Success;
}

int CharmmICDestroy(CharmmIC* pThis){
  for(int i=0;i<4;i++)
    strcpy(pThis->atomNames[i], "");
  return Success;
}

int CharmmICCopy(CharmmIC* pThis, CharmmIC* pOther){
  *pThis = *pOther;
  return Success;
}

char* CharmmICGetAtomA(CharmmIC* pThis){
  return pThis->atomNames[0];
}

char* CharmmICGetAtomB(CharmmIC* pThis){
  return pThis->atomNames[1];
}

char* CharmmICGetAtomC(CharmmIC* pThis){
  return pThis->atomNames[2];
}

char* CharmmICGetAtomD(CharmmIC* pThis){
  return pThis->atomNames[3];
}

double* CharmmICGetICParams(CharmmIC* pThis){
  return pThis->icParam;
}

int CharmmICSetTorsionAngle(CharmmIC* pThis, double newTorsionAngle){
  pThis->icParam[2] = newTorsionAngle;
  return Success;
}

BOOL CharmmICGetTorsionProperFlag(CharmmIC* pThis){
  return pThis->torsionProperFlag;
}

int CharmmICCalcXYZ(CharmmIC* pThis, AtomArray* atomArray, XYZ* pDestXYZ){
  XYZ* pAtomsABC[3];
  for(int i=0;i<3;i++){
    Atom* pAtom = AtomArrayGetByName(atomArray, pThis->atomNames[i]);
    if(pAtom == NULL){
      return DataNotExistError;
    }

    if(pAtom->isXyzValid==FALSE){
      return DataNotExistError;
    }

    pAtomsABC[i] = &pAtom->xyz;
  }
  return GetFourthAtom(pAtomsABC[0], pAtomsABC[1], pAtomsABC[2], pThis->icParam, pDestXYZ);
}

int CharmmICShow(CharmmIC* pThis){
  printf("IC  %-7.7s %-7.7s %c%-7.7s %-7.7s  ", CharmmICGetAtomA(pThis), CharmmICGetAtomB(pThis),CharmmICGetTorsionProperFlag(pThis) ? ' ' : '*',CharmmICGetAtomC(pThis), CharmmICGetAtomD(pThis));
  printf("%8.3f ", pThis->icParam[0]);
  for(int i=1;i<4;i++){
    printf("%8.3f ", RadToDeg(pThis->icParam[i]));
  }
  printf("%8.3f", pThis->icParam[4]);
  return Success;
}


int ResidueTopologyCreate(ResidueTopology* pThis){
  StringArrayCreate(&pThis->atoms);
  StringArrayCreate(&pThis->deletes);
  BondSetCreate(&pThis->bonds);
  pThis->icCount = 0;
  pThis->ics = NULL;
  return Success;
}

int ResidueTopologyCreateFromFileReader(ResidueTopology* pThis, FileReader* pFileReader){
  BOOL residueNameAlreadySet = FALSE;
  ResidueTopologyCreate(pThis);
  while(!FileReaderEndOfFile(pFileReader)){
    char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    char* keyword;
    StringArray wordsInLine;
    StringArrayCreate(&wordsInLine);

    FileReaderGetNextLine(pFileReader, line);
    StringArraySplitString(&wordsInLine, line, ' ');
    keyword = StringArrayGet(&wordsInLine, 0);

    if(strcmp(keyword, "RESI")==0 || strcmp(keyword, "PRES")==0){
      if(residueNameAlreadySet){
        FileReaderSetCurrentPos(pFileReader, FileReaderGetCurrentPos(pFileReader)-1);
        StringArrayDestroy(&wordsInLine);
        return Success;
      }
      else{
        char* residueName = StringArrayGet(&wordsInLine, 1);
        if(residueName == NULL || FAILED(ResidueTopologySetName(pThis, residueName))){
          char errMsg[MAX_LENGTH_ERR_MSG+1];
          int errorCode = ValueError;
          sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
          TraceError(errMsg, errorCode);
          return errorCode;
        }
        residueNameAlreadySet = TRUE;
      }
    }
    else if(strcmp(keyword, "ATOM")==0){
      char* atomName = StringArrayGet(&wordsInLine, 1);
      if(atomName == NULL || FAILED(StringArrayAppend(&pThis->atoms, atomName))){
        char errMsg[MAX_LENGTH_ERR_MSG+1];
        int errorCode = ValueError;
        sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
        TraceError(errMsg, errorCode);
        return errorCode;
      }
    }
    else if(strcmp(keyword, "DELETE")==0){
      char* deleteAtomName = StringArrayGet(&wordsInLine, 2);
      if(deleteAtomName == NULL || FAILED(StringArrayAppend(&pThis->deletes, deleteAtomName))){
        char errMsg[MAX_LENGTH_ERR_MSG+1];
        int errorCode = ValueError;
        sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
        TraceError(errMsg, errorCode);
        return errorCode;
      }
    }
    else if( strcmp(keyword, "BOND")==0||strcmp(keyword, "DOUBLE")==0||strcmp(keyword, "TRIPLE")==0){
      for(int i=1;i<StringArrayGetCount(&wordsInLine);i+=2){
        char* atomFromName = StringArrayGet(&wordsInLine, i);
        char* atomToName   = StringArrayGet(&wordsInLine, i+1);
        Type_Bond bondType;

        if(atomFromName == NULL || atomToName == NULL){
          char errMsg[MAX_LENGTH_ERR_MSG+1];
          int errorCode = ValueError;
          sprintf(errMsg, "in file %s function %s() line %d, atom %s or %s cannot be found", __FILE__, __FUNCTION__, __LINE__, atomFromName, atomToName);
          TraceError(errMsg, errorCode);
          return errorCode;
        }
        if(strcmp(keyword, "BOND")==0){
          bondType = Type_Bond_Single;
        }
        else if(strcmp(keyword, "DOUBLE")==0){
          bondType = Type_Bond_double;
        }
        else if(strcmp(keyword, "TRIPLE")==0){ 
          bondType = Type_Bond_Triple;
        }

        if(FAILED(BondSetAdd(&pThis->bonds, atomFromName, atomToName, bondType))){
          char errMsg[MAX_LENGTH_ERR_MSG+1];
          int errorCode = ValueError;
          sprintf(errMsg, "in file %s function %s() line %d, cannot add bond to set",  __FILE__, __FUNCTION__, __LINE__);
          TraceError(errMsg, errorCode);
          return errorCode;
        }
      }
    }
    else if(strcmp(keyword, "IC")==0){
      CharmmIC newIC;
      if(FAILED(CharmmICCreateFromStringArray(&newIC, &wordsInLine)) || FAILED(ResidueTopologyAddCharmmIC(pThis, &newIC)))
      {
        char errMsg[MAX_LENGTH_ERR_MSG+1];
        int errorCode = ValueError;
        sprintf(errMsg, "in file %s function %s() line %d, CHARMMIC creating error", __FILE__, __FUNCTION__, __LINE__);
        TraceError(errMsg, errorCode);
        return errorCode;
      }
      CharmmICDestroy(&newIC);
    }
    else{
      // neglect other unrecognized keyword in topology file
    }
    StringArrayDestroy(&wordsInLine);
  }
  return Success;
}

int ResidueTopologyDestroy(ResidueTopology* pThis){
  StringArrayDestroy(&pThis->atoms);
  StringArrayDestroy(&pThis->deletes);
  BondSetDestroy(&pThis->bonds);

  for(int i=0;i<pThis->icCount;i++){
    CharmmICDestroy(&pThis->ics[i]);
  }
  free(pThis->ics);
  pThis->ics = NULL;
  pThis->icCount = 0;
  return Success;
}

int ResidueTopologyCopy(ResidueTopology* pThis, ResidueTopology* pOther){
  ResidueTopologyDestroy(pThis);
  ResidueTopologyCreate(pThis);
  ResidueTopologySetName(pThis, ResidueTopologyGetName(pOther));
  StringArrayCopy(&pThis->atoms, &pOther->atoms);
  StringArrayCopy(&pThis->deletes, &pOther->deletes);
  BondSetCopy(&pThis->bonds, &pOther->bonds);
  pThis->icCount = pOther->icCount;
  pThis->ics = (CharmmIC*)malloc(sizeof(CharmmIC)*pThis->icCount);
  for(int i=0;i<pThis->icCount;i++){
    CharmmICCreate(&pThis->ics[i]);
    CharmmICCopy(&pThis->ics[i], &pOther->ics[i]);
  }
  return Success;
}

char* ResidueTopologyGetName(ResidueTopology* pThis){
  return pThis->residueName;
}

int ResidueTopologySetName(ResidueTopology* pThis, char* newName){
  if(strlen(newName)>MAX_LENGTH_RESIDUE_NAME){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s function %s() line %d, residue name is too long", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->residueName, newName);
  return Success;
}

StringArray* ResidueTopologyGetAtoms(ResidueTopology* pThis){
  return &pThis->atoms;
}

StringArray* ResidueTopologyGetDeletes(ResidueTopology* pThis){
  return &pThis->deletes;
}

BondSet* ResidueTopologyGetBonds(ResidueTopology* pThis){
  return &pThis->bonds;
}

int ResidueTopologyGetCharmmICCount(ResidueTopology* pThis){
  return pThis->icCount;
}

int ResidueTopologyGetCharmmIC(ResidueTopology* pThis, int index, CharmmIC* pDestIC){
  if(index<0 || index>=pThis->icCount)
    return IndexError;
  CharmmICCopy(pDestIC, &pThis->ics[index]);
  return Success;
}

int ResidueTopologyFindCharmmIC(ResidueTopology* pThis, char* atomDName, CharmmIC* pDestIC){
  int i;
  for(i=0;i<pThis->icCount;i++){
    if(strcmp(CharmmICGetAtomD(&pThis->ics[i]), atomDName)==0){
      CharmmICCopy(pDestIC, &pThis->ics[i]);
      return Success;
    }
  }
  return DataNotExistError;
}

int ResidueTopologyFindCharmmICIndex(ResidueTopology* pThis, char* atomDName, int *index){
  for(int i=0;i<pThis->icCount;i++){
    if(strcmp(CharmmICGetAtomD(&pThis->ics[i]), atomDName)==0){
      *index = i;
      return Success;
    }
  }
  return DataNotExistError;
}

int ResidueTopologyAddCharmmIC(ResidueTopology* pThis, CharmmIC* pNewIC){
  // Will not check if it already exists
  int newICCount = pThis->icCount+1;
  pThis->ics = (CharmmIC*)realloc(pThis->ics, sizeof(CharmmIC)*newICCount);
  CharmmICCreate(&pThis->ics[newICCount-1]);
  CharmmICCopy(&pThis->ics[newICCount-1], pNewIC);
  pThis->icCount = newICCount;
  return Success;
}

int ResidueTopologyShow(ResidueTopology* pThis){
  printf("RESI  %s\n", ResidueTopologyGetName(pThis));
  for(int i=0;i<StringArrayGetCount(&pThis->atoms);i++){
    printf("ATOM   %s\n", StringArrayGet(&pThis->atoms, i));
  }
  BondSetShow(&pThis->bonds);
  for(int i=0;i<pThis->icCount;i++){
    CharmmICShow(&pThis->ics[i]);
    printf("\n");
  }
  return Success;
}


int ResiTopoSetCreate(ResiTopoSet* pThis){
  pThis->count = 0;
  pThis->topos = NULL;
  return Success;
}

int ResiTopoSetDestroy(ResiTopoSet* pThis){
  for(int i=0;i<pThis->count;i++){
    ResidueTopologyDestroy(&pThis->topos[i]);
  }
  free(pThis->topos);
  pThis->topos = NULL;
  pThis->count = 0;
  return Success;
}

int ResiTopoSetCopy(ResiTopoSet* pThis, ResiTopoSet* pOther){
  ResiTopoSetDestroy(pThis);
  pThis->count = pOther->count;
  pThis->topos = (ResidueTopology*)malloc(sizeof(ResidueTopology)*pThis->count);
  for(int i=0;i<pThis->count;i++){
    ResidueTopologyCreate(&pThis->topos[i]);
    ResidueTopologyCopy(&pThis->topos[i], &pOther->topos[i]);
  }
  return Success;
}

int ResiTopoSetGet(ResiTopoSet* pThis, char* resiName, ResidueTopology* pDestTopo){
  for(int i=0;i<pThis->count;i++){
    if(strcmp(ResidueTopologyGetName(&pThis->topos[i]), resiName)==0){
      ResidueTopologyCopy(pDestTopo, &pThis->topos[i]);
      return Success;
    }
  }
  return DataNotExistError;
}

int ResiTopoCollectionGetIndex(ResiTopoSet* pThis, char* resiName, int *index){
  for(int i=0;i<pThis->count;i++){
    if(strcmp(ResidueTopologyGetName(&pThis->topos[i]), resiName)==0){
      *index = i;
      return Success;
    }
  }
  return DataNotExistError;
}

int ResiTopoSetAdd(ResiTopoSet* pThis, ResidueTopology* pNewTopo){
  // If already exist, replace the original
  for(int i=0;i<pThis->count;i++){
    if(strcmp(ResidueTopologyGetName(&pThis->topos[i]), ResidueTopologyGetName(pNewTopo))==0){
      ResidueTopologyCopy(&pThis->topos[i], pNewTopo);
      return Success;
    }
  }

  int newCount = pThis->count+1;
  pThis->topos = (ResidueTopology*)realloc(pThis->topos, sizeof(ResidueTopology)*newCount);
  ResidueTopologyCreate(&pThis->topos[newCount-1]);
  ResidueTopologyCopy(&pThis->topos[newCount-1], pNewTopo);
  pThis->count = newCount;
  return Success;
}

int ResiTopoSetAddFromFile(ResiTopoSet* pThis, char* filepath){
  FileReader file;
  int result = FileReaderCreate(&file, filepath);
  if(FAILED(result)){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, result);
    return result;
  }

  while(!FileReaderEndOfFile(&file)){
    ResidueTopology newTopo;
    result = ResidueTopologyCreateFromFileReader(&newTopo, &file);
    if(FAILED(result)){
      char errMsg[MAX_LENGTH_ERR_MSG+1];
      sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
      TraceError(errMsg, result);
      return result;
    }
    ResiTopoSetAdd(pThis, &newTopo);
    ResidueTopologyDestroy(&newTopo);
  }
  FileReaderDestroy(&file);
  return Success;
}


int ResiTopoSetRead(ResiTopoSet* pResiTopo, char* filePath){
  if(filePath == NULL  || FAILED(ResiTopoSetAddFromFile(pResiTopo, filePath))){
    return FormatError;
  }
  return Success;
}



int ResiTopoSetShow(ResiTopoSet* pThis){
  for(int i=0;i<pThis->count;i++){
    printf("Residue No. %d\n", i);
    ResidueTopologyShow(&pThis->topos[i]);
  }
  return Success;
}
