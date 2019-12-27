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

#include "AtomParamsSet.h"
#include <string.h>

int AtomParamsSetCreate(AtomParamsSet* pThis){
  StringArrayCreate(&pThis->residueNames);
  pThis->atomCount = NULL;
  pThis->atoms = NULL;
  return Success;
}

int AtomParamsSetDestroy(AtomParamsSet* pThis){
  for(int i=0;i<AtomParamsSetGetResidueCount(pThis);i++){
    for(int j=0;j<pThis->atomCount[i];j++){
      AtomDestroy(&pThis->atoms[i][j]);
    }
    free(pThis->atoms[i]);
  }
  free(pThis->atoms);
  free(pThis->atomCount);
  StringArrayDestroy(&pThis->residueNames);
  return Success;
}

int AtomParamsSetAdd(AtomParamsSet* pThis, char* residueName, Atom* pNewAtom){
  int resiIndex;
  int atomIndex;

  if(strlen(residueName)>MAX_LENGTH_RESIDUE_NAME){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s function %s() line %d, residue name is too long", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }

  int result = StringArrayFind(&pThis->residueNames, residueName, &resiIndex);
  // new residue
  if(FAILED(result)){
    int newResidueCount = AtomParamsSetGetResidueCount(pThis) + 1;
    pThis->atoms = (Atom**)realloc(pThis->atoms, sizeof(Atom*)*newResidueCount);
    pThis->atoms[newResidueCount-1] = NULL;
    StringArrayAppend(&pThis->residueNames, residueName);
    pThis->atomCount = (int*) realloc(pThis->atomCount, sizeof(int)*newResidueCount);
    pThis->atomCount[newResidueCount-1] = 0;
    resiIndex = AtomParamsSetGetResidueCount(pThis)-1;
  }

  for(atomIndex=0;atomIndex<pThis->atomCount[resiIndex];atomIndex++){
    if(strcmp(pThis->atoms[resiIndex][atomIndex].name, pNewAtom->name) == 0){
      break;
    }
  }
  // new atom
  if(atomIndex == pThis->atomCount[resiIndex]){
    int newAtomCount = pThis->atomCount[resiIndex] + 1;
    pThis->atoms[resiIndex] = (Atom*)realloc(pThis->atoms[resiIndex], sizeof(Atom)*newAtomCount);
    pThis->atomCount[resiIndex]++;
  }

  AtomCopy(&pThis->atoms[resiIndex][atomIndex], pNewAtom);
  //AtomShowParams(&pThis->atoms[resiIndex][atomIndex]);

  return Success;
}

int AtomParamsSetAddFromFile(AtomParamsSet* pThis, char* filePath){
  FileReader file;
  int result;
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  char errMsg[MAX_LENGTH_ERR_MSG+1];

  result = FileReaderCreate(&file, filePath);
  if(FAILED(result)){
    sprintf(errMsg, "in file %s function %s() line %d, when opening: \n%s", __FILE__, __FUNCTION__, __LINE__, filePath);
    TraceError(errMsg, IOError);
    return IOError;
  }

  while(!FAILED(FileReaderGetNextLine(&file, line))){
    Atom atom;
    StringArray params;

    AtomCreate(&atom);
    StringArrayCreate(&params);
    StringArraySplitString(&params, line, ' ');
    char* resiName = StringArrayGet(&params, 0);
    if( resiName == NULL || FAILED(AtomSetParamsByStringArray(&atom, &params)) || FAILED(AtomParamsSetAdd(pThis, resiName, &atom)) ){
      sprintf(errMsg, "in file %s function %s() line %d, when reading:\n%s", __FILE__, __FUNCTION__, __LINE__, line);
      TraceError(errMsg, FormatError);
      return FormatError;
    }
    StringArrayDestroy(&params);
    AtomDestroy(&atom);
  }

  FileReaderDestroy(&file);
  return Success;
}

int AtomParamsSetGetResidueCount(AtomParamsSet* pThis){
  return StringArrayGetCount(&pThis->residueNames);
}

int AtomParamsSetGetResidueName(AtomParamsSet* pThis, int index, char* residueName){
  char* name = StringArrayGet(&pThis->residueNames, index);
  if(name==NULL){
    return DataNotExistError;
  }
  else{
    strcpy(residueName, name);
    return Success;
  }
}

int AtomParamsSetGetAtomCount(AtomParamsSet* pThis,char* residueName, int* pCount){
  int resiIndex;
  int result;
  result = StringArrayFind(&pThis->residueNames, residueName, &resiIndex);
  if(FAILED(result)){
    return result;
  }
  *pCount = pThis->atomCount[resiIndex];
  return Success;
}

int AtomParamsSetGetAtomParam(AtomParamsSet* pThis,char* residueName, int index, Atom* pDestAtom){
  int resiIndex;
  int result;
  result = StringArrayFind(&pThis->residueNames, residueName, &resiIndex);
  if(FAILED(result)){
    return result;
  }
  if( index<0 || index>=pThis->atomCount[resiIndex] ){
    result = IndexError;
    return result;
  }
  AtomCopy(pDestAtom, &pThis->atoms[resiIndex][index]);
  return Success;  
}

int AtomParamsSetGetAtomParamByName(AtomParamsSet* pThis,char* residueName, char* atomName, Atom* pDestAtom){
  int resiIndex, atomIndex;
  int result = StringArrayFind(&pThis->residueNames, residueName, &resiIndex);
  if(FAILED(result)){
    return result;
  }
  for(atomIndex=0;atomIndex<pThis->atomCount[resiIndex];atomIndex++){
    if(strcmp(pThis->atoms[resiIndex][atomIndex].name, atomName)==0)
      break;
  }
  if(atomIndex==pThis->atomCount[resiIndex]){ // atom not found
    result = DataNotExistError;
    return result;
  }
  AtomCopy(pDestAtom, &pThis->atoms[resiIndex][atomIndex]);
  return Success;
}

int AtomParameterRead(AtomParamsSet* pAtomParam, char* filePath){
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  if(FAILED(AtomParamsSetAddFromFile(pAtomParam, filePath))){
    sprintf(errMsg, "in file %s function %s() line %d, when opening:\n%s",__FILE__, __FUNCTION__, __LINE__, filePath);
    TraceError(errMsg, FormatError);
    return FormatError;
  }
  return Success;
}
