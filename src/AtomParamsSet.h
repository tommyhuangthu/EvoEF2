///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#ifndef ATOMPARAMSSET_H
#define ATOMPARAMSSET_H

#include "Atom.h"

typedef struct _AtomParamsSet{
  StringArray residueNames;
  Atom** atoms;
  int*   atomCount;
} AtomParamsSet;

int AtomParamsSetCreate(AtomParamsSet* pThis);
int AtomParamsSetDestroy(AtomParamsSet* pThis);
int AtomParamsSetAdd(AtomParamsSet* pThis, char* residueName, Atom* pNewAtom);
int AtomParamsSetAddFromFile(AtomParamsSet* pThis, char* filepath);
int AtomParamsSetGetResidueCount(AtomParamsSet* pThis);
int AtomParamsSetGetResidueName(AtomParamsSet* pThis, int index, char* residueName);
int AtomParamsSetGetAtomCount(AtomParamsSet* pThis, char* residueName, int* pCount);
int AtomParamsSetGetAtomParam(AtomParamsSet* pThis, char* residueName, int index, Atom* pDestAtom);
int AtomParamsSetGetAtomParamByName(AtomParamsSet* pThis, char* residueName, char* atomName, Atom* pDestAtom);
int AtomParameterRead(AtomParamsSet* pAtomParam, char* filePath);

#endif