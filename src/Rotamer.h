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

#ifndef ROTAMER_H
#define ROTAMER_H

#include "Atom.h"
#include "Residue.h"

typedef struct _BBindRotamerLib{
  StringArray residueTypeNames;
  IntArray    rotamerCounts;
  DoubleArray** torsions;
} BBindRotamerLib;

int BBindRotamerLibCreate(BBindRotamerLib* pThis,char* rotlibFile);
int BBindRotamerLibDestroy(BBindRotamerLib* pThis);
int BBindRotamerLibGetCount(BBindRotamerLib* pThis,char* typeName);
int BBindRotamerLibGet(BBindRotamerLib* pThis,char* typeName,int index,DoubleArray* pDestTorsion);
int BBindRotamerLibShow(BBindRotamerLib* pThis);
int BBindRotamerLibTester(char* rotlibFile);


typedef struct _Rotamer{
  AtomArray atoms;
  BondSet bonds;
  XYZArray xyzs;                       
  char type[MAX_LENGTH_RESIDUE_NAME+1];
  char chainName[MAX_LENGTH_CHAIN_NAME+1];
  int  posInChain;
  double vdwBackbone;
  double vdwInternal;
  double selfenergy;
  double dunbrack;
  DoubleArray xtorsions;
} Rotamer;

int RotamerCreate(Rotamer* pThis);
int RotamerDestroy(Rotamer* pThis);
int RotamerCopy(Rotamer* pThis,Rotamer* pOther);
char* RotamerGetType(Rotamer* pThis);
int RotamerSetType(Rotamer* pThis,char* newType);
int RotamerCopyAtomXYZ(Rotamer* pThis,XYZArray* pNewXYZ);
char* RotamerGetChainName(Rotamer* pThis);
int RotamerSetChainName(Rotamer* pThis,char* newChainname);
int RotamerGetPosInChain(Rotamer* pThis);
int RotamerSetPosInChain(Rotamer* pThis,int newPosInChain);
int RotamerShow(Rotamer* pThis);
BOOL RotamerAndRotamerInSameType(Rotamer* pThis,Rotamer* pOther);
BOOL RotamerAndResidueInSameType(Rotamer* pThis,Residue* pOther);

// These functions below are valid only when Rotamer is restored
int RotamerGetAtomCount(Rotamer* pThis);
Atom* RotamerGetAtom(Rotamer* pThis,int index);
Atom* RotamerGetAtomByName(Rotamer* pThis,char* atomName);
int RotamerFindAtom(Rotamer* pThis,char* atomName,int* index);
int RotamerAddAtoms(Rotamer* pThis,AtomArray* pNewAtoms);
BondSet* RotamerGetBonds(Rotamer* pThis);

typedef enum _Type_ProteinAtomOrder{
  Type_ProteinAtomOrder_Alpha = 0,
  Type_ProteinAtomOrder_Beta,
  Type_ProteinAtomOrder_Gamme,
  Type_ProteinAtomOrder_Delta,
  Type_ProteinAtomOrder_Epsilon,
  Type_ProteinAtomOrder_Zeta,
  Type_ProteinAtomOrder_Other
} Type_ProteinAtomOrder;

int Type_ProteinAtomOrder_ToInt(Type_ProteinAtomOrder order); 
Type_ProteinAtomOrder Type_ProteinAtomOrder_FromInt(int order);
Type_ProteinAtomOrder Type_ProteinAtomOrder_JudgedByAtomName(char* atomName);

int RotamerOfProteinInitAtomsAndBonds_Charmm19(Rotamer* pThis,Residue* pResi, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int RotamerOfProteinInitAtomsAndBonds_Charmm22(Rotamer* pThis,Residue* pResi, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int RotamerOfProteinPatch(Rotamer* pThis,char* patchType, AtomParamsSet* atomParams, ResiTopoSet* resiTopo);
int RotamerOfProteinCalcXYZ(Rotamer* pThis, Residue *pResi, char* patchName, DoubleArray* torsions, ResiTopoSet* resiTopos);
int RotamerOfProteinGenerate(Rotamer* pThis,Residue* pResi,char* rotamerType, char* patchType, DoubleArray* torsions, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int RotamerShowInPDBFormat(Rotamer* pThis,char* header,char* chainName, int atomIndex,int resiIndex,BOOL showHydrogen,FILE* pFile);
int RotamerShowAtomParameter(Rotamer* pThis);


double RotamerAndRotamerSidechainRMSD(Rotamer* pThis, Rotamer* pOther);
double RotamerAndResidueSidechainRMSD(Rotamer* pThis, Residue* pOther);



typedef struct _RotamerSet{
  Rotamer* rotamers;
  Rotamer* representatives;
  int count;
  int capacity;
  int representativeCount;
} RotamerSet;

int RotamerSetCreate(RotamerSet* pThis);
int RotamerSetDestroy(RotamerSet* pThis);
int RotamerSetCopy(RotamerSet* pThis,RotamerSet* pOther);
int RotamerSetGetCount(RotamerSet* pThis);
Rotamer* RotamerSetGet(RotamerSet* pThis,int index);
Rotamer* RotamerSetGetRepresentative(RotamerSet* pThis,char* type);
int RotamerSetAdd(RotamerSet* pThis,Rotamer* pNewRotamer);
int RotamerSetOfProteinGenerate(RotamerSet* pThis,Residue* pResi, StringArray* designTypes, StringArray* patchTypes, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopo);
int RotamerSetShow(RotamerSet* pThis,FILE* pFile);

int RotamerExtract(Rotamer* pThis);
int RotamerRestore(Rotamer* pThis,RotamerSet* pRotamerSet);
int RotamerSetGetRepresentativeCount(RotamerSet* pThis);
Rotamer* RotamerSetGetRepresentativeByIndex(RotamerSet* pThis,int index);
int RotamerShowBondInformation(Rotamer* pThis);


typedef struct _RotLibPhiPsi{
  StringArray rotTypes;
  IntArray    rotamerCounts;
  DoubleArray** torsions;
  DoubleArray** deviations;
  DoubleArray* probability;
  double phi;
  double psi;
}RotLibPhiPsi;

int RotLibPhiPsiCreate(RotLibPhiPsi* pThis);
int RotLibPhiPsiDestroy(RotLibPhiPsi* pThis);
int RotLibPhiPsiGetCount(RotLibPhiPsi* pThis,char* typeName);
int RotLibPhiPsiGet(RotLibPhiPsi* pThis,char* typeName,int index,DoubleArray* pDestTorsion,double* probability);


typedef struct _BBdepRotamerLib{
  int phipsicount;
  RotLibPhiPsi* rotlibphipsis;
}BBdepRotamerLib;


int BBdepRotamerLibCreate(BBdepRotamerLib* pRotLib,char* rotlibfile);
int BBdepRotamerLibDestroy(BBdepRotamerLib* pThis);
int RotamerOfProteinGenerateByBBdepRot(Rotamer* pThis,Residue* pResi,char* rotamerType,char* patchType,DoubleArray* torsions,AtomParamsSet* atomParams,ResiTopoSet* resiTopos);
int RotamerSetOfProteinGenerateByBBdepRotLib(RotamerSet* pThis, Residue* pResi, StringArray* designTypes, StringArray* patchTypes, BBdepRotamerLib* bbrotlib,AtomParamsSet* atomParams, ResiTopoSet* resiTopo);
int RotamerCalcDunbrackEnergy(Rotamer* pThis,double probability);
BOOL RotamerAndResidueWithSimilarTorsions(Rotamer* pThis,Residue* pOther,double cutoff);
int ResidueCalcSidechainTorsion(Residue* pThis,ResiTopoSet* pResiTopos);
Rotamer* RotamerSetFindFirstGivenTypeRotamer(RotamerSet* pThis,char* type);

BOOL RotamerIsSymmetricalCheck(Rotamer* pRotamer);
int SymmetricalRotamerGenerate(Rotamer* pSymmetrical, Rotamer* pRotamer);
#endif //ROTAMER_H
