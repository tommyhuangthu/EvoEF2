///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#ifndef RESIDUE_TOPOLOGY_H
#define RESIDUE_TOPOLOGY_H

#include "Utility.h"
#include "Atom.h"

typedef enum _Type_Bond{
  Type_Bond_Single, 
  Type_Bond_double, 
  Type_Bond_Triple, 
  Type_Bond_None
}Type_Bond ;

typedef struct _Bond{
  char atomFromName[MAX_LENGTH_ATOM_NAME+1];
  char atomToName[MAX_LENGTH_ATOM_NAME+1];
  Type_Bond type;
} Bond;

int BondCreate(Bond* pThis);
int BondDestroy(Bond* pThis);
int BondCopy(Bond* pThis, Bond* pOther);
char* BondGetFromName(Bond* pThis);
int BondSetFromName(Bond* pThis, char* from);
char* BondGetToName(Bond* pThis);
int BondSetToName(Bond* pThis, char* to);
Type_Bond BondGetType(Bond* pThis);
int BondSetType(Bond* pThis, Type_Bond newType);
int BondShow(Bond* pThis);

typedef struct _BondSet{
  int count;
  Bond* bonds;
} BondSet;

int BondSetCreate(BondSet* pThis);
int BondSetDestroy(BondSet* pThis);
int BondSetCopy(BondSet* pThis, BondSet* pOther);
int BondSetGetCount(BondSet* pThis);
int BondSetAdd(BondSet* pThis, char* atom1, char* atom2, Type_Bond bondType);
int BondSetRemove(BondSet* pThis, char* atom1, char* atom2);
Type_Bond BondSetFind(BondSet* pThis, char* atom1, char* atom2);
Bond* BondSetGet(BondSet* pThis, int index);
int BondSetShow(BondSet* pThis);

typedef struct _CharmmIC{
  double icParam[5];
  char atomNames[4][MAX_LENGTH_ATOM_NAME+1];
  BOOL torsionProperFlag;                    // 1 byte
  // Proper dihedral angle (A, B, C, D): five parameters are R(AB), Theta(ABC), Phi(ABCD), Theta(BCD), R(CD)
  // Improper dihedral angle (A, B, *C, D); five parameters are R(AC), Theta(ACB), Phi(ABCD), Theta(BCD), R(CD)
}CharmmIC;

int CharmmICCreate(CharmmIC* pThis);
int CharmmICCreateFromStringArray(CharmmIC* pThis, StringArray* params);
int CharmmICDestroy(CharmmIC* pThis);
int CharmmICCopy(CharmmIC* pThis, CharmmIC* pOther);
char* CharmmICGetAtomA(CharmmIC* pThis);
char* CharmmICGetAtomB(CharmmIC* pThis);
char* CharmmICGetAtomC(CharmmIC* pThis);
char* CharmmICGetAtomD(CharmmIC* pThis);
double* CharmmICGetICParams(CharmmIC* pThis);
int CharmmICSetTorsionAngle(CharmmIC* pThis, double newTorsionAngle);
BOOL CharmmICGetTorsionProperFlag(CharmmIC* pThis);
int CharmmICCalcXYZ(CharmmIC* pThis, AtomArray* atomArray, XYZ* pDestXYZ);
int CharmmICShow(CharmmIC* pThis);

typedef struct _ResidueTopology{
  StringArray atoms;
  StringArray deletes;
  BondSet bonds;
  CharmmIC* ics;
  char residueName[MAX_LENGTH_RESIDUE_NAME+1];
  int icCount;
} ResidueTopology;

int ResidueTopologyCreate(ResidueTopology* pThis);
int ResidueTopologyCreateFromFileReader(ResidueTopology* pThis, FileReader* pFileReader);
int ResidueTopologyDestroy(ResidueTopology* pThis);
int ResidueTopologyCopy(ResidueTopology* pThis, ResidueTopology* pOther);

char* ResidueTopologyGetName(ResidueTopology* pThis);
int ResidueTopologySetName(ResidueTopology* pThis, char* newName);

StringArray* ResidueTopologyGetAtoms(ResidueTopology* pThis);
StringArray* ResidueTopologyGetDeletes(ResidueTopology* pThis);
BondSet* ResidueTopologyGetBonds(ResidueTopology* pThis);

int ResidueTopologyGetCharmmICCount(ResidueTopology* pThis);
int ResidueTopologyGetCharmmIC(ResidueTopology* pThis, int index, CharmmIC* pDestIC);
int ResidueTopologyFindCharmmIC(ResidueTopology* pThis, char* atomDName, CharmmIC* pDestIC);
int ResidueTopologyAddCharmmIC(ResidueTopology* pThis, CharmmIC* pNewIC);
int ResidueTopologyShow(ResidueTopology* pThis);



typedef struct _ResiTopoSet{
  int count;
  ResidueTopology* topos;
} ResiTopoSet;

int ResiTopoSetCreate(ResiTopoSet* pThis);
int ResiTopoSetDestroy(ResiTopoSet* pThis);
int ResiTopoSetCopy(ResiTopoSet* pThis, ResiTopoSet* pOther);
int ResiTopoSetGet(ResiTopoSet* pThis, char* resiName, ResidueTopology* pDestTopo);
int ResiTopoSetAdd(ResiTopoSet* pThis, ResidueTopology* pNewTopo);
int ResiTopoSetAddFromFile(ResiTopoSet* pThis, char* filepath);
int ResiTopoSetRead(ResiTopoSet* pResiTopo, char* filePath);
int ResiTopoSetShow(ResiTopoSet* pThis);
#endif