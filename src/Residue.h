///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#ifndef RESIDUE_H
#define RESIDUE_H

#include "AtomParamsSet.h"
#include "ResidueTopology.h"



// enum types associated with Residue
typedef enum _Type_ResiduePolarity{
  Type_ResiduePolarity_Charged, 
  Type_ResiduePolarity_Polar, 
  Type_ResiduePolarity_NonPolar
}Type_ResiduePolarity;

typedef enum _Type_ResiduePosition{
  Type_ResiduePosition_Buried, 
  Type_ResiduePosition_Intermediate, 
  Type_ResiduePosition_Exposed
}Type_ResiduePosition;


typedef enum _Type_ResidueDesign{
  Type_ResidueDesign_Catalytic, 
  Type_ResidueDesign_PrimaryContact, 
  Type_ResidueDesign_SecondaryContact
}Type_ResidueDesign;


typedef enum _Type_ResidueChainType{
  Type_ResidueChainType_Protein, 
  Type_ResidueChainType_Ligand, 
  Type_ResidueChainType_DNA, 
  Type_ResidueChainType_RNA, 
  Type_ResidueChainType_Water, 
  Type_ResidueChainType_Metal, 
} Type_ResidueChainType;

typedef enum _Type_ResidueIsTerminal{
  Type_ResidueIsNter,
  Type_ResidueIsCter,
  Type_ResidueIsNotTerminal
} Type_ResidueIsTerminal;

typedef enum _Type_ResidueDesignType{
  Type_ResidueDesignType_Fixed,
  Type_ResidueDesignType_Mutated,
  Type_ResidueDesignType_Rotameric,
  Type_ResidueDesignType_SmallMol,
  Type_ResidueDesignType_Catalytic,
}Type_ResidueDesignType;

typedef struct _Residue{
  StringArray patches;
  AtomArray atoms;
  BondSet bonds;
  char name[MAX_LENGTH_RESIDUE_NAME+1];
  char chainName[MAX_LENGTH_CHAIN_NAME+1];
  int posInChain;
  int nCbIn8A;
  Type_ResidueIsTerminal resiTerm;
  Type_ResidueDesignType designSiteType;
  double internalEnergy;
  double backboneEnergy;
  double phipsi[2];
  DoubleArray xtorsions;
} Residue;

int ResidueCreate(Residue* pThis);
int ResidueDestroy(Residue* pThis);
int ResidueCopy(Residue* pThis, Residue* pOther);
char* ResidueGetName(Residue* pThis);
int ResidueSetName(Residue* pThis, char* newName);
char* ResidueGetChainName(Residue* pThis);
int ResidueSetChainName(Residue* pThis, char* newChainName);
int ResidueGetPosInChain(Residue* pThis);
int ResidueSetPosInChain(Residue* pThis, int newPosInChain);
int ResidueSetDesignSiteFlag(Residue* pThis, Type_ResidueDesignType newFlag);
double ResidueGetCharge(Residue* pThis);
int ResidueGetPolarity(Residue* pThis, Type_ResiduePolarity* pPolarity);
int ResidueGetAtomCount(Residue* pThis);
Atom* ResidueGetAtom(Residue* pThis, int index);
Atom* ResidueGetAtomByName(Residue* pThis, char* atomName);
int ResidueFindAtom(Residue* pThis, char* atomName, int* pIndex);
int ResidueGetAtomXYZ(Residue* pThis, char* atomName, XYZ* pXYZ);
AtomArray* ResidueGetAllAtoms(Residue* pThis);
int ResidueInsertAtom(Residue* pThis, int newIndex, Atom* pNewAtom);
int ResidueAddAtom(Residue* pThis, Atom* pNewAtom);
int ResidueDeleteAtom(Residue* pThis, char* atomName);
int ResidueReadXYZFromPDB(Residue* pThis, FileReader* pPDBFileReader);
int ResidueAddAtomsFromAtomParams(Residue* pThis, AtomParamsSet* pAtomParams);


BondSet* ResidueGetBonds(Residue* pThis);
int ResidueAddBondsFromResiTopos(Residue* pThis, ResiTopoSet* pResiTopoCollection);
int ResidueShowInPDBFormat(Residue* pThis, char* header, char* chainName,int atomIndex, int resiIndex, BOOL showHydrogen, FILE* pFile);

int ResiduePatch(Residue* pThis, char* patchName,AtomParamsSet* pAtomParam,ResiTopoSet* pTopos);
int ResiduePatchCTER(Residue* pThis, char* patchName,AtomParamsSet* pAtomParam,ResiTopoSet* pTopos);
int ResiduePatchNTERorCTER(Residue* pThis, char* NTERorCTER,AtomParamsSet* pAtomParam,ResiTopoSet* pTopos);
StringArray* ResidueGetPatchingHistory(Residue* pThis);
int ResidueCalcAtomXYZ(Residue* pThis, ResiTopoSet* pResiTopos,Residue* pPrevResi, Residue* pNextResi,char* atomName, XYZ* pDestXYZ);
int ResidueCalcAllAtomXYZ(Residue* pThis, ResiTopoSet* pResiTopos,Residue* pPrevResi, Residue* pNextResi);
int ResiTopoCollectionGetIndex(ResiTopoSet* pThis, char* resiName, int *index);
int ResidueTopologyFindCharmmICIndex(ResidueTopology* pThis, char* atomDName, int *index);
int ResidueShowAtomParameter(Residue* pThis);
int ResidueShowBondInformation(Residue* pThis);
double ResidueAndResidueSidechainRMSD(Residue* pThis, Residue* pOther);
BOOL LigandResidueNameConflictWithAminoAcid(char* ligname);
int ResidueCheckAtomCoordinateValidity(Residue* pThis);




#endif //RESIDUE_H
