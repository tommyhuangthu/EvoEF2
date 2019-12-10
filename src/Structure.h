///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "Chain.h"
#include "Rotamer.h"
#include "DesignSite.h"
#include "SmallMol.h"


typedef struct _Structure{
  char name[MAX_LENGTH_STRUCTURE_NAME+1];
  Chain* chains;
  DesignSite *designSites;
  int chainNum;
  int designSiteCount;
} Structure;

int StructureCreate(Structure* pThis);
int StructureDestroy(Structure* pThis);
char* StructureGetName(Structure* pThis);
int StructureSetName(Structure* pThis, char* newName);
int StructureGetChainCount(Structure* pThis);
Chain* StructureGetChain(Structure* pThis, int index);
Chain* StructureFindChainByName(Structure* pThis, char* chainName);
int StructureFindChainIndex(Structure* pThis, char* chainName, int* index);
int StructureFindSmallMol(Structure* pThis, Residue** ppSmallMol);
int StructureAddChain(Structure* pThis, Chain* newChain);
int StructureDeleteChain(Structure* pThis, char* chainName);
int StructureShowInPDBFormat(Structure* pThis, BOOL showHydrogen, FILE* pFile);

//deal with design sites
int StructureGetDesignSiteCount(Structure* pThis);
DesignSite* StructureGetDesignSite(Structure* pThis, int index);
DesignSite* StructureFindDesignSite(Structure* pThis, int chainIndex, int resiIndex);
int StructureShowDesignSites(Structure* pThis, FILE* pFile);
int ProteinSiteAddDesignSite(Structure* pThis, int chainIndex, int resiIndex);
int StructureRemoveAllDesignSites(Structure* pThis);
int ProteinSiteRemoveDesignSite(Structure* pThis, int chainIndex, int resiIndex);

//other functions
int ChainComputeResiduePosition(Structure *pStructure, int chainIndex);
int StructureComputeResiduePosition(Structure *pStructure);
int StructureGetAminoAcidComposition(Structure* pStructure, int *aas);

//checker and debuggers
int StructureShowAtomParameter(Structure* pStructure);
int StructureShowBondInformation(Structure* pStructure);

int StructureConfig(Structure *pStructure, char* pdbFile, AtomParamsSet* pAtomParams, ResiTopoSet* pResiTopos);
int StructureConfigLigand(Structure* pStructure, char* mol2file, AtomParamsSet* pAtomParams, ResiTopoSet* pTopos);

DesignSite * StructureFindDesignSiteByChainName(Structure *pStructure, char *chainName, int posInChain);
int DesignSiteIndexFind(Structure *pStructure, char *chainName, int posInSeg);
int DesignSiteIndexGet(Structure *pStructure,int chainIndex, int resiIndex);
int StructureCalcProteinResidueSidechainTorsion(Structure* pThis, ResiTopoSet* pResiTopos);


#endif // STRUCTURE_H
