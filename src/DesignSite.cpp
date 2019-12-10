///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#include "DesignSite.h"

int DesignSiteCreate(DesignSite* pThis){
  pThis->pResidue = NULL;
  RotamerSetCreate(&pThis->rotamers);
  pThis->chainIndex = -1;
  pThis->resiIndex = -1;
  return Success;
}
int DesignSiteDestroy(DesignSite* pThis){
  pThis->pResidue->designSiteType=Type_ResidueDesignType_Fixed;
  pThis->pResidue = NULL;
  RotamerSetDestroy(&pThis->rotamers);
  pThis->chainIndex = -1;
  pThis->resiIndex = -1;
  return Success;
}
RotamerSet* DesignSiteGetRotamers(DesignSite* pThis){
  return &pThis->rotamers;
}
int DesignSiteShow(DesignSite* pThis,FILE* pFile){
  printf("Chain %s residue %s %d has %d rotamers\n",ResidueGetChainName(pThis->pResidue),ResidueGetName(pThis->pResidue),ResidueGetPosInChain(pThis->pResidue),RotamerSetGetCount(&pThis->rotamers));
  return Success;
}


char* DesignSiteGetChainName(DesignSite* pThis){
  return pThis->pResidue->chainName;
}
char* DesignSiteGetResiName(DesignSite* pThis){
  return pThis->pResidue->name;
}
int DesignSiteGetPosInChain(DesignSite* pThis){
  return pThis->pResidue->posInChain;
}

Residue* DesignSiteGetResidue(DesignSite* pThis){
  return pThis->pResidue;
}

int DesignSiteShowRepresentativeRotamerAtomParameter(DesignSite* pThis){
  RotamerSet* pSet = DesignSiteGetRotamers(pThis);
  for(int i=0; i<RotamerSetGetRepresentativeCount(pSet); ++i){
    Rotamer* pRotamer=RotamerSetGetRepresentativeByIndex(pSet,i);
    printf("Rotamer %s Atom Parameter:\n", RotamerGetType(pRotamer));
    RotamerShowAtomParameter(pRotamer);
  }
  return Success;
}


int DesignSiteShowRepresentativeRotamerBondInformation(DesignSite* pThis){
  RotamerSet* pSet = DesignSiteGetRotamers(pThis);
  for(int i=0; i<RotamerSetGetRepresentativeCount(pSet); ++i){
    Rotamer* pRotamer=RotamerSetGetRepresentativeByIndex(pSet,i);
    printf("Rotamer %s Bond information:\n", RotamerGetType(pRotamer));
    RotamerShowBondInformation(pRotamer);
  }
  return Success;
}

int DesignSiteRemoveRotamers(DesignSite* pThis){
  RotamerSetDestroy(&pThis->rotamers);
  return RotamerSetCreate(&pThis->rotamers);
}


