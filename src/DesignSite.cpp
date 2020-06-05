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


int DesignSiteCopy(DesignSite* pThis,DesignSite* pOther){
  pThis->chainIndex=pOther->chainIndex;
  pThis->resiIndex=pOther->resiIndex;
  pThis->pResidue=pOther->pResidue;
  RotamerSetCopy(&pThis->rotamers,&pOther->rotamers);
  return Success;
}