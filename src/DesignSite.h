///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#ifndef DESIGN_SITE_H
#define DESIGN_SITE_H

#include "Rotamer.h"

typedef struct _DesignSite{
  RotamerSet rotamers;
  Residue* pResidue;
  int chainIndex;
  int resiIndex;
} DesignSite;

int DesignSiteCreate(DesignSite* pThis);
int DesignSiteDestroy(DesignSite* pThis);
RotamerSet* DesignSiteGetRotamers(DesignSite* pThis);
char* DesignSiteGetChainName(DesignSite* pThis);
char* DesignSiteGetResiName(DesignSite* pThis);
int DesignSiteGetPosInChain(DesignSite* pThis);
int DesignSiteShow(DesignSite* pThis,FILE* pFile);
Residue* DesignSiteGetResidue(DesignSite* pThis);
int DesignSiteRemoveRotamers(DesignSite* pThis);


int DesignSiteShowRepresentativeRotamerAtomParameter(DesignSite* pThis);
int DesignSiteShowRepresentativeRotamerBondInformation(DesignSite* pThis);




#endif

