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
int DesignSiteCopy(DesignSite* pThis,DesignSite* pOther);



#endif

