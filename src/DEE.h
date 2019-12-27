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

#ifndef DEE_H
#define DEE_H

#include "EnergyMatrix.h"
#define DEE_PARA_THRESHOLD  0.0

typedef struct _DEENodeUnifyRecord{
  int newSiteCount;
  int* oldSiteIndexOfNewSites;
  int whichTwoOldSitesHaveBeenUnified[2];
  int totalRotCountOnOldSites[2]; 

  int rotCountOfUnifiedNewSite;
  int* oldRotIndexFromSite1;
  int* oldRotIndexFromSite2;
} DEENodeUnifyRecord;

int DEE(EnergyMatrix* pEnergyMatrix,RotamerList* pRotamerList,double deeTreshold);
int DEERotamerListAndEnergyMatrixDelete(RotamerList* pRotamerList,
                                        EnergyMatrix* pEnergyMatrix,
                                        EnergyMatrix* pRemainFlag,
                                        IntArray* pIndexOfDeletedRotamers);
double DEECalcMinEnergy(EnergyMatrix* pEnergyMatrix,RotamerList* pList);
BOOL DEEGoldsteinCriteria(EnergyMatrix* pEnergyMatrix,int designSiteI,int rotamerJ,
                          int referenceRotamerOnSiteI,double deeThreshold,EnergyMatrix* pRemainFlag);
int DEEGoldstein(EnergyMatrix* pEnergyMatrix,IntArray* pDeleteList,double deeThreshold,EnergyMatrix* pRemainFlag);
BOOL DEESplitCriteria(EnergyMatrix* pEnergyMatrix,int designSiteI,int rotamerJ,double deeTreshold,EnergyMatrix* pRemainFlag);
int DEESplit(EnergyMatrix* pEnergyMatrix,IntArray* pDeleteList,double deeTreshold,EnergyMatrix* pRemainFlag);
int DEEDouble(EnergyMatrix* pEnergyMatrix,double deeThreshold,EnergyMatrix* pRemainFlag);
BOOL DEEDoubleCriteria(EnergyMatrix* pEnergyMatrix,int designSite[2],int rot[2],int refRot[2],double deeThreshold,EnergyMatrix* pRemainFlag);
int DEENodeUnifyRecordCreate(DEENodeUnifyRecord* pThis);
void DEENodeUnifyRecordDestroy(DEENodeUnifyRecord* pThis);
void DEENodeUnifyRecordShow(DEENodeUnifyRecord* pThis);
int DEENodeUnifyUpdateRotamerList(RotamerList* pList,DEENodeUnifyRecord* pRecord);
int DEENodeUnifyUpdateEnergyMatrix(EnergyMatrix* pEnergyMatrix,DEENodeUnifyRecord* pRecord);
int DEENodeUnifyUpdateRemainFlag(EnergyMatrix* pRemainFlag,DEENodeUnifyRecord* pRecord);
int DEENodeUnify(EnergyMatrix* pEnergyMatrix,EnergyMatrix* pRemainFlag,RotamerList* pList,DEENodeUnifyRecord* pRecord);
int DEENodeDeUnify(RotamerList* pList,DEENodeUnifyRecord* pRecord);

int DEEShowDesignedStructure(Structure* pStructure, RotamerList* pRotamerList, char* structureFile);


#endif //DEE