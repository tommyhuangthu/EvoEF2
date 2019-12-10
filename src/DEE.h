///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

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