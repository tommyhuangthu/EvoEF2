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

#ifndef ENERGY_MATRIX_H
#define ENERGY_MATRIX_H


#include "Utility.h"
#include "Structure.h"
#include "ErrorTracker.h"
#include "EnergyFunction.h"

#define  SELF_ENERGY_SAME_THRESHOLD                                7.5 //20.0
#define  SELF_ENERGY_DIFFERENT_THRESHOLD                           15.0 //50.0
#define  PARTITION_SIZE_PARA                                       2000

#define  COMPUTATIONAL_AMOUNT_MAGNIFICATION                         1.0


typedef struct _EnergyMatrixBlock{
  int    DesignSiteI;
  int    DesignSiteK;
  int    RotamerCountSiteI;
  int    RotamerCountSiteK;
  double  *energyIK;
} EnergyMatrixBlock;

int     EnergyMatrixBlockCreate(EnergyMatrixBlock* pThis);
void    EnergyMatrixBlockDestroy(EnergyMatrixBlock* pThis);
int     EnergyMatrixBlockCopy(EnergyMatrixBlock* pThis,EnergyMatrixBlock* pOther);
double* EnergyMatrixBlockGet(EnergyMatrixBlock* pThis,int J,int S);
int     EnergyMatrixBlockGenerate(EnergyMatrixBlock* pThis,Structure *pStructure,
                    int designSiteI,int designSiteK,FILE* outputFile);

typedef struct _EnergyMatrix{
  int designSiteCount;
  EnergyMatrixBlock* blocks;
} EnergyMatrix;

typedef struct _Job{
  double computationAmount;
  int designSiteI;
  int designSiteK;
  int rotamerPartitionIndexOnSiteI;
  int rotamerPartitionIndexOnSiteK;
} Job;

typedef struct _SitePartition{
  int designSite;
  int partitionCount;
  int* rotamerCount;
  IntArray* rotamerIndexOld;
  IntArray* rotamerIndexNew;
} SitePartition;

typedef struct _Slot{
  Job* jobs;
  int jobCount;
  double totalComputationAmount;
} Slot;

int  EnergyMatrixCreate(EnergyMatrix* pThis);
void EnergyMatrixDestroy(EnergyMatrix* pThis);
int  EnergyMatrixCopy(EnergyMatrix* pThis,EnergyMatrix* pOther);
int  EnergyMatrixRead(EnergyMatrix* pThis,char* filepath);
int  EnergyMatrixWrite(EnergyMatrix* pThis,char* filepath);
int  EnergyMatrixCheck(char* filepath);
EnergyMatrixBlock* EnergyMatrixGetBlock(EnergyMatrix* pThis,int designSiteI,int designSiteK);
double* EnergyMatrixGet(EnergyMatrix* pThis,int designSiteI,int designSiteK,int rotamerIJ,int rotamerKS);
int  EnergyMatrixGetSiteCount(EnergyMatrix* pThis);
int  EnergyMatrixGetRotamerCount(EnergyMatrix* pThis,int designSiteI);
int EnergyMatrixShow(EnergyMatrix* pThis);
//'slotIndex' and 'slotCount' is designed for parallel computation,
//If not used in parallel computation environment, set both parameters to be 1
int EnergyMatrixGenerate(Structure* pStructure,char* energyMatrixFilePath,int slotIndex,int slotCount);



typedef struct _RotamerList{
  int designSiteCount;
  int* rotamerCount;
  int* remainRotamerCount;
  BOOL** remainFlag;
} RotamerList;

int RotamerListCreateFromStructure(RotamerList* pThis,Structure* pStructure);
int RotamerListCreateFromEnergyMatrix(RotamerList* pThis,EnergyMatrix* pEnergyMatrix);
void RotamerListDestroy(RotamerList* pThis);
int RotamerListCopy(RotamerList* pThis, RotamerList* pOther);
int RotamerListRead(RotamerList* pThis,char* filepath);
int RotamerListWrite(RotamerList* pThis,char* filepath);
int RotamerListShow(RotamerList* pThis);



int EnergyMatrixBlockUpdate(EnergyMatrixBlock* pEnergyBlock,BOOL* siteIRotamerDeletedFlag,BOOL* siteKRotamerDeletedFlag);
int RotamerListAndEnergyMatrixDelete(RotamerList* pRotamerList, EnergyMatrix* pEnergyMatrix, IntArray* pIndexOfDeletedRotamers);




int DesignSiteGetRemainRotamerCount(int designSiteI, RotamerList* pList);
int RotamerOriginalIndexGet(RotamerList* pList, int designSiteI, int rotamerIJ, int *trueIndexIJ);
int RotamerReducedIndexGet(RotamerList* pList, int designSiteI, int trueIndexIJ, int *reducedIndex);


//////////////////////////////////////////////////////////////
//Functions Used for Accelerating Energy Matrix Calculation
//////////////////////////////////////////////////////////////
int EnergyMatrixUpdateByRotamerList(EnergyMatrix* pOriginalMatrix, RotamerList* pList);
int EnergyMatrixItemRead(Structure* pStructure, RotamerList* pList, char* completeEnergyMatrixFile, double*** pppEnergyMatrixItem);
int EnergyMatrixItemDestroy(double*** pppEnergyMatrixItem, int designSiteCount);
int EnergyMatrixItemUpdateForWeight(Structure* pStructure, RotamerList* pList, double*** pppEnergyMatrixItem);
int EnergyMatrixCreateFromEnergyItemFile(EnergyMatrix* pThis, Structure* pStructure, RotamerList* pList, char* completeEnergyMatrixFile);




int EnergyMatrixGenerateBasedOnPartition(Structure* pStructure, RotamerList* pList, char* energyMatrixFilePath, int slotIndex, int slotCount);
int EnergyMatrixPartitionPairGenerate(Structure* pStructure, int designSiteI, int designSiteK, 
                    IntArray* pRotamerIndexOldOnSiteI, IntArray* pRotamerIndexOldOnSiteK, 
                    IntArray* pRotamerIndexNewOnSiteI, IntArray* pRotamerIndexNewOnSiteK, FILE* outputFile);
int EnergyMatrixPartitionPairGenerateStoreInMemory(Structure* pStructure, int designSiteI, int designSiteK, 
                           IntArray* pRotamerIndexOldOnSiteI, IntArray* pRotamerIndexOldOnSiteK, 
                           IntArray* pRotamerIndexNewOnSiteI, IntArray* pRotamerIndexNewOnSiteK, 
                           double* energyInMemory, int startIndexThisPartition);
int EnergyMatrixPartitionPairWriteEnergyToFile(Structure* pStructure, int designSiteI, int designSiteK, 
                         IntArray* pRotamerIndexOldOnSiteI, IntArray* pRotamerIndexOldOnSiteK, 
                         IntArray* pRotamerIndexNewOnSiteI, IntArray* pRotamerIndexNewOnSiteK, 
                         double* energyInMemory, int startIndexThisPartition,
                         FILE* outputFile);

int StructureShowDesignSitesAfterRotamerDelete(Structure* pStructure, RotamerList* pList,FILE* pFile);


int RotamerListUpdateByDeleteArray(RotamerList* pList, IntArray* pDeleteArray);

int RotamerDeleteBySelfEnergyCheckNew(EnergyMatrix *pMatrix, Structure *pStructure, RotamerList* pList, IntArray *pDeleteList, EnergyMatrix *pRemainFlag);

int SelfEnergyGenerate(Structure *pStructure, char* selfEnergyFilePath);
int SelfEnergyGenerateWithBBdepEnergy(Structure *pStructure,AAppTable* pAAppTable,RamaTable* pRamaTable,char* selfEnergyFilePath);
int SelfEnergyReadAndCheck(Structure* pStructure, RotamerList* pRotamerList,char* selfEnergyFile);

int EnergyMatrixReadNew(EnergyMatrix* pThis, RotamerList* pList, char* energyMatrixFile);
int DeleteArrayGenerateFromTwoRotamerList(IntArray* pDeleteRotamerArray, RotamerList* pOld, RotamerList* pNew);


int RotamerListCreateFromFile(RotamerList *pThis, char* filepath);


#endif //ENERGY_MATRIX_H

