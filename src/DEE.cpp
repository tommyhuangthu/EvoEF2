///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#include "DEE.h"
#include <string.h>
#include "SmallMol.h"



int DEERotamerListAndEnergyMatrixDelete(RotamerList* pRotamerList,
                    EnergyMatrix* pEnergyMatrix,
                    EnergyMatrix* pRemainFlag,
                    IntArray* pIndexOfDeletedRotamers)
{
  int i,j;
  int result;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  BOOL** rotamerDeletedFlag = (BOOL**)malloc(sizeof(BOOL*)*pRotamerList->designSiteCount);

  for(i=0;i<pRotamerList->designSiteCount;i++){
    rotamerDeletedFlag[i] = (BOOL*)malloc(sizeof(BOOL)*EnergyMatrixGetBlock(pEnergyMatrix,i,i)->RotamerCountSiteI);
    for(j=0;j<EnergyMatrixGetBlock(pEnergyMatrix,i,i)->RotamerCountSiteI;j++){
      rotamerDeletedFlag[i][j] = FALSE;
    }
  }

  if(IntArrayGetLength(pIndexOfDeletedRotamers)%2!=0){
    result = AssertionError;
    sprintf(errMsg,"in file %s function %s line %d, parameter #3 must be a int array with even number length",__FILE__,__FUNCTION__,__LINE__);
    TraceError(errMsg,result);
    return result;
  }
  for(i=0;i<IntArrayGetLength(pIndexOfDeletedRotamers);i+=2){
    int siteIndex = IntArrayGet(pIndexOfDeletedRotamers,i);
    int rotamerIndex = IntArrayGet(pIndexOfDeletedRotamers,i+1);
    if( siteIndex >= pRotamerList->designSiteCount ||
      rotamerIndex >= EnergyMatrixGetBlock(pEnergyMatrix,siteIndex,siteIndex)->RotamerCountSiteI ){
        result = AssertionError;
        sprintf(errMsg,"in file %s function %s line %d, rotamer #%d on site #%d does not exist",__FILE__,__FUNCTION__,__LINE__,rotamerIndex,siteIndex);
        TraceError(errMsg,result);
        return AssertionError;
    }
    rotamerDeletedFlag[siteIndex][rotamerIndex] = TRUE;
  }

  //Update the RemainRotamerList
  for(i=0;i<pRotamerList->designSiteCount;i++){
    int remainCount = 0;
    for(j=0;j<pRotamerList->rotamerCount[i];j++){
      if(pRotamerList->remainFlag[i][j] == TRUE){
        if(rotamerDeletedFlag[i][remainCount]==TRUE){
          pRotamerList->remainFlag[i][j] = FALSE;
        }
        remainCount++;      
      }
    }
    if(remainCount==0){
      result = AssertionError;
      sprintf(errMsg,"in file %s function %s line %d, all rotamers on design site %d have been deleted",__FILE__,__FUNCTION__,__LINE__,i);
      TraceError(errMsg,result);
      return result;
    }
  }

  //Update the EnergyMatrix and RemainFlag
  for(i=0;i<pRotamerList->designSiteCount;i++){
    for(j=i;j<pRotamerList->designSiteCount;j++){
      EnergyMatrixBlockUpdate(EnergyMatrixGetBlock(pEnergyMatrix,i,j),
        rotamerDeletedFlag[i],rotamerDeletedFlag[j]);
      EnergyMatrixBlockUpdate(EnergyMatrixGetBlock(pRemainFlag,i,j),
        rotamerDeletedFlag[i],rotamerDeletedFlag[j]);
    }
  }

  

  for(i=0;i<pRotamerList->designSiteCount;i++){
    free(rotamerDeletedFlag[i]);
  }
  free(rotamerDeletedFlag);
  return Success;
}


BOOL DEEGoldsteinCriteria(EnergyMatrix* pEnergyMatrix,int designSiteI,int rotamerJ,
             int referenceRotamerOnSiteI,double deeThreshold,EnergyMatrix* pRemainFlag)
{
  //Test whether RotamerJ on DesignSiteI can be eliminated based on referenceRotamer
  int otherDesignSite;
  double energyDiff;

  if(referenceRotamerOnSiteI == rotamerJ){
    return FALSE;
  }

  energyDiff = *EnergyMatrixGet(pEnergyMatrix,designSiteI,designSiteI,rotamerJ,rotamerJ)
         - *EnergyMatrixGet(pEnergyMatrix,designSiteI,designSiteI,referenceRotamerOnSiteI,referenceRotamerOnSiteI);

  for(otherDesignSite=0; otherDesignSite<pEnergyMatrix->designSiteCount; otherDesignSite++){
    int rotamerOnOtherSite;
    double minDiff = 1e8;

    if(otherDesignSite==designSiteI){
      continue;
    }

    for(rotamerOnOtherSite=0; 
      rotamerOnOtherSite<EnergyMatrixGetRotamerCount(pEnergyMatrix,otherDesignSite); 
      rotamerOnOtherSite++){
        double diff;

        if(*EnergyMatrixGet(pRemainFlag,otherDesignSite,otherDesignSite,rotamerOnOtherSite,rotamerOnOtherSite)<0.0){
          continue;
        }
        if(*EnergyMatrixGet(pRemainFlag,designSiteI,otherDesignSite,rotamerJ,rotamerOnOtherSite)<0.0){
          continue;
        }

        diff = *EnergyMatrixGet(pEnergyMatrix,designSiteI,otherDesignSite,rotamerJ,rotamerOnOtherSite)
          - *EnergyMatrixGet(pEnergyMatrix,designSiteI,otherDesignSite,referenceRotamerOnSiteI,rotamerOnOtherSite);
        if(diff < minDiff){
          minDiff = diff;
      }
    }

    energyDiff += minDiff;

  }


  if(energyDiff > deeThreshold){
    return TRUE;
  }
  else{
    return FALSE;
  }

}

int DEEGoldstein(EnergyMatrix* pEnergyMatrix,IntArray* pDeleteList,double deeThreshold,EnergyMatrix* pRemainFlag){
  int designSiteI;
  IntArrayResize(pDeleteList,0);

  for(designSiteI=0;designSiteI<pEnergyMatrix->designSiteCount;designSiteI++){
    int rotamerJ;
    int referenceRotamer;
    int rotamerCountOnSiteI = EnergyMatrixGetBlock(pEnergyMatrix,designSiteI,designSiteI)->RotamerCountSiteI;
    for(rotamerJ=0; rotamerJ<rotamerCountOnSiteI; rotamerJ++){
      for(referenceRotamer=0; referenceRotamer<rotamerCountOnSiteI; referenceRotamer++){
        if(referenceRotamer==rotamerJ){
          continue;
        }
        if(*EnergyMatrixGet(pRemainFlag,designSiteI,designSiteI,referenceRotamer,referenceRotamer)<0.0){
          continue;
        }
        if(DEEGoldsteinCriteria(pEnergyMatrix,designSiteI,rotamerJ,referenceRotamer,deeThreshold,pRemainFlag)){
          IntArrayAppend(pDeleteList,designSiteI);
          IntArrayAppend(pDeleteList,rotamerJ);
          *EnergyMatrixGet(pRemainFlag,designSiteI,designSiteI,rotamerJ,rotamerJ) = -1.0;
          break;
        }
      }
    }
  }

  return Success;
}


BOOL DEESplitCriteria(EnergyMatrix* pEnergyMatrix,int designSiteI,int rotamerJ,double deeTreshold,EnergyMatrix* pRemainFlag){
  int i;
  int refRot;
  int splitSite;
  int canBeEliminated;
  double** minEnergyDiff      = (double**)malloc(sizeof(double*)*EnergyMatrixGetRotamerCount(pEnergyMatrix,designSiteI));
  double*  sumOfminEnergyDiff = (double*)malloc(sizeof(double)* EnergyMatrixGetRotamerCount(pEnergyMatrix,designSiteI));
  

  
  for(i=0;i<EnergyMatrixGetRotamerCount(pEnergyMatrix,designSiteI);i++){
    minEnergyDiff[i] = (double*)malloc(sizeof(double)*EnergyMatrixGetSiteCount(pEnergyMatrix));
  }


  for(refRot=0; refRot<EnergyMatrixGetRotamerCount(pEnergyMatrix,designSiteI);refRot++){
    int designSiteK;
    sumOfminEnergyDiff[refRot] = 0.0;
    
    for(designSiteK=0;designSiteK<EnergyMatrixGetSiteCount(pEnergyMatrix);designSiteK++){
      int rotamerKS;
      if(refRot==rotamerJ || designSiteK==designSiteI){
        minEnergyDiff[refRot][designSiteK] = 0.0;
        continue;
      }
      else{
        minEnergyDiff[refRot][designSiteK] = 1e10;
      }

      for(rotamerKS=0; rotamerKS<EnergyMatrixGetRotamerCount(pEnergyMatrix,designSiteK); rotamerKS++){
        double diff;
        if(*EnergyMatrixGet(pRemainFlag,designSiteK,designSiteK,rotamerKS,rotamerKS)<0.0){
          continue;
        }
        if(*EnergyMatrixGet(pRemainFlag,designSiteI,designSiteK,rotamerJ,rotamerKS)<0.0){
          continue;
        }
        diff = *EnergyMatrixGet(pEnergyMatrix,designSiteI,designSiteK,rotamerJ,rotamerKS)
                - *EnergyMatrixGet(pEnergyMatrix,designSiteI,designSiteK,refRot,  rotamerKS);
        if(diff < minEnergyDiff[refRot][designSiteK]){
          minEnergyDiff[refRot][designSiteK] = diff;
        }
      }
      sumOfminEnergyDiff[refRot] += minEnergyDiff[refRot][designSiteK]; 
    }
  }

  canBeEliminated = FALSE;
  for(splitSite=0; splitSite<EnergyMatrixGetSiteCount(pEnergyMatrix);splitSite++){
    int rotOnSplitSite;
    BOOL splitSiteSuccessful = TRUE;
    if(splitSite==designSiteI){
      continue;
    }
    for(rotOnSplitSite=0; rotOnSplitSite<EnergyMatrixGetRotamerCount(pEnergyMatrix,splitSite); rotOnSplitSite++){
      BOOL foundInTheSplittedConformation = FALSE;
      for(refRot=0; refRot<EnergyMatrixGetRotamerCount(pEnergyMatrix,designSiteI); refRot++ ){
        double diff;
        if(refRot==rotamerJ){
          continue;
        }
        diff = *EnergyMatrixGet(pEnergyMatrix,designSiteI,designSiteI,rotamerJ,rotamerJ)
           - *EnergyMatrixGet(pEnergyMatrix,designSiteI,designSiteI,refRot,  refRot)
           + *EnergyMatrixGet(pEnergyMatrix,designSiteI,splitSite,  rotamerJ,rotOnSplitSite)
           - *EnergyMatrixGet(pEnergyMatrix,designSiteI,splitSite,  refRot,  rotOnSplitSite)
           + sumOfminEnergyDiff[refRot]
           - minEnergyDiff[refRot][splitSite];
        if(diff > deeTreshold){
          foundInTheSplittedConformation = TRUE;
          break;
        }
      }
      if(foundInTheSplittedConformation == FALSE){
        splitSiteSuccessful = FALSE;
        break;
      }
    }
    if(splitSiteSuccessful == FALSE){
      continue;
    }
    else{
      canBeEliminated = TRUE;
      break;
    }
    
  }




  for(i=0;i<EnergyMatrixGetRotamerCount(pEnergyMatrix,designSiteI);i++){
    free(minEnergyDiff[i]);
  }
  free(minEnergyDiff);
  free(sumOfminEnergyDiff);

  return canBeEliminated;
}

int DEEDouble(EnergyMatrix* pEnergyMatrix,double deeThreshold,EnergyMatrix* pRemainFlag){
  int eliminatedPairCount = 0;
  int site[2];
  int criteriaRunnedTimes = 0;

  for(site[0]=0; site[0]<EnergyMatrixGetSiteCount(pEnergyMatrix); site[0]++){
    for(site[1]=site[0]+1; site[1]<EnergyMatrixGetSiteCount(pEnergyMatrix); site[1]++){
      int rot[2];
      for(rot[0]=0; rot[0]<EnergyMatrixGetRotamerCount(pEnergyMatrix,site[0]); rot[0]++){
        for(rot[1]=0; rot[1]<EnergyMatrixGetRotamerCount(pEnergyMatrix,site[1]); rot[1]++){
          int refRot[2];
          BOOL rotPairDeleted = FALSE;
          for(refRot[0]=0; refRot[0]<EnergyMatrixGetRotamerCount(pEnergyMatrix,site[0]); refRot[0]++){
            for(refRot[1]=0; refRot[1]<EnergyMatrixGetRotamerCount(pEnergyMatrix,site[1]); refRot[1]++){
              int result;
              if(rotPairDeleted){
                break;
              }
              if(refRot[0]==rot[0] && refRot[1]==rot[1]){
                continue;
              }
              result = DEEDoubleCriteria(pEnergyMatrix,site,rot,refRot,deeThreshold,pRemainFlag);
              if(result){
                *EnergyMatrixGet(pRemainFlag,site[0],site[1],rot[0],rot[1]) = -1.0;
                eliminatedPairCount++;
                rotPairDeleted = TRUE;
              }
              criteriaRunnedTimes++;

            }
          }
        }
      }
    }
  }
  return eliminatedPairCount;
}
BOOL DEEDoubleCriteria(EnergyMatrix* pEnergyMatrix,int designSite[2],int rot[2],int refRot[2],
             double deeThreshold,EnergyMatrix* pRemainFlag)
{
  int otherSite;
  int rotOnOtherSite;
  double energyDiff;
  if(*EnergyMatrixGet(pRemainFlag,designSite[0],designSite[1],rot[0],rot[1])<0.0){
    return FALSE;
  }
  if(*EnergyMatrixGet(pRemainFlag,designSite[0],designSite[1],refRot[0],refRot[1])<0.0){
    return FALSE;
  }
  


  energyDiff = *EnergyMatrixGet(pEnergyMatrix,designSite[0],designSite[0],rot[0],rot[0])
         + *EnergyMatrixGet(pEnergyMatrix,designSite[1],designSite[1],rot[1],rot[1])
         + *EnergyMatrixGet(pEnergyMatrix,designSite[0],designSite[1],rot[0],rot[1])
         - *EnergyMatrixGet(pEnergyMatrix,designSite[0],designSite[0],refRot[0],refRot[0])
         - *EnergyMatrixGet(pEnergyMatrix,designSite[1],designSite[1],refRot[1],refRot[1])
         - *EnergyMatrixGet(pEnergyMatrix,designSite[0],designSite[1],refRot[0],refRot[1]);

  for(otherSite=0; otherSite<EnergyMatrixGetSiteCount(pEnergyMatrix); otherSite++){
    double minDiff = 1e8;
    int site0Increment;
    int site1Increment;
    int rotCountOnOtherSite = EnergyMatrixGetRotamerCount(pEnergyMatrix,otherSite);
    double* pRemainFlag_site0_otherSite_rot0_rotOnOtherSite = EnergyMatrixGet(pRemainFlag,designSite[0],otherSite,rot[0],0);
    double* pRemainFlag_site1_otherSite_rot1_rotOnOtherSite = EnergyMatrixGet(pRemainFlag,designSite[1],otherSite,rot[1],0);
    double* pEM_Site0_otherSite_rot0_rotOnOtherSite = EnergyMatrixGet(pEnergyMatrix,designSite[0],otherSite,rot[0],0);
    double* pEM_Site1_otherSite_rot1_rotOnOtherSite = EnergyMatrixGet(pEnergyMatrix,designSite[1],otherSite,rot[1],0);
    double* pEM_Site0_otherSite_refRot0_rotOnOtherSite = EnergyMatrixGet(pEnergyMatrix,designSite[0],otherSite,refRot[0],0);
    double* pEM_Site1_otherSite_refRot1_rotOnOtherSite = EnergyMatrixGet(pEnergyMatrix,designSite[1],otherSite,refRot[1],0);
    if(otherSite > designSite[0]){
      site0Increment = 1;
    }
    else{
      site0Increment = EnergyMatrixGetRotamerCount(pEnergyMatrix,designSite[0]);
    }
    if(otherSite > designSite[1]){
      site1Increment = 1;
    }
    else{
      site1Increment = EnergyMatrixGetRotamerCount(pEnergyMatrix,designSite[1]);
    }


    if(otherSite==designSite[0] || otherSite==designSite[1]){
      continue;
    }
    for(rotOnOtherSite=0; rotOnOtherSite<rotCountOnOtherSite; rotOnOtherSite++){
      double diff;
      BOOL skip;
      
      if( *pRemainFlag_site0_otherSite_rot0_rotOnOtherSite<0.0 ||
        *pRemainFlag_site1_otherSite_rot1_rotOnOtherSite<0.0){
          skip = TRUE;
      }
      else{
        skip = FALSE;
      }

      if(!skip){
        diff = *pEM_Site0_otherSite_rot0_rotOnOtherSite
           + *pEM_Site1_otherSite_rot1_rotOnOtherSite
           - *pEM_Site0_otherSite_refRot0_rotOnOtherSite
           - *pEM_Site1_otherSite_refRot1_rotOnOtherSite;
        if(diff<minDiff){
          minDiff = diff;
        }
      }

      pRemainFlag_site0_otherSite_rot0_rotOnOtherSite += site0Increment;
      pRemainFlag_site1_otherSite_rot1_rotOnOtherSite += site1Increment;
      pEM_Site0_otherSite_rot0_rotOnOtherSite += site0Increment;
      pEM_Site0_otherSite_refRot0_rotOnOtherSite += site0Increment;
      pEM_Site1_otherSite_rot1_rotOnOtherSite += site1Increment;
      pEM_Site1_otherSite_refRot1_rotOnOtherSite += site1Increment;

    }
    energyDiff += minDiff;
  }
  //printf("Site %d %d, Rot %d %d, Ref %d %d, Total : %f\n",
  //  designSite[0],designSite[1],
  //  rot[0],rot[1],refRot[0],refRot[1],energyDiff);
  
  if(energyDiff > deeThreshold){
    return TRUE;
  }
  else{
    return FALSE;
  }

}


int DEESplit(EnergyMatrix* pEnergyMatrix,IntArray* pDeleteList,double deeTreshold,EnergyMatrix* pRemainFlag){
  int designSiteI;
  int rotamerIJ;

  IntArrayResize(pDeleteList,0);

  for(designSiteI=0;designSiteI<EnergyMatrixGetSiteCount(pEnergyMatrix);designSiteI++){
    int rotamerCountOnSiteI = EnergyMatrixGetRotamerCount(pEnergyMatrix,designSiteI);
    for(rotamerIJ=0;rotamerIJ< rotamerCountOnSiteI; rotamerIJ++){
      if(DEESplitCriteria(pEnergyMatrix,designSiteI,rotamerIJ,deeTreshold,pRemainFlag)){
        IntArrayAppend(pDeleteList,designSiteI);
        IntArrayAppend(pDeleteList,rotamerIJ);
      }
    }
  }
  
  return Success;
}
double DEECalcMinEnergy(EnergyMatrix* pEnergyMatrix,RotamerList* pList){
  int i;
  double totalEnergy = 0.0;

  for(i=0;i<EnergyMatrixGetSiteCount(pEnergyMatrix);i++){
    int k;
    totalEnergy += *EnergyMatrixGet(pEnergyMatrix,i,i,0,0);
    for(k=i+1;k<EnergyMatrixGetSiteCount(pEnergyMatrix);k++){
      totalEnergy += *EnergyMatrixGet(pEnergyMatrix,i,k,0,0);
    }
  }
  return totalEnergy;
}



int DEENodeUnifyRecordCreate(DEENodeUnifyRecord* pThis){
  pThis->oldSiteIndexOfNewSites = NULL;
  pThis->oldRotIndexFromSite1 = NULL;
  pThis->oldRotIndexFromSite2 = NULL;
  return Success;
}
void DEENodeUnifyRecordDestroy(DEENodeUnifyRecord* pThis){
  free(pThis->oldSiteIndexOfNewSites);
  free(pThis->oldRotIndexFromSite1);
  free(pThis->oldRotIndexFromSite2);
}
void DEENodeUnifyRecordShow(DEENodeUnifyRecord* pThis){
  int i;
  printf("NewSiteCount: %d\n",pThis->newSiteCount);
  printf("Site %d and %d unified to Site %d\n",
    pThis->whichTwoOldSitesHaveBeenUnified[0],pThis->whichTwoOldSitesHaveBeenUnified[1],
    pThis->newSiteCount-1);
  printf("There are originally %d and %d rotamers on the old sites\n",
    pThis->totalRotCountOnOldSites[0],
    pThis->totalRotCountOnOldSites[1]);
  for(i=0;i<pThis->newSiteCount;i++){
    printf("Site %d : %d\n",i,pThis->oldSiteIndexOfNewSites[i]);
  }
  printf("New Site have %d rotamers:\n",pThis->rotCountOfUnifiedNewSite);
  for(i=0;i<pThis->rotCountOfUnifiedNewSite;i++){
    printf("Rot %d comes from %d and %d\n",i,pThis->oldRotIndexFromSite1[i],pThis->oldRotIndexFromSite2[i]);
  }
}


int DEENodeUnifyUpdateRotamerList(RotamerList* pList,DEENodeUnifyRecord* pRecord){
  int i,j,counter;
  BOOL** newRemainFlag;
  int* newRotamerCount;
  newRotamerCount = (int*)malloc(sizeof(int)*pRecord->newSiteCount);
  newRemainFlag = (BOOL**)malloc(sizeof(BOOL*)*pRecord->newSiteCount);
  counter = 0;
  for(i=0;i<pList->designSiteCount;i++){
    if(i==pRecord->whichTwoOldSitesHaveBeenUnified[0] || i==pRecord->whichTwoOldSitesHaveBeenUnified[1]){
      continue;
    }
    newRotamerCount[counter] = pList->rotamerCount[i];
    newRemainFlag[counter] = (BOOL*)malloc(sizeof(BOOL)*pList->rotamerCount[i]);
    for(j=0;j<newRotamerCount[counter];j++){
      newRemainFlag[counter][j] = pList->remainFlag[i][j];
    }
    counter++;
  }

  newRotamerCount[pRecord->newSiteCount-1] = pRecord->rotCountOfUnifiedNewSite;
  newRemainFlag[pRecord->newSiteCount-1] = (BOOL*)malloc(sizeof(BOOL)*pRecord->rotCountOfUnifiedNewSite);
  for(i=0;i<pRecord->rotCountOfUnifiedNewSite;i++){
    newRemainFlag[pRecord->newSiteCount-1][i] = TRUE;
  }
  
  RotamerListDestroy(pList);
  pList->designSiteCount = pRecord->newSiteCount;
  pList->remainFlag = newRemainFlag;
  pList->rotamerCount = newRotamerCount;
  return Success;
}
int DEENodeUnifyUpdateEnergyMatrix(EnergyMatrix* pEnergyMatrix,DEENodeUnifyRecord* pRecord){
  int i,k;
  EnergyMatrix newMatrix;
  EnergyMatrixCreate(&newMatrix);
  
  
  newMatrix.designSiteCount = pEnergyMatrix->designSiteCount-1;
  newMatrix.blocks = (EnergyMatrixBlock*)malloc(
    sizeof(EnergyMatrixBlock)*newMatrix.designSiteCount*newMatrix.designSiteCount);
  
  for(i=0;i<newMatrix.designSiteCount-1;i++){
    for(k=i;k<newMatrix.designSiteCount-1;k++){
      int oldSite1 = pRecord->oldSiteIndexOfNewSites[i];
      int oldSIte2 = pRecord->oldSiteIndexOfNewSites[k];
      EnergyMatrixBlockCreate(EnergyMatrixGetBlock(&newMatrix,i,k));
      EnergyMatrixBlockCopy(EnergyMatrixGetBlock(&newMatrix,i,k),
        EnergyMatrixGetBlock(pEnergyMatrix,oldSite1,oldSIte2));
    }
  }

  //The last column of the new energy matrix needs to be calculated
  for(i=0;i<newMatrix.designSiteCount;i++){
    EnergyMatrixBlock* pBlock;
    int oldSite1 = pRecord->whichTwoOldSitesHaveBeenUnified[0];
    int oldSite2 = pRecord->whichTwoOldSitesHaveBeenUnified[1];
    int j;
    int s;
    
    pBlock = EnergyMatrixGetBlock(&newMatrix,i,newMatrix.designSiteCount-1);
    EnergyMatrixBlockCreate(pBlock);
    pBlock->DesignSiteI = i;
    pBlock->DesignSiteK = newMatrix.designSiteCount-1;

    if(i==newMatrix.designSiteCount-1){
      pBlock->RotamerCountSiteI = pBlock->RotamerCountSiteK = pRecord->rotCountOfUnifiedNewSite;
      pBlock->energyIK = (double*)malloc(sizeof(double)*pBlock->RotamerCountSiteI*pBlock->RotamerCountSiteK);
      for(j=0;j<pBlock->RotamerCountSiteI;j++){
        int rotIndexOnSite1 = j/EnergyMatrixGetRotamerCount(pEnergyMatrix,oldSite2);
        int rotIndexOnSite2 = j%EnergyMatrixGetRotamerCount(pEnergyMatrix,oldSite2);
        *EnergyMatrixBlockGet(pBlock,j,j) = 
          *EnergyMatrixGet(pEnergyMatrix,oldSite1,oldSite1,rotIndexOnSite1,rotIndexOnSite1)
          + *EnergyMatrixGet(pEnergyMatrix,oldSite2,oldSite2,rotIndexOnSite2,rotIndexOnSite2)
          + *EnergyMatrixGet(pEnergyMatrix,oldSite1,oldSite2,rotIndexOnSite1,rotIndexOnSite2);
      }
    }
    else{
      pBlock->RotamerCountSiteI = EnergyMatrixGetRotamerCount(pEnergyMatrix,pRecord->oldSiteIndexOfNewSites[i]);
      pBlock->RotamerCountSiteK = pRecord->rotCountOfUnifiedNewSite;
      pBlock->energyIK = (double*)malloc(sizeof(double)*pBlock->RotamerCountSiteI*pBlock->RotamerCountSiteK);
      for(j=0;j<pBlock->RotamerCountSiteI;j++){
        for(s=0;s<pBlock->RotamerCountSiteK;s++){
          int rotIndexOnSite1 = s/EnergyMatrixGetRotamerCount(pEnergyMatrix,oldSite2);
          int rotIndexOnSite2 = s%EnergyMatrixGetRotamerCount(pEnergyMatrix,oldSite2);
          *EnergyMatrixBlockGet(pBlock,j,s) = 
            *EnergyMatrixGet(pEnergyMatrix,pRecord->oldSiteIndexOfNewSites[i],oldSite1,j,rotIndexOnSite1)
          +   *EnergyMatrixGet(pEnergyMatrix,pRecord->oldSiteIndexOfNewSites[i],oldSite2,j,rotIndexOnSite2);
          
          
        }
      }
    }
  }
  




  EnergyMatrixDestroy(pEnergyMatrix);
  pEnergyMatrix->designSiteCount = newMatrix.designSiteCount;
  pEnergyMatrix->blocks = newMatrix.blocks;
  return Success;
}

int DEENodeUnifyUpdateRemainFlag(EnergyMatrix* pRemainFlag,DEENodeUnifyRecord* pRecord){
  int i,k;
  EnergyMatrix newMatrix;
  EnergyMatrixCreate(&newMatrix);


  newMatrix.designSiteCount = pRemainFlag->designSiteCount-1;
  newMatrix.blocks = (EnergyMatrixBlock*)malloc(
    sizeof(EnergyMatrixBlock)*newMatrix.designSiteCount*newMatrix.designSiteCount);

  for(i=0;i<newMatrix.designSiteCount-1;i++){
    for(k=i;k<newMatrix.designSiteCount-1;k++){
      int oldSite1 = pRecord->oldSiteIndexOfNewSites[i];
      int oldSIte2 = pRecord->oldSiteIndexOfNewSites[k];
      EnergyMatrixBlockCreate(EnergyMatrixGetBlock(&newMatrix,i,k));
      EnergyMatrixBlockCopy(EnergyMatrixGetBlock(&newMatrix,i,k),
        EnergyMatrixGetBlock(pRemainFlag,oldSite1,oldSIte2));
    }
  }

  //The last column of the new energy matrix needs to be calculated
  for(i=0;i<newMatrix.designSiteCount;i++){
    EnergyMatrixBlock* pBlock;
    int oldSite1 = pRecord->whichTwoOldSitesHaveBeenUnified[0];
    int oldSite2 = pRecord->whichTwoOldSitesHaveBeenUnified[1];
    int j;
    int s;

    pBlock = EnergyMatrixGetBlock(&newMatrix,i,newMatrix.designSiteCount-1);
    EnergyMatrixBlockCreate(pBlock);
    pBlock->DesignSiteI = i;
    pBlock->DesignSiteK = newMatrix.designSiteCount-1;

    if(i==newMatrix.designSiteCount-1){
      pBlock->RotamerCountSiteI = pBlock->RotamerCountSiteK = pRecord->rotCountOfUnifiedNewSite;
      pBlock->energyIK = (double*)malloc(sizeof(double)*pBlock->RotamerCountSiteI*pBlock->RotamerCountSiteK);
      for(j=0;j<pBlock->RotamerCountSiteI;j++){
        int rotIndexOnSite1 = j/EnergyMatrixGetRotamerCount(pRemainFlag,oldSite2);
        int rotIndexOnSite2 = j%EnergyMatrixGetRotamerCount(pRemainFlag,oldSite2);
        if( *EnergyMatrixGet(pRemainFlag,oldSite1,oldSite1,rotIndexOnSite1,rotIndexOnSite1)<0.0 ||
          *EnergyMatrixGet(pRemainFlag,oldSite2,oldSite2,rotIndexOnSite2,rotIndexOnSite2)<0.0 ||
          *EnergyMatrixGet(pRemainFlag,oldSite1,oldSite2,rotIndexOnSite1,rotIndexOnSite2)<0.0 ){
            *EnergyMatrixBlockGet(pBlock,j,j) = -1.0;
        }
        else{
          *EnergyMatrixBlockGet(pBlock,j,j) = 1.0;
        }
      }
    }
    else{
      pBlock->RotamerCountSiteI = EnergyMatrixGetRotamerCount(pRemainFlag,pRecord->oldSiteIndexOfNewSites[i]);
      pBlock->RotamerCountSiteK = pRecord->rotCountOfUnifiedNewSite;
      pBlock->energyIK = (double*)malloc(sizeof(double)*pBlock->RotamerCountSiteI*pBlock->RotamerCountSiteK);
      for(j=0;j<pBlock->RotamerCountSiteI;j++){
        for(s=0;s<pBlock->RotamerCountSiteK;s++){
          int rotIndexOnSite1 = s/EnergyMatrixGetRotamerCount(pRemainFlag,oldSite2);
          int rotIndexOnSite2 = s%EnergyMatrixGetRotamerCount(pRemainFlag,oldSite2);
          if( *EnergyMatrixGet(pRemainFlag,pRecord->oldSiteIndexOfNewSites[i],oldSite1,j,rotIndexOnSite1)<0.0 ||
            *EnergyMatrixGet(pRemainFlag,pRecord->oldSiteIndexOfNewSites[i],oldSite2,j,rotIndexOnSite2)<0.0 ){
              *EnergyMatrixBlockGet(pBlock,j,s) = -1.0;
          }else{
            *EnergyMatrixBlockGet(pBlock,j,s) = 1.0;
          }
        }
      }
    }
  }

  EnergyMatrixDestroy(pRemainFlag);
  pRemainFlag->designSiteCount = newMatrix.designSiteCount;
  pRemainFlag->blocks = newMatrix.blocks;
  return Success;
}

int DEENodeUnify(EnergyMatrix* pEnergyMatrix,EnergyMatrix* pRemainFlag,RotamerList* pList,DEENodeUnifyRecord* pRecord){
  int i,j,counter;
  int unifySite[2] = {-1,-1};
  int unifySiteRotCount[2] = {1000000,1000000};
  //Find the two sites with fewest rotamers
  if(pList->designSiteCount<2){
    TraceError("In DEENodeUnify(), only 1 site remained, do not need unifying",AssertionError);
    return AssertionError;
  }
  for(i=0;i<EnergyMatrixGetSiteCount(pEnergyMatrix);i++){
    int rotCount = EnergyMatrixGetRotamerCount(pEnergyMatrix,i);
    if(rotCount<unifySiteRotCount[0]){
      unifySiteRotCount[1] = unifySiteRotCount[0];
      unifySiteRotCount[0] = rotCount;
      unifySite[1] = unifySite[0];
      unifySite[0] = i;
    }
    else if(rotCount>=unifySiteRotCount[0] && rotCount<unifySiteRotCount[1]){
      unifySiteRotCount[1] = rotCount;
      unifySite[1] = i;
    }
  }
  if(unifySite[0]>unifySite[1]){
    int swap = unifySite[1];
    unifySite[1] = unifySite[0];
    unifySite[0] = swap;      
  }
  
  //Recording
  DEENodeUnifyRecordDestroy(pRecord);
  DEENodeUnifyRecordCreate(pRecord);
  pRecord->newSiteCount = EnergyMatrixGetSiteCount(pEnergyMatrix)-1;
  pRecord->oldSiteIndexOfNewSites = (int*)malloc(sizeof(int)*pRecord->newSiteCount);
  pRecord->whichTwoOldSitesHaveBeenUnified[0] = unifySite[0];
  pRecord->whichTwoOldSitesHaveBeenUnified[1] = unifySite[1];
  pRecord->totalRotCountOnOldSites[0] = pList->rotamerCount[unifySite[0]];
  pRecord->totalRotCountOnOldSites[1] = pList->rotamerCount[unifySite[1]];

  //The new unified site is placed as the last site
  counter = 0;
  for(i=0;i<EnergyMatrixGetSiteCount(pEnergyMatrix);i++){
    if(i==unifySite[0] || i==unifySite[1]){
      continue;
    }
    pRecord->oldSiteIndexOfNewSites[counter++] = i;
  }
  pRecord->oldSiteIndexOfNewSites[pRecord->newSiteCount-1] = -1;
  pRecord->rotCountOfUnifiedNewSite = EnergyMatrixGetRotamerCount(pEnergyMatrix,unifySite[0]) *
                    EnergyMatrixGetRotamerCount(pEnergyMatrix,unifySite[1]);
  pRecord->oldRotIndexFromSite1 = (int*)malloc(sizeof(int)*pRecord->rotCountOfUnifiedNewSite);
  pRecord->oldRotIndexFromSite2 = (int*)malloc(sizeof(int)*pRecord->rotCountOfUnifiedNewSite);

  counter = 0;
  for(i=0;i<pList->rotamerCount[unifySite[0]];i++){
    if(pList->remainFlag[unifySite[0]][i]==FALSE){
      continue;
    }
    for(j=0;j<pList->rotamerCount[unifySite[1]];j++){
      if(pList->remainFlag[unifySite[1]][j]==FALSE){
        continue;
      }
      pRecord->oldRotIndexFromSite1[counter] = i;
      pRecord->oldRotIndexFromSite2[counter] = j;
      counter++;
    }
  }


  //Unify these two sites
  DEENodeUnifyUpdateRotamerList(pList,pRecord);
  DEENodeUnifyUpdateEnergyMatrix(pEnergyMatrix,pRecord);
  DEENodeUnifyUpdateRemainFlag(pRemainFlag,pRecord);

  return Success;
}

int DEENodeDeUnify(RotamerList* pList,DEENodeUnifyRecord* pRecord){
  int i,j;
  int unifySite[2];
  RotamerList listBeforeUnify;
  listBeforeUnify.designSiteCount = pList->designSiteCount+1;
  listBeforeUnify.rotamerCount = (int*)malloc(sizeof(int)*listBeforeUnify.designSiteCount);
  listBeforeUnify.remainFlag = (BOOL**)malloc(sizeof(BOOL*)*listBeforeUnify.designSiteCount);

  for(i=0;i<pList->designSiteCount-1;i++){
    int indexBeforeUnify = pRecord->oldSiteIndexOfNewSites[i];
    listBeforeUnify.rotamerCount[indexBeforeUnify] = pList->rotamerCount[i];
    listBeforeUnify.remainFlag[indexBeforeUnify] = (BOOL*)malloc(
      sizeof(BOOL)*pList->rotamerCount[i]);
    memcpy(listBeforeUnify.remainFlag[indexBeforeUnify],
      pList->remainFlag[i],sizeof(BOOL)*pList->rotamerCount[i]);
  }

  for(i=0;i<2;i++){
    unifySite[i] = pRecord->whichTwoOldSitesHaveBeenUnified[i];
    listBeforeUnify.rotamerCount[unifySite[i]] = pRecord->totalRotCountOnOldSites[i];
    listBeforeUnify.remainFlag[unifySite[i]] = (BOOL*)malloc(sizeof(BOOL)*pRecord->totalRotCountOnOldSites[i]);
    for(j=0;j<pRecord->totalRotCountOnOldSites[i];j++){
      listBeforeUnify.remainFlag[unifySite[i]][j] = FALSE;
    }
    for(j=0;j<pRecord->rotCountOfUnifiedNewSite;j++){
      int oldRotIndex = i==0? pRecord->oldRotIndexFromSite1[j] : pRecord->oldRotIndexFromSite2[j];
      if(pList->remainFlag[pList->designSiteCount-1][j]==TRUE){
        listBeforeUnify.remainFlag[unifySite[i]][ oldRotIndex ] = TRUE;
      }
    }
  }

  RotamerListDestroy(pList);
  pList->designSiteCount = listBeforeUnify.designSiteCount;
  pList->remainFlag = listBeforeUnify.remainFlag;
  pList->rotamerCount = listBeforeUnify.rotamerCount;
  return Success;
}

int DEE(EnergyMatrix* pEnergyMatrix,RotamerList* pRotamerList,double deeTreshold){
  int i;
  DEENodeUnifyRecord* unifyRecords;
  int unifyTimes;
  EnergyMatrix remainFlag;
  IntArray deleteRotamerList;
  IntArrayCreate(&deleteRotamerList,0);
  
  EnergyMatrixCreate(&remainFlag);
  EnergyMatrixCopy(&remainFlag,pEnergyMatrix);
  for(i=0;i<EnergyMatrixGetSiteCount(&remainFlag);i++){
    int k;
    for(k=0;k<EnergyMatrixGetSiteCount(&remainFlag);k++){
      int j,s;
      for(j=0;j<EnergyMatrixGetRotamerCount(&remainFlag,i);j++){
        for(s=0;s<EnergyMatrixGetRotamerCount(&remainFlag,k);s++){
          *EnergyMatrixGet(&remainFlag,i,k,j,s) = 1.0;
        }
      }
    }
  }


  unifyTimes = 0;
  unifyRecords = (DEENodeUnifyRecord*)malloc(sizeof(DEENodeUnifyRecord)*pRotamerList->designSiteCount);
  for(i=0;i<pRotamerList->designSiteCount;i++){
    DEENodeUnifyRecordCreate(&unifyRecords[i]);
  }

  do{
    int eliminatedCountThisLoop;
    do{
      eliminatedCountThisLoop = 0;
      do{
        //EnergyMatrixShow(pEnergyMatrix);
        printf("DEEGoldstein...");
        DEEGoldstein(pEnergyMatrix,&deleteRotamerList,DEE_PARA_THRESHOLD,&remainFlag);
        DEERotamerListAndEnergyMatrixDelete(pRotamerList,pEnergyMatrix,&remainFlag,&deleteRotamerList);
        printf("%d\n",IntArrayGetLength(&deleteRotamerList)/2);
        eliminatedCountThisLoop += IntArrayGetLength(&deleteRotamerList)/2;
      }while(IntArrayGetLength(&deleteRotamerList)>0);

      do{
        //EnergyMatrixShow(pEnergyMatrix);
        printf("DEESplit...");
        DEESplit(pEnergyMatrix,&deleteRotamerList,DEE_PARA_THRESHOLD,&remainFlag);
        DEERotamerListAndEnergyMatrixDelete(pRotamerList,pEnergyMatrix,&remainFlag,&deleteRotamerList);
        printf("%d\n",IntArrayGetLength(&deleteRotamerList)/2);
        eliminatedCountThisLoop += IntArrayGetLength(&deleteRotamerList)/2;

      }while(IntArrayGetLength(&deleteRotamerList)>0);
      
      if(eliminatedCountThisLoop < 10){
        break;
      }

      //EnergyMatrixShow(pEnergyMatrix);
      printf("DEEDouble...");
      i = DEEDouble(pEnergyMatrix,deeTreshold,&remainFlag);
      printf("%d\n",i);

    }while(TRUE);

    if(pRotamerList->designSiteCount>1){
      printf("Node Unifying...");
      DEENodeUnify(pEnergyMatrix,&remainFlag,pRotamerList,&unifyRecords[unifyTimes++]);
      printf(" Site %d and %d -> %d\n",
        unifyRecords[unifyTimes-1].whichTwoOldSitesHaveBeenUnified[0],
        unifyRecords[unifyTimes-1].whichTwoOldSitesHaveBeenUnified[1],
        pRotamerList->designSiteCount-1);
    }
  }while(EnergyMatrixGetSiteCount(pEnergyMatrix)>1 || EnergyMatrixGetRotamerCount(pEnergyMatrix,0)>1 );


  while(unifyTimes>0){
    DEENodeDeUnify(pRotamerList,&unifyRecords[--unifyTimes]);
  }

  RotamerListShow(pRotamerList);
  printf("Min Energy : %f\n",DEECalcMinEnergy(pEnergyMatrix,pRotamerList));

  IntArrayDestroy(&deleteRotamerList);
  for(i=0;i<pRotamerList->designSiteCount;i++){
    DEENodeUnifyRecordDestroy(&unifyRecords[i]);
  }
  free(unifyRecords);

  EnergyMatrixDestroy(&remainFlag);
  return Success;
}


int DEEShowDesignedStructure(Structure* pStructure, RotamerList* pRotamerList, char* structureFile){
  FILE* pf = fopen(structureFile, "w");
  if(pf==NULL){
    printf("In %s %s %d, cannot open file %s for writing\n", __FILE__, __FUNCTION__, __LINE__, structureFile);
    return IOError;
  }
  int atomIndex = 1;
  for(int i=0; i<StructureGetDesignSiteCount(pStructure); i++){
    DesignSite* pDesignSite = StructureGetDesignSite(pStructure,i);
    RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
    for(int j = 0; j < pRotamerList->rotamerCount[i]; j++){
      if(pRotamerList->remainFlag[i][j]==TRUE){
        Rotamer* pRotamer = RotamerSetGet(pRotamerSet, j);
        RotamerRestore(pRotamer, pRotamerSet);
        RotamerShowInPDBFormat(pRotamer, "ATOM", RotamerGetChainName(pRotamer),atomIndex,RotamerGetPosInChain(pRotamer),FALSE, pf);
        RotamerExtract(pRotamer);
        atomIndex++;
      }
    }
    
  }
  if(pf!=NULL){
    fclose(pf);
  }
  return Success;
}
