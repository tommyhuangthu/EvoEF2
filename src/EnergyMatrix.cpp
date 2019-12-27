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

#pragma warning(disable:4090)
#pragma warning(disable:4101)
#include "EnergyMatrix.h"
#include <string.h>
#include <time.h>

//functions for energy block
int EnergyMatrixBlockCreate(EnergyMatrixBlock* pThis)
{
  pThis->DesignSiteI = 0;
  pThis->DesignSiteK = 0;
  pThis->RotamerCountSiteI = 0;
  pThis->RotamerCountSiteK = 0;
  pThis->energyIK = NULL;

  return Success;
}
void        EnergyMatrixBlockDestroy(EnergyMatrixBlock* pThis){
  if(pThis->energyIK!=NULL){
    free(pThis->energyIK);
    pThis->energyIK = NULL;
  }
}
int         EnergyMatrixBlockCopy(EnergyMatrixBlock* pThis,EnergyMatrixBlock* pOther){
  EnergyMatrixBlockDestroy(pThis);
  EnergyMatrixBlockCreate(pThis);
  pThis->DesignSiteI = pOther->DesignSiteI;
  pThis->DesignSiteK = pOther->DesignSiteK;
  pThis->RotamerCountSiteI = pOther->RotamerCountSiteI;
  pThis->RotamerCountSiteK = pOther->RotamerCountSiteK;
  pThis->energyIK = (double*)malloc(sizeof(double)*pThis->RotamerCountSiteI*pThis->RotamerCountSiteK);
  memcpy(pThis->energyIK, pOther->energyIK, sizeof(double)*pThis->RotamerCountSiteI*pThis->RotamerCountSiteK);

  return Success;
}

double* EnergyMatrixBlockGet(EnergyMatrixBlock* pThis,int J,int S){
  return &pThis->energyIK[J*pThis->RotamerCountSiteK + S];
}


int EnergyMatrixBlockGenerate(EnergyMatrixBlock* pThis,Structure *pStructure,int designSiteI,int designSiteK,FILE* outputFile){
  DesignSite *pDesignSiteI = StructureGetDesignSite(pStructure, designSiteI);
  DesignSite *pDesignSiteK = StructureGetDesignSite(pStructure, designSiteK);
  if(designSiteI>designSiteK || pDesignSiteI==NULL || pDesignSiteK==NULL){
    return ValueError;
  }

  pThis->DesignSiteI = designSiteI;
  pThis->DesignSiteK = designSiteK;
  pThis->RotamerCountSiteI = RotamerSetGetCount(DesignSiteGetRotamers(pDesignSiteI));
  pThis->RotamerCountSiteK = RotamerSetGetCount(DesignSiteGetRotamers(pDesignSiteK));
  pThis->energyIK = (double*)malloc(sizeof(double)*pThis->RotamerCountSiteI*pThis->RotamerCountSiteK);

  if(designSiteI==designSiteK){
    for(int j=0;j<pThis->RotamerCountSiteI;j++){
      Rotamer *pRotamerIJ  = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteI), j);
      RotamerRestore(pRotamerIJ,DesignSiteGetRotamers(pDesignSiteI));
      double* pEnergy = EnergyMatrixBlockGet(pThis,j,j);
      *pEnergy = 0.0;
      double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
      EnergyTermInitialize(energyTerms);
      Chain *pChainI = StructureGetChain(pStructure, pDesignSiteI->chainIndex);
      if(ChainGetType(pChainI) == Type_Chain_Protein){
        EVOEF_AminoAcidReferenceEnergy(RotamerGetType(pRotamerIJ), energyTerms);
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){//same chain
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(b==pDesignSiteI->resiIndex){
                EVOEF_EnergyProteinRotamerIntraEnergy(pRotamerIJ, energyTerms);
              }
              else{
                if(ChainGetType(pChainA) != Type_Chain_SmallMol){
                  if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                    EVOEF_EnergyProteinRotamerAndFixedResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                  else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                    pResidueAB->designSiteType == Type_ResidueDesignType_Mutated){
                      EVOEF_EnergyProteinRotamerAndDesignResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                }
              }
            }
          }
          else{//different chains
            if(ChainGetType(pChainA) != Type_Chain_SmallMol){
              for(int b=0; b<ChainGetResidueCount(pChainA); b++){
                Residue* pResidueAB = ChainGetResidue(pChainA, b);
                if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                  EVOEF_EnergyProteinRotamerAndFixedResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
                }
                else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                  pResidueAB->designSiteType == Type_ResidueDesignType_Mutated){
                    EVOEF_EnergyProteinRotamerAndDesignResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
                }
              }
            }
          }
        }
      }
      else if(ChainGetType(pChainI) == Type_Chain_SmallMol){
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){
            //*pEnergy += pRotamerIJ->vdwInternal;
          }
          else{
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                EVOEF_EnergyProteinRotamerAndFixedResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
              }
              else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                pResidueAB->designSiteType == Type_ResidueDesignType_Mutated){
                  EVOEF_EnergyProteinRotamerAndDesignResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
              }
            }
          }
        }
      }

      EnergyTermWeighting(energyTerms);
      *pEnergy += energyTerms[0];
      fprintf(outputFile,"%f %d %d %d %d\n",*pEnergy,designSiteI,j,designSiteK,j);
      RotamerExtract(pRotamerIJ);
    }
    return Success;
  }

  for(int j=0; j<pThis->RotamerCountSiteI; j++){
    Rotamer *pRotamerIJ = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteI), j);
    RotamerRestore(pRotamerIJ,DesignSiteGetRotamers(pDesignSiteI));
    for(int s=0; s<pThis->RotamerCountSiteK; s++){
      Rotamer *pRotamerKS = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteK), s);
      RotamerRestore(pRotamerKS,DesignSiteGetRotamers(pDesignSiteK));
      double* pEnergy = EnergyMatrixBlockGet(pThis,j,s);
      *pEnergy = 0.0;

      double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
      EnergyTermInitialize(energyTerms);
      if(pDesignSiteI->chainIndex == pDesignSiteK->chainIndex){
        EVOEF_EnergyProteinRotamerAndRotamerSameChain(pRotamerIJ, pRotamerKS, energyTerms);
      }
      else{
        EVOEF_EnergyProteinRotamerAndRotamerDifferentChain(pRotamerIJ, pRotamerKS, energyTerms);
      }
      EnergyTermWeighting(energyTerms);
      for(int a=1; a<MAX_EVOEF_ENERGY_TERM_NUM;a++){
        *pEnergy += energyTerms[a];
      }
      fprintf(outputFile,"%f %d %d %d %d\n",*pEnergy,designSiteI,j,designSiteK,s);
      RotamerExtract(pRotamerKS);
    }
    RotamerExtract(pRotamerIJ);
  }
  return Success;
}

//functions for EnergyMatrix
int EnergyMatrixCreate(EnergyMatrix* pThis){
  pThis->designSiteCount = 0;
  pThis->blocks = NULL;
  return Success;
}

void EnergyMatrixDestroy(EnergyMatrix* pThis){
  int i,j;
  for(i=0;i<pThis->designSiteCount;i++){
    for(j=i;j<pThis->designSiteCount;j++){
      EnergyMatrixBlockDestroy(&pThis->blocks[ i*pThis->designSiteCount + j ]);
    }
  }
  free(pThis->blocks);
  pThis->blocks = NULL;
}
int EnergyMatrixCopy(EnergyMatrix* pThis,EnergyMatrix* pOther){
  int i,j;
  EnergyMatrixDestroy(pThis);
  EnergyMatrixCreate(pThis);

  pThis->designSiteCount = pOther->designSiteCount;
  pThis->blocks = (EnergyMatrixBlock*)malloc(sizeof(EnergyMatrixBlock)*
    pThis->designSiteCount * pThis->designSiteCount);

  for(i=0;i<pThis->designSiteCount;i++){
    for(j=i;j<pThis->designSiteCount;j++){
      EnergyMatrixBlockCreate(EnergyMatrixGetBlock(pThis,i,j));
      EnergyMatrixBlockCopy(EnergyMatrixGetBlock(pThis,i,j),EnergyMatrixGetBlock(pOther,i,j));
    }
  }
  return Success;
}

int  EnergyMatrixRead(EnergyMatrix* pThis,char* filepath){
  int i,j;
  int result;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  int fileLength;
  FILE* pFile = fopen(filepath,"r");
  if(pFile==NULL){
    result = IOError;
    sprintf(errMsg,"in file %s function %s line %d, cannot open file %s to read",__FILE__,__FUNCTION__,__LINE__,filepath);
    TraceError(errMsg,result);
    return result;
  }

  fseek(pFile,0,SEEK_END);
  fileLength = ftell(pFile);
  fseek(pFile,0,SEEK_SET);

  printf("Reading EnergyMatrix File %s\n",filepath);

  while( fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pFile)){
    double energy;
    int designSiteI,designSiteK,rotamerIJ,rotamerKS;
    EnergyMatrixBlock* energyIK;

    sscanf(line,"%lf %d %d %d %d",&energy,&designSiteI,&rotamerIJ,&designSiteK,&rotamerKS);

    if(designSiteI<0 || designSiteK<0 || rotamerIJ<0 || rotamerKS<0 || designSiteI>designSiteK){
      result = FormatError;
      sprintf(errMsg,"in file %s function %s line %d, bad format at this line of file %s\n%s",__FILE__,__FUNCTION__,__LINE__,filepath,line);
      TraceError(errMsg,result);
      return result;
    }

    if(designSiteI>=pThis->designSiteCount || designSiteK>=pThis->designSiteCount){
      int newDesignSiteCount = designSiteI>designSiteK? designSiteI+1 : designSiteK+1;
      EnergyMatrixBlock* newBlocks = (EnergyMatrixBlock*)malloc(
        sizeof(EnergyMatrixBlock) * newDesignSiteCount * newDesignSiteCount);
      for(i=0;i<newDesignSiteCount;i++){
        for(j=0;j<newDesignSiteCount;j++){
          EnergyMatrixBlockCreate(&newBlocks[i*newDesignSiteCount+j]);
        }
      }
      for(i=0;i<pThis->designSiteCount;i++){
        for(j=i;j<pThis->designSiteCount;j++){
          newBlocks[i*newDesignSiteCount+j] = *EnergyMatrixGetBlock(pThis,i,j);
        }
      }

      pThis->designSiteCount = newDesignSiteCount;
      free(pThis->blocks);
      pThis->blocks = newBlocks;
    }

    energyIK = EnergyMatrixGetBlock(pThis,designSiteI,designSiteK);
    if(energyIK->energyIK != NULL){
      result = AssertionError;
      sprintf(errMsg,"in file %s function %s line %d, energy block between designSite %d and %d has already "
        "been read, when reading file %s at this line:\n%s",__FILE__,__FUNCTION__,__LINE__,designSiteI,designSiteK,filepath,line);
      TraceError(errMsg,result);
      return result;
    }
    if(rotamerIJ >= energyIK->RotamerCountSiteI){
      energyIK->RotamerCountSiteI = rotamerIJ+1;
    }
    if(rotamerKS >= energyIK->RotamerCountSiteK){
      energyIK->RotamerCountSiteK = rotamerKS+1;
    }

  }

  //All energy blocks whose (RotamerCountSiteI!=0 && RotamerCountSiteJ!=0 && energyIK==NULL) are updated this time

  for(i=0;i<pThis->designSiteCount;i++){
    for(j=i;j<pThis->designSiteCount;j++){
      int k;
      EnergyMatrixBlock* pEnergyIK = EnergyMatrixGetBlock(pThis,i,j);
      if(pEnergyIK->RotamerCountSiteI!=0 && pEnergyIK->RotamerCountSiteK!=0 && pEnergyIK->energyIK==NULL){
        pEnergyIK->energyIK = (double*)malloc(
          sizeof(double)*pEnergyIK->RotamerCountSiteI*pEnergyIK->RotamerCountSiteK);
      }
      for(k=0;k<pEnergyIK->RotamerCountSiteI*pEnergyIK->RotamerCountSiteK;k++){
        pEnergyIK->energyIK[k] = 0.0;
      }
      pEnergyIK->DesignSiteI = i;
      pEnergyIK->DesignSiteK = j;
    }
  }

  //read for a second time
  fseek(pFile,0,SEEK_SET);
  while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pFile)){
    double energy;
    int designSiteI,designSiteK,rotamerIJ,rotamerKS;
    sscanf(line,"%lf %d %d %d %d",&energy,&designSiteI,&rotamerIJ,&designSiteK,&rotamerKS);
    *EnergyMatrixGet(pThis,designSiteI,designSiteK,rotamerIJ,rotamerKS) = energy;
  }

  fclose(pFile);

  return Success;
}

int  EnergyMatrixCheck(char* filepath){
  int MAX_DESIGN_COUNT = 1000;
  int siteCount;
  int rotamerCount[1000];
  int i,j,k,s;
  double value;
  int lineCounter;
  FILE* pFile;
  int fileLength;

  siteCount = 0;
  for(i=0; i<1000; i++) rotamerCount[i] = 0;

  pFile = fopen(filepath,"r");
  if(pFile==NULL){
    printf("Cannot open file %s\n",filepath);
    return IOError;
  }

  fseek(pFile,0,SEEK_END);
  fileLength = ftell(pFile);
  fseek(pFile,0,SEEK_SET);

  lineCounter = 0;
  printf("Counting design sites and rotamers...\n");
  while( fscanf(pFile,"%lf %d %d %d %d",&value,&i,&j,&k,&s)!=EOF ){
    lineCounter++;
    if(lineCounter%100000 == 0){
      printf("\r");ShowProgress(40,100.0*ftell(pFile)/fileLength);
    }
    if(i<0 || i>=MAX_DESIGN_COUNT || k<0 || k>=MAX_DESIGN_COUNT){
      printf("Error at this line#%d, i and k must between 0 and %d :\n%f %d %d %d %d\n",
        lineCounter,MAX_DESIGN_COUNT,value,i,j,k,s);
      return FormatError;
    }
    if(i>=siteCount){
      siteCount = i+1;
    }
    if(k>=siteCount){
      siteCount = k+1;
    }
    if(j>=rotamerCount[i]){
      rotamerCount[i] = j+1;
    }
    if(s>=rotamerCount[k]){
      rotamerCount[k] = s+1;
    }
  }
  printf("\n totally %d sites, %d lines in file\n",siteCount,lineCounter);;
  for(i=0;i<siteCount;i++){
    printf("Design site #%d, %d rotamers\n",i,rotamerCount[i]);
  }

  fseek(pFile,0,SEEK_SET);
  lineCounter = 0;
  for(i=0;i<siteCount;i++){
    for(k=i;k<siteCount;k++){
      printf("Checking energy block i = %d, k = %d, size= %d*%d...",
        i,k,rotamerCount[i],rotamerCount[k]);
      for(j=0;j<rotamerCount[i];j++){
        for(s=0;s<rotamerCount[k];s++){
          int iread,jread,kread,sread,result;
          if(i==k && j!=s){
            continue;
          }

          lineCounter++;
          result = fscanf(pFile,"%lf %d %d %d %d",&value,&iread,&jread,&kread,&sread);
          if(result == EOF){
            printf("\nAt the end of file, energy item i=%d, j=%d, k=%d, s=%d is missing\n",
              i,j,k,s);
            return FormatError;
          }
          if(iread!=i || jread!=j || kread!=k || sread!=s){
            printf("\nAt line #%d, expect i=%d, j=%d, k=%d,s=%d, but this line is read:\n"
              "%f %d %d %d %d",lineCounter,i,j,k,s,value,iread,jread,kread,sread);
            return FormatError;
          }
          if(fabs(value)>1e8){
            printf("\nAt line #%d, energy value is too large:\n"
              "%f %d %d %d %d",lineCounter,value,i,j,k,s);
            return FormatError;
          }
        }
      }
      printf("OK\n");
    }
  }
  printf("all blocks of energy matrix file %s are okay\n",filepath);

  fclose(pFile);
  return Success;

}

int EnergyMatrixWrite(EnergyMatrix* pThis,char* filepath){
  int designSiteI,designSiteK;
  int result;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  FILE* pFile = NULL;

  if(filepath==NULL){
    pFile = stdout;
  }else{
    pFile = fopen(filepath,"w");
  }

  if(pFile==NULL){
    result = IOError;
    sprintf(errMsg,"in file %s function %s line %d, cannot open file %s to write",__FILE__,__FUNCTION__,__LINE__,filepath);
    TraceError(errMsg,result);
    return result;
  }

  printf("Write EnergyMatrix file %s\n",filepath);
  for(designSiteI=0;designSiteI<pThis->designSiteCount;designSiteI++){
    for(designSiteK=designSiteI;designSiteK<pThis->designSiteCount;designSiteK++){
      EnergyMatrixBlock* pEnergyBlock = EnergyMatrixGetBlock(pThis,designSiteI,designSiteK);
      int rotamerIJ,rotamerKS;
      for(rotamerIJ=0; rotamerIJ<pEnergyBlock->RotamerCountSiteI; rotamerIJ++){
        for(rotamerKS=0; rotamerKS<pEnergyBlock->RotamerCountSiteK; rotamerKS++){
          if(designSiteI==designSiteK && rotamerIJ!=rotamerKS){
            continue;
          }
          fprintf(pFile,"%f  %d  %d  %d  %d\n",
            *EnergyMatrixBlockGet(pEnergyBlock,rotamerIJ,rotamerKS),
            designSiteI,rotamerIJ,designSiteK,rotamerKS);
        }
      }
    }
  }

  if(filepath!=NULL){
    fclose(pFile);
  }

  return Success;
}

EnergyMatrixBlock* EnergyMatrixGetBlock(EnergyMatrix* pThis,int designSiteI,int designSiteK){
  return &pThis->blocks[designSiteI*pThis->designSiteCount + designSiteK];
}

double* EnergyMatrixGet(EnergyMatrix* pThis,int designSiteI,int designSiteK,int rotamerIJ,int rotamerKS){
  EnergyMatrixBlock* pEnerygIK;

  if(designSiteI>designSiteK){
    int swap = designSiteI;
    designSiteI = designSiteK;
    designSiteK = swap;

    swap = rotamerIJ;
    rotamerIJ = rotamerKS;
    rotamerKS = swap;
  }
  pEnerygIK = EnergyMatrixGetBlock(pThis,designSiteI,designSiteK);
  return EnergyMatrixBlockGet(pEnerygIK,rotamerIJ,rotamerKS);
}


int  EnergyMatrixGetSiteCount(EnergyMatrix* pThis){
  return pThis->designSiteCount;
}

int  EnergyMatrixGetRotamerCount(EnergyMatrix* pThis,int designSiteI){
  return pThis->blocks[designSiteI*pThis->designSiteCount+designSiteI].RotamerCountSiteI;
}

int EnergyMatrixShow(EnergyMatrix* pThis){
  int i,k,j,s;
  for(i=0;i<pThis->designSiteCount;i++){
    for(k=i;k<pThis->designSiteCount;k++){
      printf("Design Site %d and %d:\n",i,k);
      for(j=0; j<EnergyMatrixGetRotamerCount(pThis,i);j++){
        for(s=0; s<EnergyMatrixGetRotamerCount(pThis,k);s++){
          if(i==k && j!=s){
            continue;
          }
          printf("%6.3f ",*EnergyMatrixGet(pThis,i,k,j,s));
        }
        printf("\n");
      }
    }
  }
  return Success;
}



//'jobIndex' and 'totalJobCount' is designed for parallel computation,
//If not used in parallel computation environment, set both 'totalJobCount' and 'jobIndex' to be 1
int EnergyMatrixGenerate(Structure* pStructure,char* energyMatrixFilePath,int slotIndex,int slotCount){
  typedef struct{
    int computationAmount;
    int designSiteI;
    int designSiteK;
  } Job;

  typedef struct{
    Job* jobs;
    int jobCount;
    int totalComputationAmount;
  } Slot;

  int jobCount = pStructure->designSiteCount * (pStructure->designSiteCount+1) / 2;
  Job* jobs = (Job*)malloc(sizeof(Job)*jobCount);
  Slot* slots= (Slot*)malloc(sizeof(Slot)*slotCount);

  int i,j,counter;
  for(i=0;i<slotCount;i++){
    slots[i].jobCount = 0;
    slots[i].totalComputationAmount = 0;
    slots[i].jobs = NULL;
  }

  if(slotCount<=0 || slotIndex<=0 || slotIndex>slotCount){
    return IndexError;
  }
  slotIndex--;

  counter = 0;
  for(i=0;i<pStructure->designSiteCount;i++){
    for(j=i;j<pStructure->designSiteCount;j++){
      int rotamerCountOnSiteI = RotamerSetGetCount(&pStructure->designSites[i].rotamers);
      int rotamerCountOnSiteK = RotamerSetGetCount(&pStructure->designSites[j].rotamers);
      if(i==j){
        jobs[counter].computationAmount = rotamerCountOnSiteI;
      }
      else{
        jobs[counter].computationAmount = rotamerCountOnSiteI * rotamerCountOnSiteK;
      }
      jobs[counter].designSiteI = i;
      jobs[counter].designSiteK = j;
      counter++;
    }
  }

  //Sort 'jobs' by descending order ; Bubble sorting
  for(i=0;i<jobCount-1;i++){
    for(j=jobCount-1;j>i;j--){
      if(jobs[j].computationAmount > jobs[j-1].computationAmount){
        Job tempJob = jobs[j];
        jobs[j] = jobs[j-1];
        jobs[j-1] = tempJob;
      }
    }
  }

  //Assigning jobs to each slot
  for(i=0;i<jobCount;i++){
    Slot* pMinSlot = &slots[0];
    int minSlotAmount = slots[0].totalComputationAmount;
    for(j=0;j<slotCount;j++){
      if(slots[j].totalComputationAmount < minSlotAmount){
        pMinSlot = &slots[j];
        minSlotAmount = slots[j].totalComputationAmount;
      }
    }
    (pMinSlot->jobCount)++;
    pMinSlot->jobs = (Job*)realloc(pMinSlot->jobs,sizeof(Job)* pMinSlot->jobCount );
    pMinSlot->jobs[pMinSlot->jobCount-1] = jobs[i];
    pMinSlot->totalComputationAmount += jobs[i].computationAmount;
  }

  //Show the job assignment of each slot
  for(i=0;i<slotCount;i++){
    printf("Slot #%d, %8d ",i+1,slots[i].totalComputationAmount);
    for(j=0;j<slots[i].jobCount; j++){
      printf("(%2d,%2d) ",slots[i].jobs[j].designSiteI,slots[i].jobs[j].designSiteK);
    }
    printf("\n");
  }

  //Do the computational jobs assigned to current slot
  for(i=slots[slotIndex].jobCount-1;i>=0;i--){ 
    int result;
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    Job curJob = slots[slotIndex].jobs[i];
    EnergyMatrixBlock energyBlock;
    FILE* energyBlockFile = NULL;
    char energyBlockFileName[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sprintf(energyBlockFileName,"%s_%3.3d_%3.3d.txt",energyMatrixFilePath,curJob.designSiteI,curJob.designSiteK);
    printf("%s\n",energyBlockFileName);
    energyBlockFile = fopen(energyBlockFileName,"w");
    if(energyBlockFile==NULL){
      sprintf(errMsg,"int file %s function %s line %d, cannot create file %s for writing",__FILE__,__FUNCTION__,__LINE__,
        energyBlockFileName);
      result = IOError;
      TraceError(errMsg,result);
      return result;        
    }


    EnergyMatrixBlockCreate(&energyBlock);

    result = EnergyMatrixBlockGenerate(&energyBlock,pStructure,
      curJob.designSiteI,curJob.designSiteK,energyBlockFile);
    if(FAILED(result)){
      return result;
    }

    EnergyMatrixBlockDestroy(&energyBlock);
    fclose(energyBlockFile);
  }

  for(i=0;i<slotCount;i++){
    free(slots[i].jobs);
  }
  free(jobs);
  free(slots);
  return Success;
}

//////////////////////////////////////////////////////////////////////////////////////////
//split the whole job of energy matrix computation into small jobs for parallel computing
//////////////////////////////////////////////////////////////////////////////////////////
int CompareJobAmount(const void *a, const void *b)
{
  double da = (*(Job*)a).computationAmount;
  double db = (*(Job*)b).computationAmount;

  return da < db ? -1 : 1;
}

int EnergyMatrixGenerateBasedOnPartition(Structure* pStructure, RotamerList* pList, char* energyMatrixFilePath,int slotIndex,int slotCount){
  int PARTITION_SIZE = PARTITION_SIZE_PARA;
  int jobCount;
  Job* jobs = NULL;
  Slot* slots= NULL;

  int i, j, k, s, counter;

  time_t totalTimeStart, totalTimeCurrent;
  int totalTimeElapsed;

  double cpuTimeStart, cpuTimeCurrent, cpuTimeElapsed;

  SitePartition* sitePartitions = (SitePartition*)malloc(sizeof(SitePartition)*StructureGetDesignSiteCount(pStructure));

  int * remainRotamerCount = (int*)malloc(sizeof(int)*pList->designSiteCount);
  double* energyInMemory = NULL;
  FILE* slotEnergyFile = NULL;
  char slotEnergyFileName[MAX_LENGTH_ONE_LINE_IN_FILE+1];

  for(i = 0; i < pList->designSiteCount; i++){
    remainRotamerCount[i] = 0;
    for(j = 0; j < pList->rotamerCount[i]; j++){
      if(pList->remainFlag[i][j] == TRUE) remainRotamerCount[i]++;
    }
  }
  for(i = 0; i < StructureGetDesignSiteCount(pStructure); i++){
    int rotamerIndex;
    DesignSite* pDesignSiteI = StructureGetDesignSite(pStructure, i);
    RotamerSet* pRotamerSetI = DesignSiteGetRotamers(pDesignSiteI);
    int partitionCountI = (int)(remainRotamerCount[i]/PARTITION_SIZE);
    if(remainRotamerCount[i]%PARTITION_SIZE != 0) partitionCountI++;
    printf("site: %d, partition count: %d\n", i, partitionCountI);
    sitePartitions[i].designSite = i;
    sitePartitions[i].partitionCount = partitionCountI;
    sitePartitions[i].rotamerCount = (int*)malloc(sizeof(int)*partitionCountI);
    sitePartitions[i].rotamerIndexOld = (IntArray*)malloc(sizeof(IntArray)*partitionCountI);
    sitePartitions[i].rotamerIndexNew = (IntArray*)malloc(sizeof(IntArray)*partitionCountI);

    for(j = 0; j < partitionCountI; j++){
      IntArrayCreate(&sitePartitions[i].rotamerIndexOld[j], 0);
      IntArrayCreate(&sitePartitions[i].rotamerIndexNew[j], 0);
      if(j == partitionCountI - 1){
        sitePartitions[i].rotamerCount[j] = remainRotamerCount[i] - j*PARTITION_SIZE;
        
      }
      else{
        sitePartitions[i].rotamerCount[j] = PARTITION_SIZE;
      }
      IntArrayResize(&sitePartitions[i].rotamerIndexOld[j], sitePartitions[i].rotamerCount[j]);
      IntArrayResize(&sitePartitions[i].rotamerIndexNew[j], sitePartitions[i].rotamerCount[j]);
    }

    rotamerIndex = 0;
    for(j = 0; j < RotamerSetGetCount(pRotamerSetI); j++){
      int partitionIndex;
      if(pList->remainFlag[i][j] == FALSE) continue;
      partitionIndex = rotamerIndex/PARTITION_SIZE;
      IntArraySet(&sitePartitions[i].rotamerIndexOld[partitionIndex], rotamerIndex - partitionIndex*PARTITION_SIZE, j);
      IntArraySet(&sitePartitions[i].rotamerIndexNew[partitionIndex], rotamerIndex - partitionIndex*PARTITION_SIZE, rotamerIndex);
      rotamerIndex++;
    }
  }

  jobCount = 0;
  for(i = 0; i < pStructure->designSiteCount; i++){
    int partitionCountOnSiteI = sitePartitions[i].partitionCount;
    for(k = i; k < pStructure->designSiteCount; k++){
      int partitionCountOnSiteK = sitePartitions[k].partitionCount;
      if(k == i){
        jobCount += partitionCountOnSiteI;
      }
      else{
        jobCount += partitionCountOnSiteI * partitionCountOnSiteK;
      }
    }
  }
  printf("jobCount: %d\n", jobCount);
  jobs = (Job*)malloc(sizeof(Job)*jobCount);
  slots= (Slot*)malloc(sizeof(Slot)*slotCount);


  for(i = 0; i < slotCount; i++){
    slots[i].jobCount = 0;
    slots[i].totalComputationAmount = 0.0;
    slots[i].jobs = NULL;
  }

  if(slotCount <= 0 || slotIndex <= 0 || slotIndex > slotCount){
    return ValueError;
  }
  slotIndex--;

  counter = 0;
  for(i = 0; i < pStructure->designSiteCount; i++){
    Chain* pChainI = StructureGetChain(pStructure, StructureGetDesignSite(pStructure, i)->chainIndex);
    int partitionCountOnSiteI = sitePartitions[i].partitionCount;
    for(k = i; k < pStructure->designSiteCount; k++){
      Chain* pChainK = StructureGetChain(pStructure, StructureGetDesignSite(pStructure, k)->chainIndex);
      int partitionCountOnSiteK = sitePartitions[k].partitionCount;
      for(j = 0; j < partitionCountOnSiteI; j++){
        if(k == i){
          jobs[counter].computationAmount = sitePartitions[i].rotamerCount[j];

          jobs[counter].designSiteI = i;
          jobs[counter].designSiteK = i;
          jobs[counter].rotamerPartitionIndexOnSiteI = j;
          jobs[counter].rotamerPartitionIndexOnSiteK = j;
          counter++;
        }
        else{
          for(s = 0; s < partitionCountOnSiteK; s++){
            jobs[counter].computationAmount = sitePartitions[i].rotamerCount[j] * sitePartitions[k].rotamerCount[s];
            
            jobs[counter].designSiteI = i;
            jobs[counter].designSiteK = k;
            jobs[counter].rotamerPartitionIndexOnSiteI = j;
            jobs[counter].rotamerPartitionIndexOnSiteK = s;
            counter++;
          }
        }
      }
    }
  }

  qsort(jobs, jobCount, sizeof(Job), CompareJobAmount);
  for(i = jobCount-1; i >= 0; i--){
    Slot* pMinSlot = &slots[0];
    double minSlotAmount = slots[0].totalComputationAmount;
    for(j = 0; j < slotCount; j++){
      if(slots[j].totalComputationAmount < minSlotAmount){
        pMinSlot = &slots[j];
        minSlotAmount = slots[j].totalComputationAmount;
      }
    }
    (pMinSlot->jobCount)++;
    pMinSlot->jobs = (Job*)realloc(pMinSlot->jobs,sizeof(Job)* pMinSlot->jobCount );
    pMinSlot->jobs[pMinSlot->jobCount-1] = jobs[i];
    pMinSlot->totalComputationAmount += jobs[i].computationAmount;
  }

  //Show the job assignment of each slot
  for(i=0;i<slotCount;i++){
    int index = 0;
    printf("Slot #%d, computational amount %.2f \n",i+1,slots[i].totalComputationAmount);
    for(j=0;j<slots[i].jobCount; j++){
      printf("(%2d,%2d,%4d,%4d) ",slots[i].jobs[j].designSiteI, slots[i].jobs[j].designSiteK,
        slots[i].jobs[j].rotamerPartitionIndexOnSiteI, slots[i].jobs[j].rotamerPartitionIndexOnSiteK);
      index++;
      if(index % 5 == 0){
        printf("\n");
      }
    }
    printf("\n");
  }


  //allocate memories to store the energy temporarily for each CPU, and finally output to files;
  counter = 0;
  for(i = slots[slotIndex].jobCount - 1; i >= 0; i--){
    int designSiteI;
    int designSiteK;
    int partitionIndexOnSiteI;
    int partitionIndexOnSiteK;
    Job curJob = slots[slotIndex].jobs[i];

    designSiteI = curJob.designSiteI;
    designSiteK = curJob.designSiteK;
    partitionIndexOnSiteI = curJob.rotamerPartitionIndexOnSiteI;
    partitionIndexOnSiteK = curJob.rotamerPartitionIndexOnSiteK;

    counter += sitePartitions[designSiteI].rotamerCount[partitionIndexOnSiteI]*sitePartitions[designSiteK].rotamerCount[partitionIndexOnSiteK];
  }

  energyInMemory = (double*)malloc(sizeof(double)*counter);
  for(i = 0; i < counter; i++){
    energyInMemory[i] = 0.0;
  }

  //Do the computational jobs assigned to current slot
  totalTimeStart = time(NULL);
  cpuTimeStart = clock();
  counter = 0;
  for(i = slots[slotIndex].jobCount - 1; i >= 0; i--){
    int result;
    int designSiteI;
    int designSiteK;
    int partitionIndexOnSiteI;
    int partitionIndexOnSiteK;
    time_t begin;
    time_t end;
    int elapsed;
    Job curJob = slots[slotIndex].jobs[i];
    char energyBlockFileName[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    begin = time(NULL);
    sprintf(energyBlockFileName,"%s_%3.3d_%3.3d_%4.4d_%4.4d.txt",energyMatrixFilePath, 
      curJob.designSiteI, curJob.designSiteK, curJob.rotamerPartitionIndexOnSiteI, curJob.rotamerPartitionIndexOnSiteK);
    printf("%s, ",energyBlockFileName);

    designSiteI = curJob.designSiteI;
    designSiteK = curJob.designSiteK;
    partitionIndexOnSiteI = curJob.rotamerPartitionIndexOnSiteI;
    partitionIndexOnSiteK = curJob.rotamerPartitionIndexOnSiteK;

    result = EnergyMatrixPartitionPairGenerateStoreInMemory(pStructure, designSiteI, designSiteK, 
      &sitePartitions[designSiteI].rotamerIndexOld[partitionIndexOnSiteI], &sitePartitions[designSiteK].rotamerIndexOld[partitionIndexOnSiteK], 
      &sitePartitions[designSiteI].rotamerIndexNew[partitionIndexOnSiteI], &sitePartitions[designSiteK].rotamerIndexNew[partitionIndexOnSiteK],
      energyInMemory, counter);
    counter += sitePartitions[designSiteI].rotamerCount[partitionIndexOnSiteI]*sitePartitions[designSiteK].rotamerCount[partitionIndexOnSiteK];
    if(FAILED(result)){
      return result;
    }

    end = time(NULL);
    elapsed = (int)(end - begin);
    printf("elapsed time: %d hr %d min %d sec\n", (int)(elapsed/3600), (int)(elapsed%3600/60), (int)(elapsed%3600%60));
  }

  // output energy values from memories to files;
  counter = 0;
  sprintf(slotEnergyFileName, "%s_CPU%03d.txt", energyMatrixFilePath, slotIndex+1);
  slotEnergyFile = fopen(slotEnergyFileName, "w");
  if(slotEnergyFile == NULL){
    printf("cannot open file %s for writing\n", slotEnergyFileName);
    return IOError;
  }
  for(i = slots[slotIndex].jobCount - 1; i >= 0; i--){
    int result;
    int designSiteI;
    int designSiteK;
    int partitionIndexOnSiteI;
    int partitionIndexOnSiteK;
    Job curJob = slots[slotIndex].jobs[i];

    designSiteI = curJob.designSiteI;
    designSiteK = curJob.designSiteK;
    partitionIndexOnSiteI = curJob.rotamerPartitionIndexOnSiteI;
    partitionIndexOnSiteK = curJob.rotamerPartitionIndexOnSiteK;

    result = EnergyMatrixPartitionPairWriteEnergyToFile(pStructure, designSiteI, designSiteK, 
      &sitePartitions[designSiteI].rotamerIndexOld[partitionIndexOnSiteI], &sitePartitions[designSiteK].rotamerIndexOld[partitionIndexOnSiteK], 
      &sitePartitions[designSiteI].rotamerIndexNew[partitionIndexOnSiteI], &sitePartitions[designSiteK].rotamerIndexNew[partitionIndexOnSiteK],
      energyInMemory, counter, slotEnergyFile);
    counter += sitePartitions[designSiteI].rotamerCount[partitionIndexOnSiteI]*sitePartitions[designSiteK].rotamerCount[partitionIndexOnSiteK];
    if(FAILED(result)){
      return result;
    }
  }
  fclose(slotEnergyFile);

  totalTimeCurrent = time(NULL);
  totalTimeElapsed = (int)(totalTimeCurrent - totalTimeStart);
  printf("elapsed time for calculating EnergyMatrixBlocks: %d hr %d min %d sec\n", totalTimeElapsed/3600, totalTimeElapsed%3600/60, totalTimeElapsed%3600%60);

  cpuTimeCurrent = clock();
  cpuTimeElapsed = cpuTimeCurrent - cpuTimeStart;
  printf("elapsed CPU time for calculating EnergyMatrixBlocks: %f secs\n", cpuTimeElapsed/CLOCKS_PER_SEC);


  for(i = 0; i < pStructure->designSiteCount; i++){
    free(sitePartitions[i].rotamerCount);
    for(j = 0; j < sitePartitions[i].partitionCount; j++){
      IntArrayDestroy(&sitePartitions[i].rotamerIndexOld[j]);
      IntArrayDestroy(&sitePartitions[i].rotamerIndexNew[j]);
    }
  }
  free(sitePartitions);

  for(i=0;i<slotCount;i++){
    free(slots[i].jobs);
  }
  free(jobs);
  free(slots);
  free(remainRotamerCount);
  free(energyInMemory);


  return Success;
}

int EnergyMatrixPartitionPairGenerate(Structure* pStructure, int designSiteI, int designSiteK, 
                     IntArray* pRotamerIndexOldOnSiteI, IntArray* pRotamerIndexOldOnSiteK, 
                     IntArray* pRotamerIndexNewOnSiteI, IntArray* pRotamerIndexNewOnSiteK, FILE* outputFile)
{
  DesignSite *pDesignSiteI = StructureGetDesignSite(pStructure, designSiteI);
  DesignSite *pDesignSiteK = StructureGetDesignSite(pStructure, designSiteK);
  Chain *pChainI = StructureGetChain(pStructure, pDesignSiteI->chainIndex);
  Chain* pChainK = StructureGetChain(pStructure, pDesignSiteK->chainIndex);
  double* partitionPairEnergy = NULL;

  if(designSiteI > designSiteK || pDesignSiteI == NULL || pDesignSiteK == NULL){
    return ValueError;
  }

  if(designSiteI == designSiteK){
    partitionPairEnergy = (double*)malloc(sizeof(double)*IntArrayGetLength(pRotamerIndexOldOnSiteI));
    for(int j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
      Rotamer *pRotamerIJ  = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteI), IntArrayGet(pRotamerIndexOldOnSiteI, j));
      RotamerRestore(pRotamerIJ,DesignSiteGetRotamers(pDesignSiteI));

      double total = 0;
      double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
      EnergyTermInitialize(energyTerms);
      if(ChainGetType(pChainI) == Type_Chain_Protein){
        EVOEF_AminoAcidReferenceEnergy(RotamerGetType(pRotamerIJ), energyTerms);
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){//same chain
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(b==pDesignSiteI->resiIndex){
                EVOEF_EnergyProteinRotamerIntraEnergy(pRotamerIJ, energyTerms);
              }
              else{
                if(ChainGetType(pChainA) != Type_Chain_SmallMol){
                  if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                    EVOEF_EnergyProteinRotamerAndFixedResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                  else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                    pResidueAB->designSiteType == Type_ResidueDesignType_Mutated){
                      EVOEF_EnergyProteinRotamerAndDesignResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                }
              }
            }
          }
          else{//different chains
            if(ChainGetType(pChainA) != Type_Chain_SmallMol){
              for(int b=0; b<ChainGetResidueCount(pChainA); b++){
                Residue* pResidueAB = ChainGetResidue(pChainA, b);
                if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                  EVOEF_EnergyProteinRotamerAndFixedResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
                }
                else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                  pResidueAB->designSiteType == Type_ResidueDesignType_Mutated){
                    EVOEF_EnergyProteinRotamerAndDesignResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
                }
              }
            }
          }
        }
      }
      else if(ChainGetType(pChainI) == Type_Chain_SmallMol){
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){
            total += pRotamerIJ->vdwInternal;
          }
          else{
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                EVOEF_EnergyProteinRotamerAndFixedResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
              }
              else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                pResidueAB->designSiteType == Type_ResidueDesignType_Mutated){
                  EVOEF_EnergyProteinRotamerAndDesignResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
              }
            }
          }
        }
      }

      EnergyTermWeighting(energyTerms);
      total += energyTerms[0];
      partitionPairEnergy[j] = total;
      RotamerExtract(pRotamerIJ);
    }

    for(int j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
      fprintf(outputFile,"%f %d %d %d %d\n",partitionPairEnergy[j],designSiteI,IntArrayGet(pRotamerIndexNewOnSiteI, j),designSiteI,IntArrayGet(pRotamerIndexNewOnSiteI, j));
    }
    free(partitionPairEnergy);
    return Success;
  }

  partitionPairEnergy = (double*)malloc(sizeof(double)*IntArrayGetLength(pRotamerIndexOldOnSiteI)*IntArrayGetLength(pRotamerIndexOldOnSiteK));
  for(int j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
    Rotamer *pRotamerIJ = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteI), IntArrayGet(pRotamerIndexOldOnSiteI, j));
    RotamerRestore(pRotamerIJ,DesignSiteGetRotamers(pDesignSiteI));
    for(int s = 0; s < IntArrayGetLength(pRotamerIndexOldOnSiteK); s++){
      Rotamer *pRotamerKS = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteK), IntArrayGet(pRotamerIndexOldOnSiteK, s));
      RotamerRestore(pRotamerKS,DesignSiteGetRotamers(pDesignSiteK));

      double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
      EnergyTermInitialize(energyTerms);
      if(pDesignSiteI->chainIndex == pDesignSiteK->chainIndex){
        EVOEF_EnergyProteinRotamerAndRotamerSameChain(pRotamerIJ, pRotamerKS, energyTerms);
      }
      else{
        EVOEF_EnergyProteinRotamerAndRotamerDifferentChain(pRotamerIJ, pRotamerKS, energyTerms);
      }
      EnergyTermWeighting(energyTerms);
      double total=0;
      for(int a=1; a<MAX_EVOEF_ENERGY_TERM_NUM;a++){
        total += energyTerms[a];
      }
      partitionPairEnergy[(j*IntArrayGetLength(pRotamerIndexOldOnSiteK)+s)] = total;
      RotamerExtract(pRotamerKS);
    }
    RotamerExtract(pRotamerIJ);
  }

  for(int j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
    for(int s = 0; s < IntArrayGetLength(pRotamerIndexOldOnSiteK); s++){
      fprintf(outputFile,"%f %d %d %d %d\n", 
        partitionPairEnergy[(j*IntArrayGetLength(pRotamerIndexOldOnSiteK)+s)], 
        designSiteI, 
        IntArrayGet(pRotamerIndexNewOnSiteI, j), 
        designSiteK, 
        IntArrayGet(pRotamerIndexNewOnSiteK, j));
    }
  }
  free(partitionPairEnergy);

  return Success;
}

int EnergyMatrixPartitionPairGenerateStoreInMemory(Structure* pStructure, int designSiteI, int designSiteK, 
                    IntArray* pRotamerIndexOldOnSiteI, IntArray* pRotamerIndexOldOnSiteK, 
                    IntArray* pRotamerIndexNewOnSiteI, IntArray* pRotamerIndexNewOnSiteK, 
                    double* energyInMemory, int startIndexThisPartition)
{
  int j, s;
  double    solvNonPolar,vdwAttract,vdwRepul, hbond, elecDesolvation, elecScreenedCoulomb, entropy, total;
  DesignSite *pDesignSiteI = StructureGetDesignSite(pStructure, designSiteI);
  DesignSite *pDesignSiteK = StructureGetDesignSite(pStructure, designSiteK);
  Chain* pChainI = StructureFindChainByName(pStructure, ResidueGetChainName(pDesignSiteI->pResidue));
  Chain* pChainK = StructureFindChainByName(pStructure, ResidueGetChainName(pDesignSiteK->pResidue));
  if(designSiteI > designSiteK || pDesignSiteI == NULL || pDesignSiteK == NULL){
    return ValueError;
  }

  if(designSiteI == designSiteK){
    for(j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
      Rotamer *pRotamerIJ  = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteI), IntArrayGet(pRotamerIndexOldOnSiteI, j));
      RotamerRestore(pRotamerIJ,DesignSiteGetRotamers(pDesignSiteI));

      double total = 0;
      double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
      EnergyTermInitialize(energyTerms);
      Chain *pChainI = StructureGetChain(pStructure, pDesignSiteI->chainIndex);
      if(ChainGetType(pChainI) == Type_Chain_Protein){
        EVOEF_AminoAcidReferenceEnergy(RotamerGetType(pRotamerIJ), energyTerms);
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){//same chain
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(b==pDesignSiteI->resiIndex){
                EVOEF_EnergyProteinRotamerIntraEnergy(pRotamerIJ, energyTerms);
              }
              else{
                if(ChainGetType(pChainA) != Type_Chain_SmallMol){
                  if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                    EVOEF_EnergyProteinRotamerAndFixedResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                  else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                    pResidueAB->designSiteType == Type_ResidueDesignType_Mutated){
                      EVOEF_EnergyProteinRotamerAndDesignResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                }
              }
            }
          }
          else{//different chains
            if(ChainGetType(pChainA) != Type_Chain_SmallMol){
              for(int b=0; b<ChainGetResidueCount(pChainA); b++){
                Residue* pResidueAB = ChainGetResidue(pChainA, b);
                if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                  EVOEF_EnergyProteinRotamerAndFixedResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
                }
                else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                  pResidueAB->designSiteType == Type_ResidueDesignType_Mutated){
                    EVOEF_EnergyProteinRotamerAndDesignResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
                }
              }
            }
          }
        }
      }
      else if(ChainGetType(pChainI) == Type_Chain_SmallMol){
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){
            total += pRotamerIJ->vdwInternal;
          }
          else{
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                EVOEF_EnergyProteinRotamerAndFixedResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
              }
              else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                pResidueAB->designSiteType == Type_ResidueDesignType_Mutated){
                  EVOEF_EnergyProteinRotamerAndDesignResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
              }
            }
          }
        }
      }

      EnergyTermWeighting(energyTerms);
      total += energyTerms[0];
      energyInMemory[startIndexThisPartition + (j*IntArrayGetLength(pRotamerIndexOldOnSiteK)+j)] = total;
      RotamerExtract(pRotamerIJ);
    }

    return Success;
  }

  for(j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
    Rotamer *pRotamerIJ        = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteI), IntArrayGet(pRotamerIndexOldOnSiteI, j));
    RotamerRestore(pRotamerIJ,DesignSiteGetRotamers(pDesignSiteI));
    for(s = 0; s < IntArrayGetLength(pRotamerIndexOldOnSiteK); s++){
      Rotamer *pRotamerKS        = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteK), IntArrayGet(pRotamerIndexOldOnSiteK, s));
      RotamerRestore(pRotamerKS,DesignSiteGetRotamers(pDesignSiteK));

      double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
      EnergyTermInitialize(energyTerms);
      if(pDesignSiteI->chainIndex == pDesignSiteK->chainIndex){
        EVOEF_EnergyProteinRotamerAndRotamerSameChain(pRotamerIJ, pRotamerKS, energyTerms);
      }
      else{
        EVOEF_EnergyProteinRotamerAndRotamerDifferentChain(pRotamerIJ, pRotamerKS, energyTerms);
      }
      EnergyTermWeighting(energyTerms);
      double total=0;
      for(int a=1; a<MAX_EVOEF_ENERGY_TERM_NUM;a++){
        total += energyTerms[a];
      }

      energyInMemory[startIndexThisPartition+(j*IntArrayGetLength(pRotamerIndexOldOnSiteK)+s)] = total;
      RotamerExtract(pRotamerKS);
    }
    RotamerExtract(pRotamerIJ);
  }
  return Success;
}


int EnergyMatrixPartitionPairWriteEnergyToFile(Structure* pStructure, int designSiteI, int designSiteK, 
                    IntArray* pRotamerIndexOldOnSiteI, IntArray* pRotamerIndexOldOnSiteK, 
                    IntArray* pRotamerIndexNewOnSiteI, IntArray* pRotamerIndexNewOnSiteK, 
                    double* energyInMemory, int startIndexThisPartition,
                    FILE* outputFile)
{
  if(designSiteI == designSiteK){
    for(int j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
      double total = energyInMemory[startIndexThisPartition + (j*IntArrayGetLength(pRotamerIndexOldOnSiteK)+j)];
      fprintf(outputFile,"%f %d %d %d %d\n", total, designSiteI, IntArrayGet(pRotamerIndexNewOnSiteI, j), designSiteI, IntArrayGet(pRotamerIndexNewOnSiteI, j));
    }

    return Success;
  }

  for(int j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
    for(int s = 0; s < IntArrayGetLength(pRotamerIndexOldOnSiteK); s++){
      double total = energyInMemory[startIndexThisPartition + (j*IntArrayGetLength(pRotamerIndexOldOnSiteK)+s)];
      fprintf(outputFile,"%f %d %d %d %d\n", total,designSiteI, IntArrayGet(pRotamerIndexNewOnSiteI, j), designSiteK, IntArrayGet(pRotamerIndexNewOnSiteK, s));
    }
  }

  return Success;
}



int SelfEnergyGenerate(Structure *pStructure,char* selfEnergyFilePath){
  FILE* outputFile = NULL;
  int result=Success;
  outputFile = fopen(selfEnergyFilePath,"w");
  if(outputFile==NULL){
    char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sprintf(errMsg,"In StructureGetAllRotamerAreasBaseOnPartition(), Cannot create file %s to write",
      selfEnergyFilePath);
    result = IOError;
    TraceError(errMsg,result);
    return result;        
  }

  for(int designSiteI=0; designSiteI<StructureGetDesignSiteCount(pStructure); designSiteI++){
    DesignSite *pDesignSiteI = StructureGetDesignSite(pStructure, designSiteI);
    RotamerSet *pSetI = DesignSiteGetRotamers(pDesignSiteI);
    Chain *pChainI = StructureGetChain(pStructure, pDesignSiteI->chainIndex);
    for(int j = 0; j < RotamerSetGetCount(pSetI); j++){
      Rotamer *pRotamerIJ  = RotamerSetGet(pSetI, j);
      RotamerRestore(pRotamerIJ,pSetI);
      double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
      EnergyTermInitialize(energyTerms);
      if(ChainGetType(pChainI) == Type_Chain_Protein){
        //reference energy
        EVOEF_AminoAcidReferenceEnergy(RotamerGetType(pRotamerIJ), energyTerms);
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){//same chain
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(b==pDesignSiteI->resiIndex){
                //intra energy
                EVOEF_EnergyProteinRotamerIntraEnergy(pRotamerIJ, energyTerms);
              }
              else{
                if(ChainGetType(pChainA) != Type_Chain_SmallMol){
                  if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                    EVOEF_EnergyProteinRotamerAndFixedResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                  else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                    pResidueAB->designSiteType == Type_ResidueDesignType_Mutated ||
                    pResidueAB->designSiteType == Type_ResidueDesignType_Catalytic){
                      EVOEF_EnergyProteinRotamerAndDesignResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                }
              }
            }
          }
          else{//different chains
            if(ChainGetType(pChainA) != Type_Chain_SmallMol){
              for(int b=0; b<ChainGetResidueCount(pChainA); b++){
                Residue* pResidueAB = ChainGetResidue(pChainA, b);
                if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                  EVOEF_EnergyProteinRotamerAndFixedResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
                }
                else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                  pResidueAB->designSiteType == Type_ResidueDesignType_Mutated ||
                  pResidueAB->designSiteType == Type_ResidueDesignType_Catalytic){
                    EVOEF_EnergyProteinRotamerAndDesignResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
                }
              }
            }
            else{//fixed small molecule
              for(int b=0; b<ChainGetResidueCount(pChainA); b++){
                Residue* pResidueAB = ChainGetResidue(pChainA, b);
                if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                  EVOEF_EnergyRotamerAndFixedLigandResidue(pRotamerIJ,pResidueAB,energyTerms);
                }
                else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                  pResidueAB->designSiteType == Type_ResidueDesignType_Mutated ||
                  pResidueAB->designSiteType == Type_ResidueDesignType_Catalytic){
                    EVOEF_EnergyProteinRotamerAndDesignResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
                }
              }
            }
          }
        }
      }

      else if(ChainGetType(pChainI) == Type_Chain_SmallMol){
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){//same chain
            pRotamerIJ->selfenergy += pRotamerIJ->vdwInternal;
          }
          else{//different chain
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                EVOEF_EnergyLigandRotamerAndFixedResidue(pRotamerIJ,pResidueAB,energyTerms);
              }
              else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                pResidueAB->designSiteType == Type_ResidueDesignType_Mutated ||
                pResidueAB->designSiteType == Type_ResidueDesignType_Catalytic){
                  EVOEF_EnergyLigandRotamerAndDesignResidue(pRotamerIJ,pResidueAB,energyTerms);
              }
            }
          }
        }
      }
      EnergyTermWeighting(energyTerms);
      pRotamerIJ->selfenergy += energyTerms[0];
      fprintf(outputFile, "%d %d %f\n", designSiteI, j, pRotamerIJ->selfenergy);
      RotamerExtract(pRotamerIJ);
    }
  }

  if(FAILED(result)){
    return result;
  }
  fclose(outputFile);

  return Success;
}


int SelfEnergyGenerateWithBBdepEnergy(Structure *pStructure,AAppTable* pAAppTable,RamaTable* pRamaTable,char* selfEnergyFilePath){
  FILE* outputFile = NULL;
  int result=Success;
  outputFile = fopen(selfEnergyFilePath,"w");
  if(outputFile==NULL){
    char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sprintf(errMsg,"In StructureGetAllRotamerAreasBaseOnPartition(), Cannot create file %s to write",
      selfEnergyFilePath);
    result = IOError;
    TraceError(errMsg,result);
    return result;        
  }

  for(int designSiteI=0; designSiteI<StructureGetDesignSiteCount(pStructure); designSiteI++){
    DesignSite *pDesignSiteI = StructureGetDesignSite(pStructure, designSiteI);
    RotamerSet *pSetI = DesignSiteGetRotamers(pDesignSiteI);
    Chain *pChainI = StructureGetChain(pStructure, pDesignSiteI->chainIndex);
    for(int j = 0; j < RotamerSetGetCount(pSetI); j++){
      Rotamer *pRotamerIJ  = RotamerSetGet(pSetI, j);
      RotamerRestore(pRotamerIJ,pSetI);
      double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
      EnergyTermInitialize(energyTerms);
      if(ChainGetType(pChainI) == Type_Chain_Protein){
        //reference energy
        EVOEF_AminoAcidReferenceEnergy(RotamerGetType(pRotamerIJ), energyTerms);
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){//same chain
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(b==pDesignSiteI->resiIndex){
                //physical intra energy
                EVOEF_EnergyProteinRotamerIntraEnergy(pRotamerIJ, energyTerms);
                EVOEF_RotamerPropensityAndRamachandranEnergy(pRotamerIJ,pResidueAB,pAAppTable,pRamaTable,energyTerms);
                EVOEF_RotamerDunbrackEnergy(pRotamerIJ,energyTerms);
              }
              else{
                if(ChainGetType(pChainA) != Type_Chain_SmallMol){
                  if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                    EVOEF_EnergyProteinRotamerAndFixedResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                  else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                    pResidueAB->designSiteType == Type_ResidueDesignType_Mutated ||
                    pResidueAB->designSiteType == Type_ResidueDesignType_Catalytic){
                      EVOEF_EnergyProteinRotamerAndDesignResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                }
              }
            }
          }
          else{//different chains
            if(ChainGetType(pChainA) != Type_Chain_SmallMol){
              for(int b=0; b<ChainGetResidueCount(pChainA); b++){
                Residue* pResidueAB = ChainGetResidue(pChainA, b);
                if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                  EVOEF_EnergyProteinRotamerAndFixedResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
                }
                else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                  pResidueAB->designSiteType == Type_ResidueDesignType_Mutated ||
                  pResidueAB->designSiteType == Type_ResidueDesignType_Catalytic){
                    EVOEF_EnergyProteinRotamerAndDesignResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
                }
              }
            }
            else{//fixed small molecule
              for(int b=0; b<ChainGetResidueCount(pChainA); b++){
                Residue* pResidueAB = ChainGetResidue(pChainA, b);
                if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                  EVOEF_EnergyRotamerAndFixedLigandResidue(pRotamerIJ,pResidueAB,energyTerms);
                }
                else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                  pResidueAB->designSiteType == Type_ResidueDesignType_Mutated ||
                  pResidueAB->designSiteType == Type_ResidueDesignType_Catalytic){
                    EVOEF_EnergyProteinRotamerAndDesignResidueDifferentChain(pRotamerIJ,pResidueAB,energyTerms);
                }
              }
            }
          }
        }
      }

      else if(ChainGetType(pChainI) == Type_Chain_SmallMol){
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){//same chain
            pRotamerIJ->selfenergy += pRotamerIJ->vdwInternal;
          }
          else{//different chain
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(pResidueAB->designSiteType == Type_ResidueDesignType_Fixed){
                EVOEF_EnergyLigandRotamerAndFixedResidue(pRotamerIJ,pResidueAB,energyTerms);
              }
              else if(pResidueAB->designSiteType == Type_ResidueDesignType_Rotameric ||
                pResidueAB->designSiteType == Type_ResidueDesignType_Mutated ||
                pResidueAB->designSiteType == Type_ResidueDesignType_Catalytic){
                  EVOEF_EnergyLigandRotamerAndDesignResidue(pRotamerIJ,pResidueAB,energyTerms);
              }
            }
          }
        }
      }
      EnergyTermWeighting(energyTerms);
      pRotamerIJ->selfenergy += energyTerms[0];
      fprintf(outputFile, "%d %d %f\n", designSiteI, j, pRotamerIJ->selfenergy);
      RotamerExtract(pRotamerIJ);
    }
  }

  if(FAILED(result)){
    return result;
  }
  fclose(outputFile);

  return Success;
}


int SelfEnergyReadAndCheck(Structure* pStructure, RotamerList* pRotamerList, char* selfEnergyFile){
  int i, j, k;
  char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  FILE* pFile = fopen(selfEnergyFile, "r");

  if(pFile == NULL){
    printf("In SelfEnergyReadAndCheck(), cannot open file %s\n", selfEnergyFile);
    return IOError;
  }

  while(fgets(buffer, MAX_LENGTH_ONE_LINE_IN_FILE, pFile)){
    int designSiteI=-1;
    int rotamerIJ=-1;
    double energyIJ=1000.0;
    sscanf(buffer, "%d %d %lf\n",&designSiteI,&rotamerIJ,&energyIJ);
    DesignSite* pSite=StructureGetDesignSite(pStructure,designSiteI);
    RotamerSet* pSet=DesignSiteGetRotamers(pSite);
    Rotamer* pRotamer=RotamerSetGet(pSet,rotamerIJ);
    pRotamer->selfenergy=energyIJ;
  }

  for(i = 0; i < StructureGetDesignSiteCount(pStructure); i++){
    DesignSite* pDesignSite = StructureGetDesignSite(pStructure, i);
    RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
    
    for(j = 0; j < RotamerSetGetCount(pRotamerSet); j++){
      Rotamer* pRotamerJ = RotamerSetGet(pRotamerSet, j);
      if(pRotamerList->remainFlag[i][j] == FALSE) continue;
      for(k = 0; k < RotamerSetGetCount(pRotamerSet); k++){
        Rotamer* pRotamerK = RotamerSetGet(pRotamerSet, k);
        if(pRotamerList->remainFlag[i][k] == FALSE) continue;
        if((pRotamerJ->selfenergy - pRotamerK->selfenergy > SELF_ENERGY_DIFFERENT_THRESHOLD && RotamerAndRotamerInSameType(pRotamerJ,pRotamerK)==FALSE) || 
          (pRotamerJ->selfenergy - pRotamerK->selfenergy > SELF_ENERGY_SAME_THRESHOLD && RotamerAndRotamerInSameType(pRotamerJ,pRotamerK)==TRUE)){
            pRotamerList->remainFlag[i][j] = FALSE;
            break;
        }
      }
    }
  }
  fclose(pFile);

  return Success;
}






//////////////////////////////////////////////
//functions for RotamerList
int RotamerListCreateFromStructure(RotamerList* pThis,Structure* pStructure){
  int i,j;
  pThis->designSiteCount = StructureGetDesignSiteCount(pStructure);
  pThis->rotamerCount = (int*)malloc(sizeof(int)*pThis->designSiteCount);
  pThis->remainFlag = (BOOL**)malloc(sizeof(BOOL*)*pThis->designSiteCount);
  for(i=0;i<pThis->designSiteCount;i++){
    pThis->rotamerCount[i] = RotamerSetGetCount(DesignSiteGetRotamers(StructureGetDesignSite(pStructure,i)));
    pThis->remainFlag[i] = (BOOL*)malloc(sizeof(BOOL)*pThis->rotamerCount[i]);
    for(j=0;j<pThis->rotamerCount[i];j++){
      pThis->remainFlag[i][j] = TRUE;
    }
  }
  pThis->remainRotamerCount = (int*)malloc(sizeof(int)*pThis->designSiteCount);
  for(int i=0;i<pThis->designSiteCount;i++){
    pThis->remainRotamerCount[i]=0;
    for(int j=0;j<pThis->rotamerCount[i];j++){
      if(pThis->remainFlag[i][j]==TRUE){
        pThis->remainRotamerCount[i]++;
      }
    }
  }
  return Success;
}

int RotamerListCreateFromEnergyMatrix(RotamerList* pThis,EnergyMatrix* pEnergyMatrix){
  int i,j;
  pThis->designSiteCount = pEnergyMatrix->designSiteCount;
  pThis->rotamerCount = (int*)malloc(sizeof(int)*EnergyMatrixGetSiteCount(pEnergyMatrix));
  pThis->remainFlag = (BOOL**)malloc(sizeof(BOOL*)*EnergyMatrixGetSiteCount(pEnergyMatrix));
  for(i=0;i<pThis->designSiteCount;i++){
    pThis->rotamerCount[i] = EnergyMatrixGetRotamerCount(pEnergyMatrix,i);
    pThis->remainFlag[i] = (BOOL*)malloc(sizeof(BOOL)*pThis->rotamerCount[i]);
    for(j=0;j<pThis->rotamerCount[i];j++){
      pThis->remainFlag[i][j] = TRUE;
    }
  }
  pThis->rotamerCount = (int*)malloc(sizeof(int)*pThis->designSiteCount);
  for(int i=0;i<pThis->designSiteCount;i++){
    pThis->remainRotamerCount[i]=0;
    for(int j=0;j<pThis->rotamerCount[i];j++){
      if(pThis->remainFlag[i][j]==TRUE){
        pThis->remainRotamerCount[i]++;
      }
    }
  }
  return Success;
}

void RotamerListDestroy(RotamerList* pThis){
  int i;
  for(i=0;i<pThis->designSiteCount;i++){
    free(pThis->remainFlag[i]);
  }
  free(pThis->remainFlag);
  free(pThis->rotamerCount);
  free(pThis->remainRotamerCount);
}

int RotamerListCopy(RotamerList* pThis, RotamerList* pOther)
{
  int i;
  RotamerListDestroy(pThis);
  pThis->designSiteCount=pOther->designSiteCount;
  pThis->rotamerCount=(int*)malloc(pThis->designSiteCount*sizeof(int));
  memcpy(pThis->rotamerCount,pOther->rotamerCount,sizeof(int)*pThis->designSiteCount);
  pThis->remainFlag=(BOOL**)malloc(pThis->designSiteCount*sizeof(BOOL*));
  for(i=0; i<pThis->designSiteCount; i++){
    pThis->remainFlag[i]=(BOOL*)malloc(pThis->rotamerCount[i]*sizeof(BOOL));
    memcpy(pThis->remainFlag[i],pOther->remainFlag[i],sizeof(BOOL)*pThis->rotamerCount[i]);
  }
  pThis->remainRotamerCount=(int*)malloc(pThis->designSiteCount*sizeof(int));
  memcpy(pThis->remainRotamerCount,pOther->remainRotamerCount,sizeof(int)*pThis->designSiteCount);
  return Success;
}

int RotamerListRead(RotamerList* pThis,char* filepath){
  int i,j;
  FILE* pFile;
  int result;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  char flag[MAX_LENGTH_ONE_LINE_IN_FILE+1];

  pFile = fopen(filepath,"r");
  if(pFile==NULL){
    result = IOError;
    sprintf(errMsg,"In RemainRotamerListRead(), cannot open file %s to read",filepath);
    TraceError(errMsg,result);
    return result;
  }

  while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pFile)){
    sscanf(line,"%d %d %s",&i,&j,flag);
    if(i>=pThis->designSiteCount){
      sprintf(errMsg,"In RemainRotamerListRead(), "
        "Structure does not have designSite #%d. At this line:\n%s",i,line);
      result = AssertionError;
      TraceError(errMsg,result);
      return result;
    }
    if(j>=pThis->rotamerCount[i]){
      sprintf(errMsg,"In RemainRotamerListRead(), "
        "Structure does not have rotamer #%d on designSite #%d. At this line:\n%s",j,i,line);
      result = AssertionError;
      TraceError(errMsg,result);
      return result;
    }
    if(strcmp(flag,"TRUE")==0){
      pThis->remainFlag[i][j] = TRUE;
    }
    else if(strcmp(flag,"FALSE")==0){
      pThis->remainFlag[i][j] = FALSE;
    }
    else{
      sprintf(errMsg,"In RemainRotamerListRead(), "
        "File format error. At this line:\n%s",line);
      result = FormatError;
      TraceError(errMsg,result);
      return result;
    }
  }

  for(int i=0;i<pThis->designSiteCount;i++){
    pThis->remainRotamerCount[i]=0;
    for(int j=0;j<pThis->rotamerCount[i];j++){
      if(pThis->remainFlag[i][j]==TRUE){
        pThis->remainRotamerCount[i]++;
      }
    }
  }
  return Success;
}

int RotamerListWrite(RotamerList* pThis,char* filepath){
  int i,j;
  int result;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  FILE* pFile = fopen(filepath,"w");
  if(pFile==NULL){
    result = IOError;
    sprintf(errMsg,"In RemainRotamerListWrite(), cannot open file %s to read",filepath);
    TraceError(errMsg,result);
    return result;
  }
  for(i=0;i<pThis->designSiteCount;i++){
    for(j=0;j<pThis->rotamerCount[i];j++){
      fprintf(pFile,"%d\t%d\t%s\n",i,j,  
        pThis->remainFlag[i][j]? "TRUE" : "FALSE"  );
    }
  }
  fclose(pFile);
  return Success;
}



int RotamerListShow(RotamerList* pThis){
  int i,j;
  int totalRotamerCount=0;
  for(i=0;i<pThis->designSiteCount;i++){
    int remainRotamerCount=0;
    printf("DesignSite %d, %d Rotamers:\n",i,pThis->rotamerCount[i]);
    for(j=0;j<pThis->rotamerCount[i];j++){
      if(pThis->remainFlag[i][j]==TRUE){
        printf("Rotamer %d\n",j);
        remainRotamerCount++;
      }
    }
    totalRotamerCount+=remainRotamerCount;
    printf("%d rotamers left\n",remainRotamerCount);
  }
  printf("total rotamer count:\t%d\n",totalRotamerCount);
  return Success;
}


int RotamerListAndEnergyMatrixDelete(RotamerList* pRotamerList, 
                                     EnergyMatrix* pEnergyMatrix, 
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
    sprintf(errMsg,"In RotamerListAndEnergyMatrixDelete(), "
      "parameter #3 must be a int array with odd number length");
    TraceError(errMsg,result);
    return result;
  }
  for(i=0;i<IntArrayGetLength(pIndexOfDeletedRotamers);i+=2){
    int siteIndex = IntArrayGet(pIndexOfDeletedRotamers,i);
    int rotamerIndex = IntArrayGet(pIndexOfDeletedRotamers,i+1);
    if( siteIndex >= pRotamerList->designSiteCount ||
      rotamerIndex >= EnergyMatrixGetBlock(pEnergyMatrix,siteIndex,siteIndex)->RotamerCountSiteI ){
        result = AssertionError;
        sprintf(errMsg,"In RotamerListAndEnergyMatrixDelete(), rotamer #%d on site #%d does not exist",
          rotamerIndex,siteIndex);
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
      sprintf(errMsg,"In RotamerListAndEnergyMatrixDelete(), "
        "all rotamers on design site %d have been deleted",i);
      TraceError(errMsg,result);
      return result;
    }
  }

  //Update the EnergyMatrix
  for(i=0;i<pRotamerList->designSiteCount;i++){
    for(j=i;j<pRotamerList->designSiteCount;j++){
      EnergyMatrixBlockUpdate(EnergyMatrixGetBlock(pEnergyMatrix,i,j),
        rotamerDeletedFlag[i],rotamerDeletedFlag[j]);
    }
  }

  for(i=0;i<pRotamerList->designSiteCount;i++){
    free(rotamerDeletedFlag[i]);
  }
  free(rotamerDeletedFlag);
  return Success;
}



int EnergyMatrixBlockUpdate(EnergyMatrixBlock* pEnergyBlock,BOOL* siteIRotamerDeletedFlag,BOOL* siteKRotamerDeletedFlag){
  int j,s;
  int newRotamerCountSiteI = 0;
  int newRotamerCountSiteK = 0;
  double* newEnergyIK;
  double* pPositionInNewEnergyIK;

  for(j=0;j<pEnergyBlock->RotamerCountSiteI;j++){
    if(siteIRotamerDeletedFlag[j]==FALSE){
      newRotamerCountSiteI++;
    }
  }
  for(s=0;s<pEnergyBlock->RotamerCountSiteK;s++){
    if(siteKRotamerDeletedFlag[s]==FALSE){
      newRotamerCountSiteK++;
    }
  }
  newEnergyIK = (double*)malloc(sizeof(double)*newRotamerCountSiteI*newRotamerCountSiteK);
  pPositionInNewEnergyIK = &newEnergyIK[0];
  for(j=0;j<pEnergyBlock->RotamerCountSiteI;j++){
    if(siteIRotamerDeletedFlag[j] == TRUE){
      continue;
    }
    for(s=0;s<pEnergyBlock->RotamerCountSiteK;s++){
      if(siteKRotamerDeletedFlag[s] == TRUE){
        continue;
      }
      *pPositionInNewEnergyIK = *EnergyMatrixBlockGet(pEnergyBlock,j,s);
      pPositionInNewEnergyIK++;
    }
  }

  pEnergyBlock->RotamerCountSiteI = newRotamerCountSiteI;
  pEnergyBlock->RotamerCountSiteK = newRotamerCountSiteK;
  free(pEnergyBlock->energyIK);
  pEnergyBlock->energyIK = newEnergyIK;

  return Success;
}




//////////////////////////////////////////////////////////////////////////
//functions used for reducing energy matrix size and energy calculation
//////////////////////////////////////////////////////////////////////////



int EnergyMatrixUpdateByRotamerList(EnergyMatrix* pOriginalMatrix, RotamerList* pList)
{
  int i, j;
  IntArray  deleteRotamerIndex;
  RotamerList  rotamerList;

  IntArrayCreate(&deleteRotamerIndex, 0);
  RotamerListCreateFromEnergyMatrix(&rotamerList, pOriginalMatrix);
  for(i=0; i<pList->designSiteCount; i++){
    for(j=0; j<pList->rotamerCount[i]; j++){
      if(pList->remainFlag[i][j]==FALSE){
        IntArrayAppend(&deleteRotamerIndex,i);
        IntArrayAppend(&deleteRotamerIndex,j);
      }
    }
  }

  RotamerListAndEnergyMatrixDelete(&rotamerList, pOriginalMatrix, &deleteRotamerIndex);
  RotamerListDestroy(&rotamerList);
  IntArrayDestroy(&deleteRotamerIndex);
  return Success;
}


int StructureShowDesignSitesAfterRotamerDelete(Structure* pStructure, RotamerList* pList,FILE* pFile){
  int i, j, totalRotamerCount;
  double conformationSpace=0.0;

  totalRotamerCount = 0;
  if(StructureGetDesignSiteCount(pStructure) != pList->designSiteCount){
    printf("In StructureShowDesignSitesAfterRotamerDelete(), the sizes of Structure and RotamerList are different, please check!\n");
    printf("Structure design site count: %d, RotamerList design site count: %d\n", StructureGetDesignSiteCount(pStructure), pList->designSiteCount);
    return ValueError;
  }
  for(i = 0; i < StructureGetDesignSiteCount(pStructure); i++){
    DesignSite* pSiteI = StructureGetDesignSite(pStructure, i);
    RotamerSet* pSetI = DesignSiteGetRotamers(pSiteI);
    int remainRotamerCount = 0;

    if(RotamerSetGetCount(pSetI) != pList->rotamerCount[i]){
      printf("In StructureShowDesignSitesAfterRotamerDelete(), the rotamer numbers on design site %d of Structure and RotamerList are different, please check!\n", i);
      printf("rotamer number on Structure: %d, rotamer number on RotamerList: %d\n", RotamerSetGetCount(pSetI), pList->rotamerCount[i]);
      return ValueError;
    }
    for(j = 0; j < RotamerSetGetCount(pSetI); j++){
      Rotamer* pRotamer = RotamerSetGet(pSetI, j);
      if(pList->remainFlag[i][j] == TRUE){
        remainRotamerCount++;
      }
    }
    fprintf(stdout,"design site %3d : %3s %s %4d, %6d rotamers.  ",
      i,ResidueGetName(pSiteI->pResidue),ResidueGetChainName(pSiteI->pResidue),ResidueGetPosInChain(pSiteI->pResidue),remainRotamerCount);
    totalRotamerCount += remainRotamerCount;
    if(remainRotamerCount>0){
      conformationSpace += log((double)remainRotamerCount)/log(10.0);
    }
    switch(pSiteI->pResidue->designSiteType){
      case Type_ResidueDesignType_Catalytic:
        fprintf(pFile,"catalytic\n"); break;
      case Type_ResidueDesignType_Fixed:
        fprintf(pFile,"fixed\n"); break;
      case Type_ResidueDesignType_Mutated:
        fprintf(pFile,"mutated\n"); break;
      case Type_ResidueDesignType_SmallMol:
        fprintf(pFile,"smallmol\n"); break;
      case Type_ResidueDesignType_Rotameric:
        fprintf(pFile,"rotameric\n"); break;
      default:
        break;
    }
  }
  fprintf(stdout, "total rotamer count: %d, conformation space: 1e^%.1f\n", totalRotamerCount, conformationSpace);

  return Success;
}




int RotamerDeleteBySelfEnergyCheckNew(EnergyMatrix *pMatrix, Structure *pStructure, RotamerList* pList, IntArray *pDeleteList, EnergyMatrix *pRemainFlag){
  int i,j,u;

  IntArrayResize(pDeleteList, 0);
  for(i=0; i<pMatrix->designSiteCount; i++){
    EnergyMatrixBlock *pBlock = EnergyMatrixGetBlock(pMatrix, i, i);
    for(j=0; j<pBlock->RotamerCountSiteI; j++){
      Rotamer *pRotamerIJ;
      int trueIndexJ;
      RotamerOriginalIndexGet(pList, i, j, &trueIndexJ);
      pRotamerIJ = RotamerSetGet(DesignSiteGetRotamers(StructureGetDesignSite(pStructure, i)), trueIndexJ);
      for(u=0; u<pBlock->RotamerCountSiteI; u++){
        Rotamer *pRotamerIU;
        int trueIndexU;
        RotamerOriginalIndexGet(pList, i, u, &trueIndexU);
        pRotamerIU = RotamerSetGet(DesignSiteGetRotamers(StructureGetDesignSite(pStructure, i)), trueIndexU);
        if(*EnergyMatrixGet(pRemainFlag, i, i, u, u) < 0.0) continue;
        //for different rotamer types
        if(*EnergyMatrixGet(pMatrix, i, i, j, j) > *EnergyMatrixGet(pMatrix, i, i, u, u) + SELF_ENERGY_DIFFERENT_THRESHOLD){
          IntArrayAppend(pDeleteList, i);
          IntArrayAppend(pDeleteList, j);
          *EnergyMatrixGet(pRemainFlag, i, i, j, j) = -1.0;
          break;
        }
        //for same rotamer types
        else if(RotamerAndRotamerInSameType(pRotamerIJ,pRotamerIU)==TRUE &&
          *EnergyMatrixGet(pMatrix, i, i, j, j) > *EnergyMatrixGet(pMatrix, i, i, u, u) + SELF_ENERGY_SAME_THRESHOLD){
            IntArrayAppend(pDeleteList, i);
            IntArrayAppend(pDeleteList, j);
            *EnergyMatrixGet(pRemainFlag, i, i, j, j) = -1.0;
            break;
        }
      }
    }
  }

  return 0;
}


int DesignSiteGetRemainRotamerCount(int designSiteI, RotamerList* pList)
{
  int j;
  int remainRotamerCount = 0;

  if(designSiteI < 0 || designSiteI >= pList->designSiteCount){
    printf("In DesignSiteGetRemainRotamerCount(), Invalid parameter!\n");
    return ValueError;
  }

  for(j = 0; j < pList->rotamerCount[designSiteI]; j++){
    if(pList->remainFlag[designSiteI][j] == TRUE){
      remainRotamerCount++;
    }
  }

  return remainRotamerCount;
}


int RotamerListUpdateByDeleteArray(RotamerList* pList, IntArray* pDeleteArray)
{
  int i;
  for(i = 0; i < IntArrayGetLength(pDeleteArray); i += 2){
    int designSiteI = IntArrayGet(pDeleteArray, i);
    int rotamerIJ = IntArrayGet(pDeleteArray, i+1);
    if(designSiteI < 0 || designSiteI >= pList->designSiteCount|| 
      rotamerIJ < 0 || rotamerIJ >= pList->rotamerCount[designSiteI]){
        printf("In RotamerListUpdateByDeleteArray(), invalid parameter!\n");
        return ValueError;
    }

    pList->remainFlag[designSiteI][rotamerIJ] = FALSE;
  }

  return Success;
}


int RotamerOriginalIndexGet(RotamerList* pList, int designSiteI, int rotamerIJ, int *trueIndexIJ)
{
  int jj, flag_j = 0;

  for(jj = 0; jj < pList->rotamerCount[designSiteI]; jj++){
    if(pList->remainFlag[designSiteI][jj] == FALSE) continue;
    if(flag_j == rotamerIJ) break;
    flag_j++;
  }
  if(jj == pList->rotamerCount[designSiteI]){
    return DataNotExistError;
  }
  else{
    *trueIndexIJ = jj;
    return Success;
  }
}

int RotamerReducedIndexGet(RotamerList* pList, int designSiteI, int trueIndexIJ, int *reducedIndex)
{
  int jj, flag_j = 0;
  for(jj = 0; jj < pList->rotamerCount[designSiteI]; jj++){
    if(pList->remainFlag[designSiteI][jj] == FALSE) continue;
    if(jj == trueIndexIJ) break;
    flag_j++;
  }
  if(jj == pList->rotamerCount[designSiteI]){
    return DataNotExistError;
  }
  else{
    *reducedIndex = flag_j;
    return Success;
  }
}


int EnergyMatrixReadNew(EnergyMatrix* pThis, RotamerList* pList, char* energyMatrixFile){
  int i, j, k;
  int *remainRotamerCount = (int*)malloc(sizeof(int)*pList->designSiteCount);

  for(i = 0; i < pList->designSiteCount; i++){
    remainRotamerCount[i] = 0;
    for(j = 0; j < pList->rotamerCount[i]; j++){
      if(pList->remainFlag[i][j] == TRUE){
        remainRotamerCount[i]++;
      }
    }
  }

  pThis->designSiteCount = pList->designSiteCount;
  pThis->blocks = (EnergyMatrixBlock*)malloc(sizeof(EnergyMatrixBlock)*pThis->designSiteCount*pThis->designSiteCount);
  for(i = 0; i < pThis->designSiteCount*pThis->designSiteCount; i++){
    EnergyMatrixBlockCreate(&pThis->blocks[i]);
  }

  for(i = 0; i < pThis->designSiteCount; i++){
    for(k = i; k < pThis->designSiteCount; k++){
      EnergyMatrixBlock* pBlockIK = EnergyMatrixGetBlock(pThis, i, k);
      pBlockIK->DesignSiteI = i;
      pBlockIK->DesignSiteK = k;
      pBlockIK->RotamerCountSiteI = remainRotamerCount[i];
      pBlockIK->RotamerCountSiteK = remainRotamerCount[k];
      pBlockIK->energyIK = (double*)malloc(sizeof(double)*pBlockIK->RotamerCountSiteI*pBlockIK->RotamerCountSiteK);
      for(j = 0; j < pBlockIK->RotamerCountSiteI*pBlockIK->RotamerCountSiteK; j++){
        pBlockIK->energyIK[j] = 1000.0;
      }
    }
  }

  FILE* pFile = fopen(energyMatrixFile, "r");
  if(pFile == NULL){
    printf("In EnergyMatrixReaderNew(), cannot open file %s", energyMatrixFile);
    return IOError;
  }

  long int lineCounter = 0;
  char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(fgets(buffer, MAX_LENGTH_ONE_LINE_IN_FILE, pFile)){
    double energy;
    int siteI, siteK, rotamerIJ, rotamerKS;

    sscanf(buffer, "%lf %d %d %d %d\n", &energy, &siteI, &rotamerIJ, &siteK, &rotamerKS);
    *EnergyMatrixGet(pThis, siteI, siteK, rotamerIJ, rotamerKS) = energy;
  }
  fclose(pFile);

  free(remainRotamerCount);

  return Success;
}

int DeleteArrayGenerateFromTwoRotamerList(IntArray* pDeleteRotamerArray, RotamerList* pOld, RotamerList* pNew)
{
  int i, j;
  
  IntArrayResize(pDeleteRotamerArray, 0);
  for(i = 0; i < pOld->designSiteCount; i++){
    int index = 0;
    for(j = 0; j < pOld->rotamerCount[i]; j++){
      if(pOld->remainFlag[i][j] == TRUE && pNew->remainFlag[i][j] == TRUE) index++;
      if(pOld->remainFlag[i][j] == TRUE && pNew->remainFlag[i][j] == FALSE){
        IntArrayAppend(pDeleteRotamerArray, i);
        IntArrayAppend(pDeleteRotamerArray, index);
        index++;
      }
    }
  }

  return Success;
}


int RotamerListCreateFromFile(RotamerList *pThis, char* filepath)
{
  int i,j;
  FILE* pFile;
  int result;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  char flag[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  int max_design_site_count = 0;
  IntArray max_remain_rotamer_counts;

  pFile = fopen(filepath,"r");
  if(pFile==NULL){
    printf("In RotamerListCreateFromFile(), cannot open file %s to read\n", filepath);
    return IOError;
  }

  while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pFile)){
    if(sscanf(line,"%d %d %s",&i,&j,flag) != 3){
      printf("In RotamerListCreateFromFile(), cannot open file %s to read\n", filepath);
      return IOError;
    }
    if(i >= max_design_site_count) max_design_site_count++;
  }
  fseek(pFile, 0, SEEK_SET);

  IntArrayCreate(&max_remain_rotamer_counts, max_design_site_count);
  for(i = 0; i < IntArrayGetLength(&max_remain_rotamer_counts); i++){
    IntArraySet(&max_remain_rotamer_counts, i, 0);
  }

  while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pFile)){
    sscanf(line,"%d %d %s",&i,&j,flag);
    if(j >= IntArrayGet(&max_remain_rotamer_counts, i)){
      IntArraySet(&max_remain_rotamer_counts, i, j+1);
    }
  }
  fseek(pFile, 0, SEEK_SET);

  pThis->designSiteCount = max_design_site_count;
  pThis->rotamerCount = (int*)malloc(sizeof(int)*pThis->designSiteCount);
  pThis->remainFlag = (BOOL**)malloc(sizeof(BOOL*)*pThis->designSiteCount);
  for(i = 0; i < pThis->designSiteCount; i++){
    pThis->rotamerCount[i] = IntArrayGet(&max_remain_rotamer_counts, i);
    pThis->remainFlag[i] = (BOOL*)malloc(sizeof(BOOL)*pThis->rotamerCount[i]);
    for(j = 0;j < pThis->rotamerCount[i]; j++){
      pThis->remainFlag[i][j] = TRUE;
    }
  }
  while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pFile)){
    sscanf(line,"%d %d %s",&i,&j,flag);
    if(strcmp(flag, "TRUE") == 0){
      pThis->remainFlag[i][j] = TRUE;
    }
    else if(strcmp(flag, "FALSE") == 0){
      pThis->remainFlag[i][j] = FALSE;
    }
  }

  fclose(pFile);
  IntArrayDestroy(&max_remain_rotamer_counts);

  return Success;
}

