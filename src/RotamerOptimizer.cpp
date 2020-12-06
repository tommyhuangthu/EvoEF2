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

#include "RotamerOptimizer.h"
#include "EnergyFunction.h"
#include <string.h>

int ProteinSiteOptimizeRotamer(Structure *pStructure, int chainIndex, int resiIndex){
  DesignSite *pDesignSite = StructureFindDesignSite(pStructure,chainIndex,resiIndex);
  if(pDesignSite==NULL) return Success;
  RotamerSet *pRotSet = DesignSiteGetRotamers(pDesignSite);
  Residue *pDesign = pDesignSite->pResidue;
  double *energyArrayOfRotamers = (double *)malloc(sizeof(double)*RotamerSetGetCount(pRotSet));
  for(int i=0; i<RotamerSetGetCount(pRotSet); ++i) energyArrayOfRotamers[i]=0.0;

  //step 1: find out residues within 5 angstroms to the design site of interest;
  int surroundingResiNum = 0;
  Residue **ppSurroundingResidues = NULL;
  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    for(int j = 0; j < pChainI->residueNum; j++){
      Residue *pResi2 = pChainI->residues + j;
      if(strcmp(pDesign->name, pResi2->name) == 0 && pDesign->posInChain == pResi2->posInChain){
        continue;
      }
      else if(AtomArrayCalcMinDistance(&pDesign->atoms, &pResi2->atoms)<ENERGY_DISTANCE_CUTOFF){
        surroundingResiNum++;
        ppSurroundingResidues = (Residue **)realloc(ppSurroundingResidues, sizeof(Residue*)*surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum-1] = pResi2;
      }
    }
  }

  // step 2: calculate the energy between the rotamers of the design site
  double minEnergy = 1000.0;
  int minEnergyRotIndex = -1;
  int minEnergyRotTypeIndex[22];
  double minEnergyRotTypeEnergy[22];
  StringArray rotTypes;
  StringArrayCreate(&rotTypes);
  for(int ir = 0; ir < RotamerSetGetCount(pRotSet); ir++){
    double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
    EnergyTermInitialize(energyTerms);
    Rotamer *pRotIR = RotamerSetGet(pRotSet, ir);
    RotamerRestore(pRotIR, pRotSet);

    Residue tempResidue;
    ResidueCreate(&tempResidue);
    ResidueCopy(&tempResidue, pDesign);
    ResidueSetName(&tempResidue, pRotIR->type);
    AtomArrayCopy(&tempResidue.atoms, &pRotIR->atoms);
    BondSetCopy(&tempResidue.bonds, &pRotIR->bonds);
    EVOEF_AminoAcidReferenceEnergy(tempResidue.name,energyTerms);
    EVOEF_EnergyResidueIntraEnergy(&tempResidue,energyTerms);
    for(int is = 0; is < surroundingResiNum; is++){
      Residue *pResIS = ppSurroundingResidues[is];
      if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain-1){
        EVOEF_EnergyResidueAndNextResidue(&tempResidue,pResIS,energyTerms);
      }
      else if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain+1){
        EVOEF_EnergyResidueAndNextResidue(pResIS,&tempResidue,energyTerms);
      }
      else{
        if(strcmp(ResidueGetChainName(&tempResidue),ResidueGetChainName(pResIS))==0){
          EVOEF_EnergyResidueAndOtherResidueSameChain(&tempResidue,pResIS,energyTerms);
        }
        else{
          EVOEF_EnergyResidueAndOtherResidueDifferentChain(&tempResidue,pResIS,energyTerms);
        }
      }
    }

    if(FALSE){
      //for debug: energy terms, not weighted
      printf("ROTAMER %4d: %1s %4d %3s, ENERGY: %6.2f\n", ir, RotamerGetChainName(pRotIR), RotamerGetPosInChain(pRotIR), RotamerGetType(pRotIR), energyArrayOfRotamers[ir]);
      printf("reference_ALA         =            %8.2f\n", energyTerms[ 1]);  
      printf("reference_CYS         =            %8.2f\n", energyTerms[ 2]);
      printf("reference_ASP         =            %8.2f\n", energyTerms[ 3]);
      printf("reference_GLU         =            %8.2f\n", energyTerms[ 4]);
      printf("reference_PHE         =            %8.2f\n", energyTerms[ 5]);
      printf("reference_GLY         =            %8.2f\n", energyTerms[ 6]);
      printf("reference_HIS         =            %8.2f\n", energyTerms[ 7]);
      printf("reference_ILE         =            %8.2f\n", energyTerms[ 8]);
      printf("reference_LYS         =            %8.2f\n", energyTerms[ 9]);
      printf("reference_LEU         =            %8.2f\n", energyTerms[10]);
      printf("reference_MET         =            %8.2f\n", energyTerms[11]);
      printf("reference_ASN         =            %8.2f\n", energyTerms[12]);
      printf("reference_PRO         =            %8.2f\n", energyTerms[13]);
      printf("reference_GLN         =            %8.2f\n", energyTerms[14]);
      printf("reference_ARG         =            %8.2f\n", energyTerms[15]);
      printf("reference_SER         =            %8.2f\n", energyTerms[16]);
      printf("reference_THR         =            %8.2f\n", energyTerms[17]);
      printf("reference_VAL         =            %8.2f\n", energyTerms[18]);
      printf("reference_TRP         =            %8.2f\n", energyTerms[19]);
      printf("reference_TYR         =            %8.2f\n", energyTerms[20]);

      printf("intraR_vdwatt         =            %8.2f\n", energyTerms[21]);
      printf("intraR_vdwrep         =            %8.2f\n", energyTerms[22]);
      printf("intraR_electr         =            %8.2f\n", energyTerms[23]);
      printf("intraR_deslvP         =            %8.2f\n", energyTerms[24]);
      printf("intraR_deslvH         =            %8.2f\n", energyTerms[25]);
      printf("intraR_hbscbb_dis     =            %8.2f\n", energyTerms[26]);
      printf("intraR_hbscbb_the     =            %8.2f\n", energyTerms[27]);
      printf("intraR_hbscbb_phi     =            %8.2f\n", energyTerms[28]);
      printf("aapropensity          =            %8.2f\n", energyTerms[91]);
      printf("ramachandran          =            %8.2f\n", energyTerms[92]);
      printf("dunbrack              =            %8.2f\n", energyTerms[93]);

      printf("interS_vdwatt         =            %8.2f\n", energyTerms[31]);
      printf("interS_vdwrep         =            %8.2f\n", energyTerms[32]);
      printf("interS_electr         =            %8.2f\n", energyTerms[33]);
      printf("interS_deslvP         =            %8.2f\n", energyTerms[34]);
      printf("interS_deslvH         =            %8.2f\n", energyTerms[35]);
      printf("interS_ssbond         =            %8.2f\n", energyTerms[36]);
      printf("interS_hbbbbb_dis     =            %8.2f\n", energyTerms[41]);
      printf("interS_hbbbbb_the     =            %8.2f\n", energyTerms[42]);
      printf("interS_hbbbbb_phi     =            %8.2f\n", energyTerms[43]);
      printf("interS_hbscbb_dis     =            %8.2f\n", energyTerms[44]);
      printf("interS_hbscbb_the     =            %8.2f\n", energyTerms[45]);
      printf("interS_hbscbb_phi     =            %8.2f\n", energyTerms[46]);
      printf("interS_hbscsc_dis     =            %8.2f\n", energyTerms[47]);
      printf("interS_hbscsc_the     =            %8.2f\n", energyTerms[48]);
      printf("interS_hbscsc_phi     =            %8.2f\n", energyTerms[49]);

      printf("interD_vdwatt         =            %8.2f\n", energyTerms[51]);
      printf("interD_vdwrep         =            %8.2f\n", energyTerms[52]);
      printf("interD_electr         =            %8.2f\n", energyTerms[53]);
      printf("interD_deslvP         =            %8.2f\n", energyTerms[54]);
      printf("interD_deslvH         =            %8.2f\n", energyTerms[55]);
      printf("interD_ssbond         =            %8.2f\n", energyTerms[56]);
      printf("interD_hbbbbb_dis     =            %8.2f\n", energyTerms[61]);
      printf("interD_hbbbbb_the     =            %8.2f\n", energyTerms[62]);
      printf("interD_hbbbbb_phi     =            %8.2f\n", energyTerms[63]);
      printf("interD_hbscbb_dis     =            %8.2f\n", energyTerms[64]);
      printf("interD_hbscbb_the     =            %8.2f\n", energyTerms[65]);
      printf("interD_hbscbb_phi     =            %8.2f\n", energyTerms[66]);
      printf("interD_hbscsc_dis     =            %8.2f\n", energyTerms[67]);
      printf("interD_hbscsc_the     =            %8.2f\n", energyTerms[68]);
      printf("interD_hbscsc_phi     =            %8.2f\n", energyTerms[69]);      
    }

    //total energy: weighted
    EnergyTermWeighting(energyTerms);
    energyArrayOfRotamers[ir] = energyTerms[0];

    if(energyArrayOfRotamers[ir] < minEnergy){
      minEnergy = energyArrayOfRotamers[ir];
      minEnergyRotIndex = ir;
      ResidueCopy(pDesign, &tempResidue);
    }

    if(ir == 0){
      StringArrayAppend(&rotTypes, tempResidue.name);
      minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
      minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
    }
    else{
      int pos = -1;
      StringArrayFind(&rotTypes, tempResidue.name, &pos);
      if(pos == -1){
        StringArrayAppend(&rotTypes, tempResidue.name);
        minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
        minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
      }
      else{
        if(energyArrayOfRotamers[ir] < minEnergyRotTypeEnergy[pos]){
          minEnergyRotTypeIndex[pos] = ir;
          minEnergyRotTypeEnergy[pos] = energyArrayOfRotamers[ir];
        }
      }
    }

    RotamerExtract(pRotIR);
    ResidueDestroy(&tempResidue);
  }

  if(energyArrayOfRotamers != NULL){
    free(energyArrayOfRotamers);
    energyArrayOfRotamers = NULL;
  }
  if(ppSurroundingResidues != NULL){
    free(ppSurroundingResidues);
    ppSurroundingResidues = NULL;
  }
  StringArrayDestroy(&rotTypes);

  return Success;
}


int ProteinSiteOptimizeRotamerLocally(Structure *pStructure, int chainIndex, int resiIndex, double rmsdcutoff){
  DesignSite *pDesignSite = StructureFindDesignSite(pStructure,chainIndex,resiIndex);
  if(pDesignSite==NULL) return Success;
  RotamerSet *pRotSet = DesignSiteGetRotamers(pDesignSite);
  Residue *pDesign = pDesignSite->pResidue;
  double *energyArrayOfRotamers = (double *)malloc(sizeof(double)*RotamerSetGetCount(pRotSet));
  for(int i=0; i<RotamerSetGetCount(pRotSet); ++i) energyArrayOfRotamers[i]=0.0;

  //step 1: find out residues within 6 angstroms to the design site of interest;
  int surroundingResiNum = 0;
  Residue **ppSurroundingResidues = NULL;
  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    for(int j = 0; j < pChainI->residueNum; j++){
      Residue *pResi2 = pChainI->residues + j;
      if(strcmp(pDesign->name, pResi2->name) == 0 && pDesign->posInChain == pResi2->posInChain){
        continue;
      }
      else if(AtomArrayCalcMinDistance(&pDesign->atoms, &pResi2->atoms)<ENERGY_DISTANCE_CUTOFF){
        surroundingResiNum++;
        ppSurroundingResidues = (Residue **)realloc(ppSurroundingResidues, sizeof(Residue*)*surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum-1] = pResi2;
      }
    }
  }

  Residue original;
  ResidueCreate(&original);
  ResidueCopy(&original,pDesign);

  // step 2: calculate the energy between the rotamers of the design site
  double minEnergy = 1000.0;
  int minEnergyRotIndex = -1;
  int minEnergyRotTypeIndex[22];
  double minEnergyRotTypeEnergy[22];
  StringArray rotTypes;
  StringArrayCreate(&rotTypes);
  for(int ir = 0; ir < RotamerSetGetCount(pRotSet); ir++){
    double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
    EnergyTermInitialize(energyTerms);
    Rotamer *pRotIR = RotamerSetGet(pRotSet, ir);
    RotamerRestore(pRotIR, pRotSet);

    Residue tempResidue;
    ResidueCreate(&tempResidue);
    ResidueCopy(&tempResidue, pDesign);
    ResidueSetName(&tempResidue, pRotIR->type);
    AtomArrayCopy(&tempResidue.atoms, &pRotIR->atoms);
    BondSetCopy(&tempResidue.bonds, &pRotIR->bonds);
    //skip the rotamers that has large side-chain shift
    double realrmsd=ResidueAndResidueSidechainRMSD(&tempResidue,&original);
    if(realrmsd>rmsdcutoff && realrmsd<1e4){
      ResidueDestroy(&tempResidue);
      continue;
    }
    EVOEF_AminoAcidReferenceEnergy(tempResidue.name,energyTerms);
    EVOEF_EnergyResidueIntraEnergy(&tempResidue,energyTerms);
    for(int is = 0; is < surroundingResiNum; is++){
      Residue *pResIS = ppSurroundingResidues[is];
      if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain-1){
        EVOEF_EnergyResidueAndNextResidue(&tempResidue,pResIS,energyTerms);
      }
      else if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain+1){
        EVOEF_EnergyResidueAndNextResidue(pResIS,&tempResidue,energyTerms);
      }
      else{
        if(strcmp(tempResidue.chainName, pResIS->chainName)==0){
          EVOEF_EnergyResidueAndOtherResidueSameChain(&tempResidue,pResIS,energyTerms);
        }
        else{
          EVOEF_EnergyResidueAndOtherResidueDifferentChain(&tempResidue,pResIS,energyTerms);
        }
      }
    }

    if(FALSE){
      printf("ROTAMER %4d: %1s %4d %3s, ENERGY: %6.2f\n", ir, RotamerGetChainName(pRotIR), RotamerGetPosInChain(pRotIR), RotamerGetType(pRotIR), energyArrayOfRotamers[ir]);
      printf("reference_ALA         =            %8.2f\n", energyTerms[ 1]);
      printf("reference_CYS         =            %8.2f\n", energyTerms[ 2]);
      printf("reference_ASP         =            %8.2f\n", energyTerms[ 3]);
      printf("reference_GLU         =            %8.2f\n", energyTerms[ 4]);
      printf("reference_PHE         =            %8.2f\n", energyTerms[ 5]);
      printf("reference_GLY         =            %8.2f\n", energyTerms[ 6]);
      printf("reference_HIS         =            %8.2f\n", energyTerms[ 7]);
      printf("reference_ILE         =            %8.2f\n", energyTerms[ 8]);
      printf("reference_LYS         =            %8.2f\n", energyTerms[ 9]);
      printf("reference_LEU         =            %8.2f\n", energyTerms[10]);
      printf("reference_MET         =            %8.2f\n", energyTerms[11]);
      printf("reference_ASN         =            %8.2f\n", energyTerms[12]);
      printf("reference_PRO         =            %8.2f\n", energyTerms[13]);
      printf("reference_GLN         =            %8.2f\n", energyTerms[14]);
      printf("reference_ARG         =            %8.2f\n", energyTerms[15]);
      printf("reference_SER         =            %8.2f\n", energyTerms[16]);
      printf("reference_THR         =            %8.2f\n", energyTerms[17]);
      printf("reference_VAL         =            %8.2f\n", energyTerms[18]);
      printf("reference_TRP         =            %8.2f\n", energyTerms[19]);
      printf("reference_TYR         =            %8.2f\n", energyTerms[20]);
      printf("intraR_vdwatt         =            %8.2f\n", energyTerms[21]);
      printf("intraR_vdwrep         =            %8.2f\n", energyTerms[22]);
      printf("intraR_electr         =            %8.2f\n", energyTerms[23]);
      printf("intraR_deslvP         =            %8.2f\n", energyTerms[24]);
      printf("intraR_deslvH         =            %8.2f\n", energyTerms[25]);
      printf("intraR_hbscbb_dis     =            %8.2f\n", energyTerms[26]);
      printf("intraR_hbscbb_the     =            %8.2f\n", energyTerms[27]);
      printf("intraR_hbscbb_phi     =            %8.2f\n", energyTerms[28]);
      printf("interS_vdwatt         =            %8.2f\n", energyTerms[31]);
      printf("interS_vdwrep         =            %8.2f\n", energyTerms[32]);
      printf("interS_electr         =            %8.2f\n", energyTerms[33]);
      printf("interS_deslvP         =            %8.2f\n", energyTerms[34]);
      printf("interS_deslvH         =            %8.2f\n", energyTerms[35]);
      printf("interS_ssbond         =            %8.2f\n", energyTerms[36]);
      printf("interS_hbbbbb_dis     =            %8.2f\n", energyTerms[41]);
      printf("interS_hbbbbb_the     =            %8.2f\n", energyTerms[42]);
      printf("interS_hbbbbb_phi     =            %8.2f\n", energyTerms[43]);
      printf("interS_hbscbb_dis     =            %8.2f\n", energyTerms[44]);
      printf("interS_hbscbb_the     =            %8.2f\n", energyTerms[45]);
      printf("interS_hbscbb_phi     =            %8.2f\n", energyTerms[46]);
      printf("interS_hbscsc_dis     =            %8.2f\n", energyTerms[47]);
      printf("interS_hbscsc_the     =            %8.2f\n", energyTerms[48]);
      printf("interS_hbscsc_phi     =            %8.2f\n", energyTerms[49]);
      printf("interD_vdwatt         =            %8.2f\n", energyTerms[51]);
      printf("interD_vdwrep         =            %8.2f\n", energyTerms[52]);
      printf("interD_electr         =            %8.2f\n", energyTerms[53]);
      printf("interD_deslvP         =            %8.2f\n", energyTerms[54]);
      printf("interD_deslvH         =            %8.2f\n", energyTerms[55]);
      printf("interD_ssbond         =            %8.2f\n", energyTerms[56]);
      printf("interD_hbbbbb_dis     =            %8.2f\n", energyTerms[61]);
      printf("interD_hbbbbb_the     =            %8.2f\n", energyTerms[62]);
      printf("interD_hbbbbb_phi     =            %8.2f\n", energyTerms[63]);
      printf("interD_hbscbb_dis     =            %8.2f\n", energyTerms[64]);
      printf("interD_hbscbb_the     =            %8.2f\n", energyTerms[65]);
      printf("interD_hbscbb_phi     =            %8.2f\n", energyTerms[66]);
      printf("interD_hbscsc_dis     =            %8.2f\n", energyTerms[67]);
      printf("interD_hbscsc_the     =            %8.2f\n", energyTerms[68]);
      printf("interD_hbscsc_phi     =            %8.2f\n", energyTerms[69]);
    }

    //total energy, weighted
    EnergyTermWeighting(energyTerms);
    energyArrayOfRotamers[ir] = energyTerms[0];
    if(energyArrayOfRotamers[ir] < minEnergy){
      minEnergy = energyArrayOfRotamers[ir];
      minEnergyRotIndex = ir;
      ResidueCopy(pDesign, &tempResidue);
    }

    if(ir == 0){
      StringArrayAppend(&rotTypes, tempResidue.name);
      minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
      minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
    }
    else{
      int pos = -1;
      StringArrayFind(&rotTypes, tempResidue.name, &pos);
      if(pos == -1){
        StringArrayAppend(&rotTypes, tempResidue.name);
        minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
        minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
      }
      else{
        if(energyArrayOfRotamers[ir] < minEnergyRotTypeEnergy[pos]){
          minEnergyRotTypeIndex[pos] = ir;
          minEnergyRotTypeEnergy[pos] = energyArrayOfRotamers[ir];
        }
      }
    }

    RotamerExtract(pRotIR);
    ResidueDestroy(&tempResidue);
  }

  if(energyArrayOfRotamers != NULL){
    free(energyArrayOfRotamers);
    energyArrayOfRotamers = NULL;
  }
  if(ppSurroundingResidues != NULL){
    free(ppSurroundingResidues);
    ppSurroundingResidues = NULL;
  }
  StringArrayDestroy(&rotTypes);
  ResidueDestroy(&original);

  return Success;
}



int ProteinSiteOptimizeRotamerHBondEnergy(Structure *pStructure, int chainIndex, int resiIndex){
  DesignSite *pDesignSite = StructureFindDesignSite(pStructure,chainIndex,resiIndex);
  if(pDesignSite==NULL) return Success;
  RotamerSet *pRotSet = DesignSiteGetRotamers(pDesignSite);
  Residue *pDesign = pDesignSite->pResidue;
  double *energyArrayOfRotamers = (double *)malloc(sizeof(double)*RotamerSetGetCount(pRotSet));
  for(int i=0; i<RotamerSetGetCount(pRotSet); ++i) energyArrayOfRotamers[i]=0.0;

  //step 1: find out residues within 5 angstroms to the design site of interest;
  int surroundingResiNum = 0;
  Residue **ppSurroundingResidues = NULL;
  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    for(int j = 0; j < pChainI->residueNum; j++){
      Residue *pResi2 = pChainI->residues + j;
      if(strcmp(pDesign->name, pResi2->name) == 0 && pDesign->posInChain == pResi2->posInChain){
        continue;
      }
      else if(AtomArrayCalcMinDistance(&pDesign->atoms, &pResi2->atoms)<ENERGY_DISTANCE_CUTOFF){
        surroundingResiNum++;
        ppSurroundingResidues = (Residue **)realloc(ppSurroundingResidues, sizeof(Residue*)*surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum-1] = pResi2;
      }
    }
  }

  // step 2: calculate the energy between the rotamers of the design site
  double minEnergy = 1000.0;
  int minEnergyRotIndex = -1;
  int minEnergyRotTypeIndex[22];
  double minEnergyRotTypeEnergy[22];
  StringArray rotTypes;
  StringArrayCreate(&rotTypes);
  for(int ir = 0; ir < RotamerSetGetCount(pRotSet); ir++){
    double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
    for(int i=0; i<MAX_EVOEF_ENERGY_TERM_NUM; i++){
      energyTerms[i]=0.0;
    }
    Rotamer *pRotIR = RotamerSetGet(pRotSet, ir);
    RotamerRestore(pRotIR, pRotSet);

    Residue tempResidue;
    ResidueCreate(&tempResidue);
    ResidueCopy(&tempResidue, pDesign);
    ResidueSetName(&tempResidue, pRotIR->type);
    AtomArrayCopy(&tempResidue.atoms, &pRotIR->atoms);
    BondSetCopy(&tempResidue.bonds, &pRotIR->bonds);
    EVOEF_AminoAcidReferenceEnergy(tempResidue.name,energyTerms);
    EVOEF_EnergyResidueIntraEnergy(&tempResidue,energyTerms);
    for(int is = 0; is < surroundingResiNum; is++){
      Residue *pResIS = ppSurroundingResidues[is];
      if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain-1){
        EVOEF_EnergyResidueAndNextResidue(&tempResidue,pResIS,energyTerms);
      }
      else if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain+1){
        EVOEF_EnergyResidueAndNextResidue(pResIS,&tempResidue,energyTerms);
      }
      else{
        if(strcmp(tempResidue.chainName, pResIS->chainName)==0){
          EVOEF_EnergyResidueAndOtherResidueSameChain(&tempResidue,pResIS,energyTerms);
        }
        else{
          EVOEF_EnergyResidueAndOtherResidueDifferentChain(&tempResidue,pResIS,energyTerms);
        }
      }
    }

    if(FALSE){
      printf("ROTAMER %4d: %1s %4d %3s, ENERGY: %6.2f\n", ir, RotamerGetChainName(pRotIR), RotamerGetPosInChain(pRotIR), RotamerGetType(pRotIR), energyArrayOfRotamers[ir]);
      printf("reference_ALA         =            %8.2f\n", energyTerms[ 1]);
      printf("reference_CYS         =            %8.2f\n", energyTerms[ 2]);
      printf("reference_ASP         =            %8.2f\n", energyTerms[ 3]);
      printf("reference_GLU         =            %8.2f\n", energyTerms[ 4]);
      printf("reference_PHE         =            %8.2f\n", energyTerms[ 5]);
      printf("reference_GLY         =            %8.2f\n", energyTerms[ 6]);
      printf("reference_HIS         =            %8.2f\n", energyTerms[ 7]);
      printf("reference_ILE         =            %8.2f\n", energyTerms[ 8]);
      printf("reference_LYS         =            %8.2f\n", energyTerms[ 9]);
      printf("reference_LEU         =            %8.2f\n", energyTerms[10]);
      printf("reference_MET         =            %8.2f\n", energyTerms[11]);
      printf("reference_ASN         =            %8.2f\n", energyTerms[12]);
      printf("reference_PRO         =            %8.2f\n", energyTerms[13]);
      printf("reference_GLN         =            %8.2f\n", energyTerms[14]);
      printf("reference_ARG         =            %8.2f\n", energyTerms[15]);
      printf("reference_SER         =            %8.2f\n", energyTerms[16]);
      printf("reference_THR         =            %8.2f\n", energyTerms[17]);
      printf("reference_VAL         =            %8.2f\n", energyTerms[18]);
      printf("reference_TRP         =            %8.2f\n", energyTerms[19]);
      printf("reference_TYR         =            %8.2f\n", energyTerms[20]);
      printf("intraR_vdwatt         =            %8.2f\n", energyTerms[21]);
      printf("intraR_vdwrep         =            %8.2f\n", energyTerms[22]);
      printf("intraR_electr         =            %8.2f\n", energyTerms[23]);
      printf("intraR_deslvP         =            %8.2f\n", energyTerms[24]);
      printf("intraR_deslvH         =            %8.2f\n", energyTerms[25]);
      printf("intraR_hbscbb_dis     =            %8.2f\n", energyTerms[26]);
      printf("intraR_hbscbb_the     =            %8.2f\n", energyTerms[27]);
      printf("intraR_hbscbb_phi     =            %8.2f\n", energyTerms[28]);
      printf("interS_vdwatt         =            %8.2f\n", energyTerms[31]);
      printf("interS_vdwrep         =            %8.2f\n", energyTerms[32]);
      printf("interS_electr         =            %8.2f\n", energyTerms[33]);
      printf("interS_deslvP         =            %8.2f\n", energyTerms[34]);
      printf("interS_deslvH         =            %8.2f\n", energyTerms[35]);
      printf("interS_ssbond         =            %8.2f\n", energyTerms[36]);
      printf("interS_hbbbbb_dis     =            %8.2f\n", energyTerms[41]);
      printf("interS_hbbbbb_the     =            %8.2f\n", energyTerms[42]);
      printf("interS_hbbbbb_phi     =            %8.2f\n", energyTerms[43]);
      printf("interS_hbscbb_dis     =            %8.2f\n", energyTerms[44]);
      printf("interS_hbscbb_the     =            %8.2f\n", energyTerms[45]);
      printf("interS_hbscbb_phi     =            %8.2f\n", energyTerms[46]);
      printf("interS_hbscsc_dis     =            %8.2f\n", energyTerms[47]);
      printf("interS_hbscsc_the     =            %8.2f\n", energyTerms[48]);
      printf("interS_hbscsc_phi     =            %8.2f\n", energyTerms[49]);
      printf("interD_vdwatt         =            %8.2f\n", energyTerms[51]);
      printf("interD_vdwrep         =            %8.2f\n", energyTerms[52]);
      printf("interD_electr         =            %8.2f\n", energyTerms[53]);
      printf("interD_deslvP         =            %8.2f\n", energyTerms[54]);
      printf("interD_deslvH         =            %8.2f\n", energyTerms[55]);
      printf("interD_ssbond         =            %8.2f\n", energyTerms[56]);
      printf("interD_hbbbbb_dis     =            %8.2f\n", energyTerms[61]);
      printf("interD_hbbbbb_the     =            %8.2f\n", energyTerms[62]);
      printf("interD_hbbbbb_phi     =            %8.2f\n", energyTerms[63]);
      printf("interD_hbscbb_dis     =            %8.2f\n", energyTerms[64]);
      printf("interD_hbscbb_the     =            %8.2f\n", energyTerms[65]);
      printf("interD_hbscbb_phi     =            %8.2f\n", energyTerms[66]);
      printf("interD_hbscsc_dis     =            %8.2f\n", energyTerms[67]);
      printf("interD_hbscsc_the     =            %8.2f\n", energyTerms[68]);
      printf("interD_hbscsc_phi     =            %8.2f\n", energyTerms[69]);
    }


    EnergyTermWeighting(energyTerms);
    //only consider the non-intra residue hbond energy
    for(int i = 41; i <= 49; i++){
      energyArrayOfRotamers[ir] += energyTerms[i];
    }
    for(int i = 61; i <= 69; i++){
      energyArrayOfRotamers[ir] += energyTerms[i];
    }

    if(energyArrayOfRotamers[ir] < minEnergy){
      minEnergy = energyArrayOfRotamers[ir];
      minEnergyRotIndex = ir;
      ResidueCopy(pDesign, &tempResidue);
    }

    if(ir == 0){
      StringArrayAppend(&rotTypes, tempResidue.name);
      minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
      minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
    }
    else{
      int pos = -1;
      StringArrayFind(&rotTypes, tempResidue.name, &pos);
      if(pos == -1){
        StringArrayAppend(&rotTypes, tempResidue.name);
        minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
        minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
      }
      else{
        if(energyArrayOfRotamers[ir] < minEnergyRotTypeEnergy[pos]){
          minEnergyRotTypeIndex[pos] = ir;
          minEnergyRotTypeEnergy[pos] = energyArrayOfRotamers[ir];
        }
      }
    }

    RotamerExtract(pRotIR);
    ResidueDestroy(&tempResidue);
  }

  if(energyArrayOfRotamers != NULL){
    free(energyArrayOfRotamers);
    energyArrayOfRotamers = NULL;
  }
  if(ppSurroundingResidues != NULL){
    free(ppSurroundingResidues);
    ppSurroundingResidues = NULL;
  }
  StringArrayDestroy(&rotTypes);

  return Success;
}

int ProteinSiteCalcRotamersEnergy(Structure *pStructure,AAppTable* pAAppTable,RamaTable* pRama,int chainIndex,int resiIndex,FILE* fp){
  DesignSite *pDesignSite = StructureFindDesignSite(pStructure,chainIndex,resiIndex);
  if(pDesignSite==NULL) return Success;
  RotamerSet *pRotSet = DesignSiteGetRotamers(pDesignSite);
  Residue *pResidue = pDesignSite->pResidue;

  //step 1: find out residues within 5 angstroms to the design site of interest;
  int surroundingResiNum = 0;
  Residue **ppSurroundingResidues = NULL;
  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    for(int j = 0; j < pChainI->residueNum; j++){
      Residue *pResi2 = pChainI->residues + j;
      if(strcmp(pResidue->name, pResi2->name) == 0 && pResidue->posInChain == pResi2->posInChain){
        continue;
      }
      else if(AtomArrayCalcMinDistance(&pResidue->atoms, &pResi2->atoms)<ENERGY_DISTANCE_CUTOFF){
        surroundingResiNum++;
        ppSurroundingResidues = (Residue **)realloc(ppSurroundingResidues, sizeof(Residue*)*surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum-1] = pResi2;
      }
    }
  }

  // step 2: calculate the energy between the rotamers of the design site
  for(int ir = 0; ir < RotamerSetGetCount(pRotSet); ir++){
    double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
    EnergyTermInitialize(energyTerms);
    Rotamer *pRotIR = RotamerSetGet(pRotSet, ir);
    RotamerRestore(pRotIR, pRotSet);

    Residue tempResidue;
    ResidueCreate(&tempResidue);
    ResidueCopy(&tempResidue, pResidue);
    ResidueSetName(&tempResidue, pRotIR->type);
    AtomArrayCopy(&tempResidue.atoms, &pRotIR->atoms);
    BondSetCopy(&tempResidue.bonds, &pRotIR->bonds);
    
    //intra energy, add the statistical energy
    EVOEF_AminoAcidReferenceEnergy(tempResidue.name,energyTerms);
    EVOEF_EnergyResidueIntraEnergy(&tempResidue,energyTerms);
    EVOEF_RotamerPropensityAndRamachandranEnergy(pRotIR,pResidue,pAAppTable,pRama,energyTerms);
    EVOEF_RotamerDunbrackEnergy(pRotIR,energyTerms);

    for(int is = 0; is < surroundingResiNum; is++){
      Residue *pResIS = ppSurroundingResidues[is];
      if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain-1){
        EVOEF_EnergyResidueAndNextResidue(&tempResidue,pResIS,energyTerms);
      }
      else if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain+1){
        EVOEF_EnergyResidueAndNextResidue(pResIS,&tempResidue,energyTerms);
      }
      else{
        if(strcmp(ResidueGetChainName(&tempResidue),ResidueGetChainName(pResIS))==0){
          EVOEF_EnergyResidueAndOtherResidueSameChain(&tempResidue,pResIS,energyTerms);
        }
        else{
          if(ChainGetType(StructureFindChainByName(pStructure,pResIS->chainName))==Type_Chain_SmallMol){
            EVOEF_EnergyResidueAndLigandResidue(&tempResidue,pResIS,energyTerms);
          }
          else{
            EVOEF_EnergyResidueAndOtherResidueDifferentChain(&tempResidue,pResIS,energyTerms);
          }
        }
      }
    }

    //disable the weighting
    EnergyTermWeighting(energyTerms);
    if(RotamerAndResidueInSameType(pRotIR,pResidue)==TRUE){
      double rmsd=RotamerAndResidueSidechainRMSD(pRotIR,pResidue);
      if(rmsd<0.005){//this should be a crystal rotamer
        fprintf(fp,"POS %3d XAL",RotamerGetPosInChain(pRotIR));
        //if the energy for monomer & PPI is too large, skip the crystal rotamer
        BOOL energyTooBig=FALSE;
        for(int index=1;index<=70;index++){
          if(energyTerms[index]>30.0 || energyTerms[index]<-30.0){
            energyTooBig=TRUE;
            break;
          }
        }
        if(energyTooBig==TRUE){
          char errMsg[MAX_LENGTH_ERR_MSG+1];
          sprintf(errMsg,"in file %s function %s line %d, the absolute energy value of crystal rotamer on residue %s%d%s is too large",__FILE__,__FUNCTION__,__LINE__,
            ResidueGetChainName(&tempResidue),ResidueGetPosInChain(&tempResidue),ResidueGetName(&tempResidue));
          TraceError(errMsg,Warning);
          //RotamerExtract(pRotIR);
          //ResidueDestroy(&tempResidue);
          //continue;
        }

      }
      else{
        fprintf(fp,"POS %3d NAT",RotamerGetPosInChain(pRotIR));
      }
    }
    else{
      fprintf(fp,"POS %3d %s",RotamerGetPosInChain(pRotIR),RotamerGetType(pRotIR));
    }
    //20 reference energy
    for(int index=1; index<=20; index++){
      fprintf(fp," %8.2f",energyTerms[index]);
    }
    //intraR vdw->desol
    for(int index=21; index<=25; index++){
      fprintf(fp," %8.2f",energyTerms[index]);
    }
    //intraR hbonds
    for(int index=26; index<=28; index++){
      fprintf(fp," %8.2f",energyTerms[index]);
    }
    //intraR aa propensity, rama and dunbrack
    for(int index=91; index<=93; index++){
      fprintf(fp," %8.2f",energyTerms[index]);
    }
    //interS vdw->desol
    for(int index=31; index<=36; index++){
      fprintf(fp," %8.2f",energyTerms[index]);
    }
    //interS hbonds
    for(int index=41; index<=49; index++){
      fprintf(fp," %8.2f",energyTerms[index]);
    }
    //interD vdw->desol
    for(int index=51; index<=56; index++){
      fprintf(fp," %8.2f",energyTerms[index]);
    }
    //interD hbonds
    for(int index=61; index<=69; index++){
      fprintf(fp," %8.2f",energyTerms[index]);
    }
    //prot-lig vdw->desol
    for(int index=71; index<=75; index++){
      fprintf(fp," %8.2f",energyTerms[index]);
    }
    //prot-lig hbonds
    for(int index=81; index<=86; index++){
      fprintf(fp," %8.2f",energyTerms[index]);
    }
    fprintf(fp,"\n");

    RotamerExtract(pRotIR);
    ResidueDestroy(&tempResidue);
  }

  if(ppSurroundingResidues != NULL){
    free(ppSurroundingResidues);
    ppSurroundingResidues = NULL;
  }

  return Success;
}


int ProteinSiteOptimizeRotamerWithBBdepRotLib(Structure *pStructure, int chainIndex, int resiIndex,BBdepRotamerLib *pBBdepRotLib){
  DesignSite *pDesignSite = StructureFindDesignSite(pStructure,chainIndex,resiIndex);
  if(pDesignSite==NULL) return Success;
  RotamerSet *pRotSet = DesignSiteGetRotamers(pDesignSite);
  Residue *pDesign = pDesignSite->pResidue;
  double *energyArrayOfRotamers = (double *)malloc(sizeof(double)*RotamerSetGetCount(pRotSet));
  for(int i=0; i<RotamerSetGetCount(pRotSet); ++i) energyArrayOfRotamers[i]=0.0;

  //step 1: find out residues within 5 angstroms to the design site of interest;
  int surroundingResiNum = 0;
  Residue **ppSurroundingResidues = NULL;
  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    for(int j = 0; j < pChainI->residueNum; j++){
      Residue *pResi2 = pChainI->residues + j;
      if(strcmp(pDesign->name, pResi2->name) == 0 && pDesign->posInChain == pResi2->posInChain){
        continue;
      }
      else if(AtomArrayCalcMinDistance(&pDesign->atoms, &pResi2->atoms)<ENERGY_DISTANCE_CUTOFF){
        surroundingResiNum++;
        ppSurroundingResidues = (Residue **)realloc(ppSurroundingResidues, sizeof(Residue*)*surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum-1] = pResi2;
      }
    }
  }

  // step 2: calculate the energy between the rotamers of the design site
  double minEnergy = 1000.0;
  int minEnergyRotIndex = -1;
  int minEnergyRotTypeIndex[22];
  double minEnergyRotTypeEnergy[22];
  StringArray rotTypes;
  StringArrayCreate(&rotTypes);
  for(int ir = 0; ir < RotamerSetGetCount(pRotSet); ir++){
    double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
    EnergyTermInitialize(energyTerms);
    Rotamer *pRotIR = RotamerSetGet(pRotSet, ir);
    RotamerRestore(pRotIR, pRotSet);

    Residue tempResidue;
    ResidueCreate(&tempResidue);
    ResidueCopy(&tempResidue, pDesign);
    ResidueSetName(&tempResidue, pRotIR->type);
    AtomArrayCopy(&tempResidue.atoms, &pRotIR->atoms);
    BondSetCopy(&tempResidue.bonds, &pRotIR->bonds);
    EVOEF_AminoAcidReferenceEnergy(tempResidue.name,energyTerms);
    EVOEF_EnergyResidueIntraEnergy(&tempResidue,energyTerms);
    //AminoAcidDunbrackEnergy(&tempResidue,pBBdepRotLib,energyTerms);
    EVOEF_RotamerDunbrackEnergy(pRotIR,energyTerms);
    tempResidue.dunbrack = pRotIR->dunbrack;
    DoubleArrayCopy(&tempResidue.xtorsions,&pRotIR->xtorsions);
    for(int is = 0; is < surroundingResiNum; is++){
      Residue *pResIS = ppSurroundingResidues[is];
      if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain-1){
        EVOEF_EnergyResidueAndNextResidue(&tempResidue,pResIS,energyTerms);
      }
      else if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain+1){
        EVOEF_EnergyResidueAndNextResidue(pResIS,&tempResidue,energyTerms);
      }
      else{
        if(strcmp(ResidueGetChainName(&tempResidue),ResidueGetChainName(pResIS))==0){
          EVOEF_EnergyResidueAndOtherResidueSameChain(&tempResidue,pResIS,energyTerms);
        }
        else{
          EVOEF_EnergyResidueAndOtherResidueDifferentChain(&tempResidue,pResIS,energyTerms);
        }
      }
    }

    if(FALSE){
      //for debug: energy terms, not weighted
      printf("ROTAMER %4d: %1s %4d %3s, ENERGY: %6.2f\n", ir, RotamerGetChainName(pRotIR), RotamerGetPosInChain(pRotIR), RotamerGetType(pRotIR), energyArrayOfRotamers[ir]);
      printf("reference_ALA         =            %8.2f\n", energyTerms[ 1]);  
      printf("reference_CYS         =            %8.2f\n", energyTerms[ 2]);
      printf("reference_ASP         =            %8.2f\n", energyTerms[ 3]);
      printf("reference_GLU         =            %8.2f\n", energyTerms[ 4]);
      printf("reference_PHE         =            %8.2f\n", energyTerms[ 5]);
      printf("reference_GLY         =            %8.2f\n", energyTerms[ 6]);
      printf("reference_HIS         =            %8.2f\n", energyTerms[ 7]);
      printf("reference_ILE         =            %8.2f\n", energyTerms[ 8]);
      printf("reference_LYS         =            %8.2f\n", energyTerms[ 9]);
      printf("reference_LEU         =            %8.2f\n", energyTerms[10]);
      printf("reference_MET         =            %8.2f\n", energyTerms[11]);
      printf("reference_ASN         =            %8.2f\n", energyTerms[12]);
      printf("reference_PRO         =            %8.2f\n", energyTerms[13]);
      printf("reference_GLN         =            %8.2f\n", energyTerms[14]);
      printf("reference_ARG         =            %8.2f\n", energyTerms[15]);
      printf("reference_SER         =            %8.2f\n", energyTerms[16]);
      printf("reference_THR         =            %8.2f\n", energyTerms[17]);
      printf("reference_VAL         =            %8.2f\n", energyTerms[18]);
      printf("reference_TRP         =            %8.2f\n", energyTerms[19]);
      printf("reference_TYR         =            %8.2f\n", energyTerms[20]);

      printf("intraR_vdwatt         =            %8.2f\n", energyTerms[21]);
      printf("intraR_vdwrep         =            %8.2f\n", energyTerms[22]);
      printf("intraR_electr         =            %8.2f\n", energyTerms[23]);
      printf("intraR_deslvP         =            %8.2f\n", energyTerms[24]);
      printf("intraR_deslvH         =            %8.2f\n", energyTerms[25]);
      printf("intraR_hbscbb_dis     =            %8.2f\n", energyTerms[26]);
      printf("intraR_hbscbb_the     =            %8.2f\n", energyTerms[27]);
      printf("intraR_hbscbb_phi     =            %8.2f\n", energyTerms[28]);
      printf("aapropensity          =            %8.2f\n", energyTerms[91]);
      printf("ramachandran          =            %8.2f\n", energyTerms[92]);
      printf("dunbrack              =            %8.2f\n", energyTerms[93]);

      printf("interS_vdwatt         =            %8.2f\n", energyTerms[31]);
      printf("interS_vdwrep         =            %8.2f\n", energyTerms[32]);
      printf("interS_electr         =            %8.2f\n", energyTerms[33]);
      printf("interS_deslvP         =            %8.2f\n", energyTerms[34]);
      printf("interS_deslvH         =            %8.2f\n", energyTerms[35]);
      printf("interS_ssbond         =            %8.2f\n", energyTerms[36]);
      printf("interS_hbbbbb_dis     =            %8.2f\n", energyTerms[41]);
      printf("interS_hbbbbb_the     =            %8.2f\n", energyTerms[42]);
      printf("interS_hbbbbb_phi     =            %8.2f\n", energyTerms[43]);
      printf("interS_hbscbb_dis     =            %8.2f\n", energyTerms[44]);
      printf("interS_hbscbb_the     =            %8.2f\n", energyTerms[45]);
      printf("interS_hbscbb_phi     =            %8.2f\n", energyTerms[46]);
      printf("interS_hbscsc_dis     =            %8.2f\n", energyTerms[47]);
      printf("interS_hbscsc_the     =            %8.2f\n", energyTerms[48]);
      printf("interS_hbscsc_phi     =            %8.2f\n", energyTerms[49]);

      printf("interD_vdwatt         =            %8.2f\n", energyTerms[51]);
      printf("interD_vdwrep         =            %8.2f\n", energyTerms[52]);
      printf("interD_electr         =            %8.2f\n", energyTerms[53]);
      printf("interD_deslvP         =            %8.2f\n", energyTerms[54]);
      printf("interD_deslvH         =            %8.2f\n", energyTerms[55]);
      printf("interD_ssbond         =            %8.2f\n", energyTerms[56]);
      printf("interD_hbbbbb_dis     =            %8.2f\n", energyTerms[61]);
      printf("interD_hbbbbb_the     =            %8.2f\n", energyTerms[62]);
      printf("interD_hbbbbb_phi     =            %8.2f\n", energyTerms[63]);
      printf("interD_hbscbb_dis     =            %8.2f\n", energyTerms[64]);
      printf("interD_hbscbb_the     =            %8.2f\n", energyTerms[65]);
      printf("interD_hbscbb_phi     =            %8.2f\n", energyTerms[66]);
      printf("interD_hbscsc_dis     =            %8.2f\n", energyTerms[67]);
      printf("interD_hbscsc_the     =            %8.2f\n", energyTerms[68]);
      printf("interD_hbscsc_phi     =            %8.2f\n", energyTerms[69]);      
    }

    //total energy: weighted
    EnergyTermWeighting(energyTerms);
    energyArrayOfRotamers[ir] = energyTerms[0];

    if(energyArrayOfRotamers[ir] < minEnergy){
      minEnergy = energyArrayOfRotamers[ir];
      minEnergyRotIndex = ir;
      ResidueCopy(pDesign, &tempResidue);
    }

    if(ir == 0){
      StringArrayAppend(&rotTypes, tempResidue.name);
      minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
      minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
    }
    else{
      int pos = -1;
      StringArrayFind(&rotTypes, tempResidue.name, &pos);
      if(pos == -1){
        StringArrayAppend(&rotTypes, tempResidue.name);
        minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
        minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
      }
      else{
        if(energyArrayOfRotamers[ir] < minEnergyRotTypeEnergy[pos]){
          minEnergyRotTypeIndex[pos] = ir;
          minEnergyRotTypeEnergy[pos] = energyArrayOfRotamers[ir];
        }
      }
    }

    RotamerExtract(pRotIR);
    ResidueDestroy(&tempResidue);
  }

  if(energyArrayOfRotamers != NULL){
    free(energyArrayOfRotamers);
    energyArrayOfRotamers = NULL;
  }
  if(ppSurroundingResidues != NULL){
    free(ppSurroundingResidues);
    ppSurroundingResidues = NULL;
  }
  StringArrayDestroy(&rotTypes);

  return Success;
}

int ProteinSiteOptimizeRotamerLocallyWithBBdepRotLib(Structure *pStructure, int chainIndex, int resiIndex,double rmsdcutoff,BBdepRotamerLib *pBBdepRotLib){
  DesignSite *pDesignSite = StructureFindDesignSite(pStructure,chainIndex,resiIndex);
  if(pDesignSite==NULL) return Success;
  RotamerSet *pRotSet = DesignSiteGetRotamers(pDesignSite);
  Residue *pDesign = pDesignSite->pResidue;
  double *energyArrayOfRotamers = (double *)malloc(sizeof(double)*RotamerSetGetCount(pRotSet));
  for(int i=0; i<RotamerSetGetCount(pRotSet); ++i) energyArrayOfRotamers[i]=0.0;

  //step 1: find out residues within 6 angstroms to the design site of interest;
  int surroundingResiNum = 0;
  Residue **ppSurroundingResidues = NULL;
  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    for(int j = 0; j < pChainI->residueNum; j++){
      Residue *pResi2 = pChainI->residues + j;
      if(strcmp(pDesign->name, pResi2->name) == 0 && pDesign->posInChain == pResi2->posInChain){
        continue;
      }
      else if(AtomArrayCalcMinDistance(&pDesign->atoms, &pResi2->atoms)<ENERGY_DISTANCE_CUTOFF){
        surroundingResiNum++;
        ppSurroundingResidues = (Residue **)realloc(ppSurroundingResidues, sizeof(Residue*)*surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum-1] = pResi2;
      }
    }
  }

  Residue original;
  ResidueCreate(&original);
  ResidueCopy(&original,pDesign);

  // step 2: calculate the energy between the rotamers of the design site
  double minEnergy = 1000.0;
  int minEnergyRotIndex = -1;
  int minEnergyRotTypeIndex[22];
  double minEnergyRotTypeEnergy[22];
  StringArray rotTypes;
  StringArrayCreate(&rotTypes);
  for(int ir = 0; ir < RotamerSetGetCount(pRotSet); ir++){
    double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
    EnergyTermInitialize(energyTerms);
    Rotamer *pRotIR = RotamerSetGet(pRotSet, ir);
    RotamerRestore(pRotIR, pRotSet);

    Residue tempResidue;
    ResidueCreate(&tempResidue);
    ResidueCopy(&tempResidue, pDesign);
    ResidueSetName(&tempResidue, pRotIR->type);
    AtomArrayCopy(&tempResidue.atoms, &pRotIR->atoms);
    BondSetCopy(&tempResidue.bonds, &pRotIR->bonds);
    //skip the rotamers that has large side-chain shift
    double realrmsd=ResidueAndResidueSidechainRMSD(&tempResidue,&original);
    if(realrmsd>rmsdcutoff && realrmsd<1e4){
      ResidueDestroy(&tempResidue);
      continue;
    }
    EVOEF_AminoAcidReferenceEnergy(tempResidue.name,energyTerms);
    EVOEF_EnergyResidueIntraEnergy(&tempResidue,energyTerms);
    AminoAcidDunbrackEnergy(&tempResidue,pBBdepRotLib,energyTerms);
    for(int is = 0; is < surroundingResiNum; is++){
      Residue *pResIS = ppSurroundingResidues[is];
      if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain-1){
        EVOEF_EnergyResidueAndNextResidue(&tempResidue,pResIS,energyTerms);
      }
      else if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain+1){
        EVOEF_EnergyResidueAndNextResidue(pResIS,&tempResidue,energyTerms);
      }
      else{
        if(strcmp(tempResidue.chainName, pResIS->chainName)==0){
          EVOEF_EnergyResidueAndOtherResidueSameChain(&tempResidue,pResIS,energyTerms);
        }
        else{
          EVOEF_EnergyResidueAndOtherResidueDifferentChain(&tempResidue,pResIS,energyTerms);
        }
      }
    }

    if(FALSE){
      printf("ROTAMER %4d: %1s %4d %3s, ENERGY: %6.2f\n", ir, RotamerGetChainName(pRotIR), RotamerGetPosInChain(pRotIR), RotamerGetType(pRotIR), energyArrayOfRotamers[ir]);
      printf("reference_ALA         =            %8.2f\n", energyTerms[ 1]);
      printf("reference_CYS         =            %8.2f\n", energyTerms[ 2]);
      printf("reference_ASP         =            %8.2f\n", energyTerms[ 3]);
      printf("reference_GLU         =            %8.2f\n", energyTerms[ 4]);
      printf("reference_PHE         =            %8.2f\n", energyTerms[ 5]);
      printf("reference_GLY         =            %8.2f\n", energyTerms[ 6]);
      printf("reference_HIS         =            %8.2f\n", energyTerms[ 7]);
      printf("reference_ILE         =            %8.2f\n", energyTerms[ 8]);
      printf("reference_LYS         =            %8.2f\n", energyTerms[ 9]);
      printf("reference_LEU         =            %8.2f\n", energyTerms[10]);
      printf("reference_MET         =            %8.2f\n", energyTerms[11]);
      printf("reference_ASN         =            %8.2f\n", energyTerms[12]);
      printf("reference_PRO         =            %8.2f\n", energyTerms[13]);
      printf("reference_GLN         =            %8.2f\n", energyTerms[14]);
      printf("reference_ARG         =            %8.2f\n", energyTerms[15]);
      printf("reference_SER         =            %8.2f\n", energyTerms[16]);
      printf("reference_THR         =            %8.2f\n", energyTerms[17]);
      printf("reference_VAL         =            %8.2f\n", energyTerms[18]);
      printf("reference_TRP         =            %8.2f\n", energyTerms[19]);
      printf("reference_TYR         =            %8.2f\n", energyTerms[20]);
      printf("intraR_vdwatt         =            %8.2f\n", energyTerms[21]);
      printf("intraR_vdwrep         =            %8.2f\n", energyTerms[22]);
      printf("intraR_electr         =            %8.2f\n", energyTerms[23]);
      printf("intraR_deslvP         =            %8.2f\n", energyTerms[24]);
      printf("intraR_deslvH         =            %8.2f\n", energyTerms[25]);
      printf("intraR_hbscbb_dis     =            %8.2f\n", energyTerms[26]);
      printf("intraR_hbscbb_the     =            %8.2f\n", energyTerms[27]);
      printf("intraR_hbscbb_phi     =            %8.2f\n", energyTerms[28]);
      printf("aapropensity          =            %8.2f\n", energyTerms[91]);
      printf("ramachandran          =            %8.2f\n", energyTerms[92]);
      printf("dunbrack              =            %8.2f\n", energyTerms[93]);
      printf("interS_vdwatt         =            %8.2f\n", energyTerms[31]);
      printf("interS_vdwrep         =            %8.2f\n", energyTerms[32]);
      printf("interS_electr         =            %8.2f\n", energyTerms[33]);
      printf("interS_deslvP         =            %8.2f\n", energyTerms[34]);
      printf("interS_deslvH         =            %8.2f\n", energyTerms[35]);
      printf("interS_ssbond         =            %8.2f\n", energyTerms[36]);
      printf("interS_hbbbbb_dis     =            %8.2f\n", energyTerms[41]);
      printf("interS_hbbbbb_the     =            %8.2f\n", energyTerms[42]);
      printf("interS_hbbbbb_phi     =            %8.2f\n", energyTerms[43]);
      printf("interS_hbscbb_dis     =            %8.2f\n", energyTerms[44]);
      printf("interS_hbscbb_the     =            %8.2f\n", energyTerms[45]);
      printf("interS_hbscbb_phi     =            %8.2f\n", energyTerms[46]);
      printf("interS_hbscsc_dis     =            %8.2f\n", energyTerms[47]);
      printf("interS_hbscsc_the     =            %8.2f\n", energyTerms[48]);
      printf("interS_hbscsc_phi     =            %8.2f\n", energyTerms[49]);
      printf("interD_vdwatt         =            %8.2f\n", energyTerms[51]);
      printf("interD_vdwrep         =            %8.2f\n", energyTerms[52]);
      printf("interD_electr         =            %8.2f\n", energyTerms[53]);
      printf("interD_deslvP         =            %8.2f\n", energyTerms[54]);
      printf("interD_deslvH         =            %8.2f\n", energyTerms[55]);
      printf("interD_ssbond         =            %8.2f\n", energyTerms[56]);
      printf("interD_hbbbbb_dis     =            %8.2f\n", energyTerms[61]);
      printf("interD_hbbbbb_the     =            %8.2f\n", energyTerms[62]);
      printf("interD_hbbbbb_phi     =            %8.2f\n", energyTerms[63]);
      printf("interD_hbscbb_dis     =            %8.2f\n", energyTerms[64]);
      printf("interD_hbscbb_the     =            %8.2f\n", energyTerms[65]);
      printf("interD_hbscbb_phi     =            %8.2f\n", energyTerms[66]);
      printf("interD_hbscsc_dis     =            %8.2f\n", energyTerms[67]);
      printf("interD_hbscsc_the     =            %8.2f\n", energyTerms[68]);
      printf("interD_hbscsc_phi     =            %8.2f\n", energyTerms[69]);
    }

    //total energy, weighted
    EnergyTermWeighting(energyTerms);
    energyArrayOfRotamers[ir] = energyTerms[0];
    if(energyArrayOfRotamers[ir] < minEnergy){
      minEnergy = energyArrayOfRotamers[ir];
      minEnergyRotIndex = ir;
      ResidueCopy(pDesign, &tempResidue);
    }

    if(ir == 0){
      StringArrayAppend(&rotTypes, tempResidue.name);
      minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
      minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
    }
    else{
      int pos = -1;
      StringArrayFind(&rotTypes, tempResidue.name, &pos);
      if(pos == -1){
        StringArrayAppend(&rotTypes, tempResidue.name);
        minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
        minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
      }
      else{
        if(energyArrayOfRotamers[ir] < minEnergyRotTypeEnergy[pos]){
          minEnergyRotTypeIndex[pos] = ir;
          minEnergyRotTypeEnergy[pos] = energyArrayOfRotamers[ir];
        }
      }
    }

    RotamerExtract(pRotIR);
    ResidueDestroy(&tempResidue);
  }

  if(energyArrayOfRotamers != NULL){
    free(energyArrayOfRotamers);
    energyArrayOfRotamers = NULL;
  }
  if(ppSurroundingResidues != NULL){
    free(ppSurroundingResidues);
    ppSurroundingResidues = NULL;
  }
  StringArrayDestroy(&rotTypes);
  ResidueDestroy(&original);

  return Success;
}