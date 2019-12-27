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

#include "EnergyFunction.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>


double weights[MAX_EVOEF_ENERGY_TERM_NUM];

#define DEAL_WITH_ENERGY_WEIGHTING
int EnergyTermInitialize(double *energyTerms){
  memset(energyTerms,0,sizeof(double)*MAX_EVOEF_ENERGY_TERM_NUM);
  return Success;
}

int EnergyTermWeighting(double *energyTerms){
  for(int i = 1; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    energyTerms[i] *= weights[i];
    energyTerms[0]+=energyTerms[i];
  }
  return Success;
}


int EnergyWeightRead(char* weightfile){
  for(int i=0;i<MAX_EVOEF_ENERGY_TERM_NUM;i++){
    weights[i]=1.0;
  }
  FILE* pf=fopen(weightfile,"r");
  if(pf!=NULL){
    char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pf)){
      char term[MAX_LENGTH_ONE_LINE_IN_FILE+1];
      double val=0.0;
      sscanf(line,"%s %lf",term,&val);
      if(!strcmp(term,"reference_ALA")) weights[ 1]=val;
      else if(!strcmp(term,"reference_CYS")) weights[ 2]=val;
      else if(!strcmp(term,"reference_ASP")) weights[ 3]=val;
      else if(!strcmp(term,"reference_GLU")) weights[ 4]=val;
      else if(!strcmp(term,"reference_PHE")) weights[ 5]=val;
      else if(!strcmp(term,"reference_GLY")) weights[ 6]=val;
      else if(!strcmp(term,"reference_HIS")) weights[ 7]=val;
      else if(!strcmp(term,"reference_ILE")) weights[ 8]=val;
      else if(!strcmp(term,"reference_LYS")) weights[ 9]=val;
      else if(!strcmp(term,"reference_LEU")) weights[10]=val;
      else if(!strcmp(term,"reference_MET")) weights[11]=val;
      else if(!strcmp(term,"reference_ASN")) weights[12]=val;
      else if(!strcmp(term,"reference_PRO")) weights[13]=val;
      else if(!strcmp(term,"reference_GLN")) weights[14]=val;
      else if(!strcmp(term,"reference_ARG")) weights[15]=val;
      else if(!strcmp(term,"reference_SER")) weights[16]=val;
      else if(!strcmp(term,"reference_THR")) weights[17]=val;
      else if(!strcmp(term,"reference_VAL")) weights[18]=val;
      else if(!strcmp(term,"reference_TRP")) weights[19]=val;
      else if(!strcmp(term,"reference_TYR")) weights[20]=val;

      else if(!strcmp(term,"intraR_vdwatt")) weights[21]=val;
      else if(!strcmp(term,"intraR_vdwrep")) weights[22]=val;
      else if(!strcmp(term,"intraR_electr")) weights[23]=val;
      else if(!strcmp(term,"intraR_deslvP")) weights[24]=val;
      else if(!strcmp(term,"intraR_deslvH")) weights[25]=val;
      else if(!strcmp(term,"intraR_hbscbb_dis")) weights[26]=val;
      else if(!strcmp(term,"intraR_hbscbb_the")) weights[27]=val;
      else if(!strcmp(term,"intraR_hbscbb_phi")) weights[28]=val;
      else if(!strcmp(term,"aapropensity")) weights[91]=val;
      else if(!strcmp(term,"ramachandran")) weights[92]=val;
      else if(!strcmp(term,"dunbrack"))     weights[93]=val;

      else if(!strcmp(term,"interS_vdwatt")) weights[31]=val;
      else if(!strcmp(term,"interS_vdwrep")) weights[32]=val;
      else if(!strcmp(term,"interS_electr")) weights[33]=val;
      else if(!strcmp(term,"interS_deslvP")) weights[34]=val;
      else if(!strcmp(term,"interS_deslvH")) weights[35]=val;
      else if(!strcmp(term,"interS_ssbond")) weights[36]=val;
      else if(!strcmp(term,"interS_hbbbbb_dis")) weights[41]=val;
      else if(!strcmp(term,"interS_hbbbbb_the")) weights[42]=val;
      else if(!strcmp(term,"interS_hbbbbb_phi")) weights[43]=val;
      else if(!strcmp(term,"interS_hbscbb_dis")) weights[44]=val;
      else if(!strcmp(term,"interS_hbscbb_the")) weights[45]=val;
      else if(!strcmp(term,"interS_hbscbb_phi")) weights[46]=val;
      else if(!strcmp(term,"interS_hbscsc_dis")) weights[47]=val;
      else if(!strcmp(term,"interS_hbscsc_the")) weights[48]=val;
      else if(!strcmp(term,"interS_hbscsc_phi")) weights[49]=val;

      else if(!strcmp(term,"interD_vdwatt")) weights[51]=val;
      else if(!strcmp(term,"interD_vdwrep")) weights[52]=val;
      else if(!strcmp(term,"interD_electr")) weights[53]=val;
      else if(!strcmp(term,"interD_deslvP")) weights[54]=val;
      else if(!strcmp(term,"interD_deslvH")) weights[55]=val;
      else if(!strcmp(term,"interD_ssbond")) weights[56]=val;
      else if(!strcmp(term,"interD_hbbbbb_dis")) weights[61]=val;
      else if(!strcmp(term,"interD_hbbbbb_the")) weights[62]=val;
      else if(!strcmp(term,"interD_hbbbbb_phi")) weights[63]=val;
      else if(!strcmp(term,"interD_hbscbb_dis")) weights[64]=val;
      else if(!strcmp(term,"interD_hbscbb_the")) weights[65]=val;
      else if(!strcmp(term,"interD_hbscbb_phi")) weights[66]=val;
      else if(!strcmp(term,"interD_hbscsc_dis")) weights[67]=val;
      else if(!strcmp(term,"interD_hbscsc_the")) weights[68]=val;
      else if(!strcmp(term,"interD_hbscsc_phi")) weights[69]=val;

      else if(!strcmp(term,"ligand_vdwatt")) weights[71]=val;
      else if(!strcmp(term,"ligand_vdwrep")) weights[72]=val;
      else if(!strcmp(term,"ligand_electr")) weights[73]=val;
      else if(!strcmp(term,"ligand_deslvP")) weights[74]=val;
      else if(!strcmp(term,"ligand_deslvH")) weights[75]=val;
      else if(!strcmp(term,"ligand_hbscbb_dis")) weights[84]=val;
      else if(!strcmp(term,"ligand_hbscbb_the")) weights[85]=val;
      else if(!strcmp(term,"ligand_hbscbb_phi")) weights[86]=val;
      else if(!strcmp(term,"ligand_hbscsc_dis")) weights[87]=val;
      else if(!strcmp(term,"ligand_hbscsc_the")) weights[88]=val;
      else if(!strcmp(term,"ligand_hbscsc_phi")) weights[89]=val;
    }
    fclose(pf);
  }
  else{
    char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sprintf(errMsg,"in file %s function %s line %d, cannot open weight file %s",__FILE__,__FUNCTION__,__LINE__,weightfile);
    TraceError(errMsg,IOError);
    return IOError;
  }

  return Success;
}

int EnergyWeightWrite(char* weightfile){
  FILE* fo=fopen(weightfile,"w");
  fprintf(fo,"reference_ALA      %8.3f\n",weights[ 1]);
  fprintf(fo,"reference_CYS      %8.3f\n",weights[ 2]);
  fprintf(fo,"reference_ASP      %8.3f\n",weights[ 3]);
  fprintf(fo,"reference_GLU      %8.3f\n",weights[ 4]);
  fprintf(fo,"reference_PHE      %8.3f\n",weights[ 5]);
  fprintf(fo,"reference_GLY      %8.3f\n",weights[ 6]);
  fprintf(fo,"reference_HIS      %8.3f\n",weights[ 7]);
  fprintf(fo,"reference_ILE      %8.3f\n",weights[ 8]);
  fprintf(fo,"reference_LYS      %8.3f\n",weights[ 9]);
  fprintf(fo,"reference_LEU      %8.3f\n",weights[10]);
  fprintf(fo,"reference_MET      %8.3f\n",weights[11]);
  fprintf(fo,"reference_ASN      %8.3f\n",weights[12]);
  fprintf(fo,"reference_PRO      %8.3f\n",weights[13]);
  fprintf(fo,"reference_GLN      %8.3f\n",weights[14]);
  fprintf(fo,"reference_ARG      %8.3f\n",weights[15]);
  fprintf(fo,"reference_SER      %8.3f\n",weights[16]);
  fprintf(fo,"reference_THR      %8.3f\n",weights[17]);
  fprintf(fo,"reference_VAL      %8.3f\n",weights[18]);
  fprintf(fo,"reference_TRP      %8.3f\n",weights[19]);
  fprintf(fo,"reference_TYR      %8.3f\n",weights[20]);

  fprintf(fo,"intraR_vdwatt      %8.3f\n",weights[21]);
  fprintf(fo,"intraR_vdwrep      %8.3f\n",weights[22]);
  fprintf(fo,"intraR_electr      %8.3f\n",weights[23]);
  fprintf(fo,"intraR_deslvP      %8.3f\n",weights[24]);
  fprintf(fo,"intraR_deslvH      %8.3f\n",weights[25]);
  fprintf(fo,"intraR_hbscbb_dis  %8.3f\n",weights[26]);
  fprintf(fo,"intraR_hbscbb_the  %8.3f\n",weights[27]);
  fprintf(fo,"intraR_hbscbb_phi  %8.3f\n",weights[28]);
  fprintf(fo,"aapropensity       %8.3f\n",weights[91]);
  fprintf(fo,"ramachandran       %8.3f\n",weights[92]);
  fprintf(fo,"dunbrack           %8.3f\n",weights[93]);

  fprintf(fo,"interS_vdwatt      %8.3f\n",weights[31]);
  fprintf(fo,"interS_vdwrep      %8.3f\n",weights[32]);
  fprintf(fo,"interS_electr      %8.3f\n",weights[33]);
  fprintf(fo,"interS_deslvP      %8.3f\n",weights[34]);
  fprintf(fo,"interS_deslvH      %8.3f\n",weights[35]);
  fprintf(fo,"interS_ssbond      %8.3f\n",weights[36]);
  fprintf(fo,"interS_hbbbbb_dis  %8.3f\n",weights[41]);
  fprintf(fo,"interS_hbbbbb_the  %8.3f\n",weights[42]);
  fprintf(fo,"interS_hbbbbb_phi  %8.3f\n",weights[43]);
  fprintf(fo,"interS_hbscbb_dis  %8.3f\n",weights[44]);
  fprintf(fo,"interS_hbscbb_the  %8.3f\n",weights[45]);
  fprintf(fo,"interS_hbscbb_phi  %8.3f\n",weights[46]);
  fprintf(fo,"interS_hbscsc_dis  %8.3f\n",weights[47]);
  fprintf(fo,"interS_hbscsc_the  %8.3f\n",weights[48]);
  fprintf(fo,"interS_hbscsc_phi  %8.3f\n",weights[49]);

  fprintf(fo,"interD_vdwatt      %8.3f\n",weights[51]);
  fprintf(fo,"interD_vdwrep      %8.3f\n",weights[52]);
  fprintf(fo,"interD_electr      %8.3f\n",weights[53]);
  fprintf(fo,"interD_deslvP      %8.3f\n",weights[54]);
  fprintf(fo,"interD_deslvH      %8.3f\n",weights[55]);
  fprintf(fo,"interD_ssbond      %8.3f\n",weights[56]);
  fprintf(fo,"interD_hbbbbb_dis  %8.3f\n",weights[61]);
  fprintf(fo,"interD_hbbbbb_the  %8.3f\n",weights[62]);
  fprintf(fo,"intraD_hbbbbb_phi  %8.3f\n",weights[63]);
  fprintf(fo,"interD_hbscbb_dis  %8.3f\n",weights[64]);
  fprintf(fo,"interD_hbscbb_the  %8.3f\n",weights[65]);
  fprintf(fo,"interD_hbscbb_phi  %8.3f\n",weights[66]);
  fprintf(fo,"interD_hbscsc_dis  %8.3f\n",weights[67]);
  fprintf(fo,"interD_hbscsc_the  %8.3f\n",weights[68]);
  fprintf(fo,"interD_hbscsc_phi  %8.3f\n",weights[69]);

  fprintf(fo,"ligand_vdwatt      %8.3f\n",weights[71]);
  fprintf(fo,"ligand_vdwrep      %8.3f\n",weights[72]);
  fprintf(fo,"ligand_electr      %8.3f\n",weights[73]);
  fprintf(fo,"ligand_deslvP      %8.3f\n",weights[74]);
  fprintf(fo,"ligand_deslvH      %8.3f\n",weights[75]);
  fprintf(fo,"ligand_hbscbb_dis  %8.3f\n",weights[84]);
  fprintf(fo,"ligand_hbscbb_the  %8.3f\n",weights[85]);
  fprintf(fo,"ligand_hbscbb_phi  %8.3f\n",weights[86]);
  fprintf(fo,"ligand_hbscsc_dis  %8.3f\n",weights[87]);
  fprintf(fo,"ligand_hbscsc_the  %8.3f\n",weights[88]);
  fprintf(fo,"ligand_hbscsc_phi  %8.3f\n",weights[89]);
  fclose(fo);
  return Success;
}


double CalcResidueBuriedRatio(Residue* pResi1){
  if(pResi1->nCbIn8A >= 15) return 1.0;
  else if(pResi1->nCbIn8A <= 8) return 0.0;
  else return (pResi1->nCbIn8A-8.0)/7.0;
}

double CalcAverageBuriedRatio(double ratio1, double ratio2){
  return (ratio1+ratio2)/2.0;
}

#define DEAL_WITH_BOND_CONNECTION
//these functions are used to check the 12, 13, 14 and 15 bond connectivity
BOOL ResidueIntraBond12Check(char *atom1, char *atom2, BondSet* pBondSet){
  for(int i = 0; i < pBondSet->count; i++){
    Bond *pBond = pBondSet->bonds+i;
    if( (strcmp(atom1, pBond->atomFromName) == 0 && strcmp(atom2, pBond->atomToName) == 0) ||
      (strcmp(atom2, pBond->atomFromName) == 0 && strcmp(atom1, pBond->atomToName) == 0) ){
      return TRUE;
    }
  }
  return FALSE;
}

BOOL ResidueIntraBond13Check(char *atom1, char *atom2, BondSet* pBondSet){
  for(int i = 0; i < pBondSet->count; i++){
    Bond *pBond = pBondSet->bonds+i;
    if( strcmp(atom1, pBond->atomFromName) == 0){
      if(ResidueIntraBond12Check(pBond->atomToName, atom2, pBondSet)){
        return TRUE;
      }
    }
    else if(strcmp(atom1, pBond->atomToName) == 0){
      if(ResidueIntraBond12Check(pBond->atomFromName, atom2, pBondSet)){
        return TRUE;
      }
    }
  }
  return FALSE;
}

BOOL ResidueIntraBond14Check(char *atom1, char *atom2, BondSet* pBondSet){
  for(int i = 0; i < pBondSet->count; i++){
    Bond *pBond = pBondSet->bonds+i;
    if( strcmp(atom1, pBond->atomFromName) == 0){
      if( ResidueIntraBond13Check(pBond->atomToName, atom2, pBondSet) ){
        return TRUE;
      }
    }
    else if( strcmp(atom1, pBond->atomToName) == 0 ){
      if( ResidueIntraBond13Check(pBond->atomFromName, atom2, pBondSet) ){
        return TRUE;
      }
    }
  }
  return FALSE;
}

int ResidueIntraBondConnectionCheck(char *atom1, char *atom2, BondSet* pBondSet){
  if(ResidueIntraBond12Check(atom1, atom2, pBondSet)){
    return 12;
  }
  else if(ResidueIntraBond13Check(atom1, atom2, pBondSet)){
    return 13;
  }
  else if(ResidueIntraBond14Check(atom1, atom2, pBondSet)){
    return 14;
  }
  else{
    return 15;
  }
}

int ResidueAndNextResidueInterBondConnectionCheck_charmm22(char *atomOnPreResi, char *atomOnNextResi, Residue *pPreResi, Residue *pNextResi){
  if( strcmp(atomOnPreResi, "C") == 0 ){
    if( strcmp(atomOnNextResi, "N") == 0 ) return 12;
    else if( strcmp(atomOnNextResi, "CA") == 0 || strcmp(atomOnNextResi, "HN") == 0 || (strcmp(atomOnNextResi, "CD") == 0 && strcmp(pNextResi->name, "PRO") == 0) ){
      return 13;
    }
    else if( strcmp(atomOnNextResi, "CB") == 0|| strcmp(atomOnNextResi, "C") == 0 || 
      strcmp(atomOnNextResi, "HA") == 0 || strcmp(atomOnNextResi, "HA1") == 0 || strcmp(atomOnNextResi, "HA2") == 0 || 
      (strcmp(atomOnNextResi, "HD1") == 0 && strcmp(pNextResi->name, "PRO") == 0) || 
      (strcmp(atomOnNextResi, "HD2") == 0 && strcmp(pNextResi->name, "PRO") == 0)|| 
      (strcmp(atomOnNextResi, "CG") == 0 && strcmp(pNextResi->name, "PRO") == 0)){
      return 14;
    }
  }
  else if( strcmp(atomOnPreResi, "O") == 0 || strcmp(atomOnPreResi, "CA") == 0 ){
    if(strcmp(atomOnNextResi, "N") == 0) return 13;
    else if( strcmp(atomOnNextResi, "CA") == 0 || strcmp(atomOnNextResi, "HN") == 0|| 
      (strcmp(atomOnNextResi, "CD") == 0 && strcmp(pNextResi->name, "PRO") == 0) ){
      return 14;
    }
  }
  else if( strcmp(atomOnPreResi, "CB") == 0 || strcmp(atomOnPreResi, "HA") == 0 || 
    strcmp(atomOnPreResi, "HA1") == 0 || strcmp(atomOnPreResi, "HA2") == 0 || strcmp(atomOnPreResi, "N") == 0 ){
    if(strcmp(atomOnNextResi, "N") == 0) return 14;
  }
  return 15;
}

int ResidueAndNextResidueInterBondConnectionCheck_charmm19(char *atomOnPreResi, char *atomOnNextResi, char *nextResiName){
  if( strcmp(atomOnPreResi, "C") == 0 ){
    if( strcmp(atomOnNextResi, "N") == 0 ) return 12;
    else if( strcmp(atomOnNextResi, "CA") == 0|| strcmp(atomOnNextResi, "H") == 0|| (strcmp(atomOnNextResi, "CD") == 0 && strcmp(nextResiName, "PRO") == 0) )
      return 13;
    else if( strcmp(atomOnNextResi, "CB") == 0|| strcmp(atomOnNextResi, "C") == 0 || (strcmp(atomOnNextResi, "CG") == 0 && strcmp(nextResiName, "PRO") == 0))
      return 14;
  }
  else if( strcmp(atomOnPreResi, "O") == 0 || strcmp(atomOnPreResi, "CA") == 0 ){
    if( strcmp(atomOnNextResi, "N") == 0 ) return 13;
    else if( strcmp(atomOnNextResi, "CA") == 0|| strcmp(atomOnNextResi, "H") == 0||
      (strcmp(atomOnNextResi, "CD") == 0 && strcmp(nextResiName, "PRO") == 0) ){
      return 14;
    }
  }
  else if( strcmp(atomOnPreResi, "CB") == 0 || strcmp(atomOnPreResi, "N") == 0 ){
    if( strcmp(atomOnNextResi, "N") == 0 ) return 14;
  }
  return 15;
}




#define DEAL_WITH_BASIC_ENERGY_TERMS
//////////////////////////////////////////////////////////////
//Atomic pairwise energy for residue
//////////////////////////////////////////////////////////////
int VdwAttEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance, int bondType, double *vdwAtt){
  if(distance>=ENERGY_DISTANCE_CUTOFF) return Success;
  if(bondType==12||bondType==13) return Success;
  if(AtomIsHydrogen(pAtom1) || AtomIsHydrogen(pAtom2)) return Success;
  double rmin = RADIUS_SCALE_FOR_VDW * (pAtom1->vdw_radius + pAtom2->vdw_radius);
  double ratio = distance/rmin;
  double energy=0.0;
  double scale=0.0;
  if(bondType==14) scale=ENERGY_SCALE_FACTOR_BOND_14;
  else if(bondType==15) scale=ENERGY_SCALE_FACTOR_BOND_15;

  if(ratio < 0.8909){
    energy=0.0;
  }
  else if(distance <= 5.0){
    double epsilon  = sqrt(pAtom1->vdw_epsilon*pAtom2->vdw_epsilon);
    double B6 = pow(1/ratio, 6.0);
    double A12 = B6*B6;
    energy = epsilon * (A12 - 2.0 * B6);
  }
  else if(distance > 5.0 && distance < ENERGY_DISTANCE_CUTOFF){
    double epsilon  = sqrt(pAtom1->vdw_epsilon*pAtom2->vdw_epsilon);
    double B6 = pow((double)rmin/5.0, 6.0);
    double A12 = B6 * B6;
    double M = epsilon * ( A12 - 2.0 * B6);
    double N = 2.4 * epsilon * (B6 - A12);
    double a = 2 * M + N;
    double b = -33 * M - 17 * N;
    double c = 180 * M + 96 * N;
    double d = -324 * M -180 * N;
    energy = a * distance * distance * distance + b * distance * distance + c * distance + d;
  }
  energy*=scale;
  *vdwAtt=energy;
  if(ENERGY_DEBUG_MODE_VDW_ATT){
    printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, bondType: %2d, dist: %5.2f, ratio: %5.2f, vdwAtt: %5.2f\n", 
      AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1), 
      AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
      bondType,distance, ratio, energy);
  }
  return Success;
}



int VdwRepEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance,int bondType, double *vdwRep){
  if(bondType==12||bondType==13) return Success;
  //if(AtomIsHydrogen(pAtom1) || AtomIsHydrogen(pAtom2)) return Success;
  double rmin = RADIUS_SCALE_FOR_VDW * (pAtom1->vdw_radius + pAtom2->vdw_radius);
  double ratio = distance/rmin;
  double epsilon  = sqrt(pAtom1->vdw_epsilon*pAtom2->vdw_epsilon);
  double RATIO_CUTOFF = 0.70; // can be adjusted
  double energy=0.0;
  double scale=0.0;
  if(bondType==14) scale=ENERGY_SCALE_FACTOR_BOND_14;
  else if(bondType==15) scale=ENERGY_SCALE_FACTOR_BOND_15;

  if(ratio > 0.8909){ // 0.8909 ~ inf
    energy=0.0;
  }
  else if(ratio >= RATIO_CUTOFF){ // [0.70, 0.8909]
    double B6 = pow(1/ratio, 6.0);
    double A12 = B6*B6;
    energy = epsilon * (A12 - 2.0 * B6);
  }
  else{ // [0, 0.70]
    double B6_0 = pow(1/RATIO_CUTOFF, 6.0);
    double a = epsilon * (B6_0 * B6_0 - 2.0 * B6_0);
    double b = epsilon * 12.0 * (B6_0 / RATIO_CUTOFF - B6_0 * B6_0 / RATIO_CUTOFF);
    double y0 = a * epsilon;
    double k = b * epsilon;
    energy = k * (ratio - RATIO_CUTOFF) + y0;
  }
  energy*=scale;
  //set a cutoff for maximum clash
  //disable the cap to achieve better design results using EvoEF2 weights(../wread/weight_EvoEF2.txt)
  //double MAX_CLASH=5.0*epsilon;
  //if(energy>MAX_CLASH) energy=MAX_CLASH;
  *vdwRep=energy;

  if(ENERGY_DEBUG_MODE_VDW_REP){
    printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, bondType: %2d, dist: %5.2f, ratio: %5.2f, vdwRep: %5.2f\n", 
      AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1),
      AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
      bondType,distance, ratio, energy);
  }
  return Success;
}


int HBondEnergyAtomAndAtom(Atom *atomH, Atom *atomA, Atom*atomD, Atom *atomB,double distanceHA, int bondType, double *etotal, double *edist, double *etheta, double *ephi){
  if(bondType==12||bondType==13) return Success;
  if(distanceHA > HBOND_DISTANCE_CUTOFF_MAX) return Success;
  XYZ xyzDH = XYZDifference(&atomD->xyz, &atomH->xyz);
  XYZ xyzHA = XYZDifference(&atomH->xyz, &atomA->xyz);
  XYZ xyzAAB = XYZDifference(&atomA->xyz, &atomB->xyz);
  //if(fabs(XYZAngle(&xyzDH, &xyzHA)-1000.0)<1.0){
  //  printf("atomD: %s, atomH: %s, atomA: %s\n", AtomGetName(atomD),AtomGetName(atomH),AtomGetName(atomA));
  //}
  //if(fabs(XYZAngle(&xyzHA, &xyzAAB)-1000.0)<1.0){
  //  printf("atomH: %s, atomA: %s, atomAB: %s\n", AtomGetName(atomH),AtomGetName(atomA),AtomGetName(atomB));
  //}
  double angleTheta = PI-XYZAngle(&xyzDH, &xyzHA);
  if(RadToDeg(angleTheta)<90) return Success;
  double anglePhi   = PI-XYZAngle(&xyzHA, &xyzAAB);
  if(RadToDeg(anglePhi)<80) return Success;

  double energyR=0.0;
  if(distanceHA<HBOND_OPTIMAL_DISTANCE){
    energyR=-1.0*HBOND_WELL_DEPTH*cos((distanceHA-HBOND_OPTIMAL_DISTANCE)*PI);
  }
  else{
    energyR=-0.5*cos(PI/(HBOND_DISTANCE_CUTOFF_MAX-HBOND_OPTIMAL_DISTANCE)*(distanceHA-HBOND_OPTIMAL_DISTANCE))-0.5;
  }
  if(energyR > 0.0) energyR = 0.0;

  double energyTheta = -1.0*cos(angleTheta)*cos(angleTheta)*cos(angleTheta)*cos(angleTheta);
  double energyPhi = 0.0;
  if(atomH->isBBAtom==TRUE && atomA->isBBAtom==TRUE){
    energyPhi=-1.0*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150));
  }
  else{
    if(atomA->hybridType == Type_AtomHybridType_SP3){
      energyPhi = -1.0*cos(anglePhi-DegToRad(135))*cos(anglePhi-DegToRad(135))*cos(anglePhi-DegToRad(135))*cos(anglePhi-DegToRad(135));
    }
    else if(atomA->hybridType == Type_AtomHybridType_SP2){
      energyPhi = -1.0*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150));
    }
  }
  
  // original function is error, we should add angle and restriction to calculate energy
  double energy = 0.0;
  if(RadToDeg(angleTheta) >= 90 && RadToDeg(anglePhi) >= 80 && distanceHA < HBOND_DISTANCE_CUTOFF_MAX){
    energy = energyR + energyTheta + energyPhi;
    if(energy>0.0) energy=0.0;
    *etotal = energy;
    *edist = energyR;
    *etheta = energyTheta;
    *ephi = energyPhi;
  }

  if(ENERGY_DEBUG_MODE_HBOND){
    printf("AtomH: %1s %4d %4s, AtomA: %1s %4d %4s, dist: %5.2f, theta: %5.1f, phi: %5.1f, edist: %5.2f, ethe: %5.2f, ephi: %5.2f\n", 
      AtomGetChainName(atomH), AtomGetPosInChain(atomH), AtomGetName(atomH),
      AtomGetChainName(atomA), AtomGetPosInChain(atomA), AtomGetName(atomA),
      distanceHA, RadToDeg(angleTheta),RadToDeg(anglePhi),energyR,energyTheta,energyPhi);
  }
  return Success;
}


int ElecEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance12,int bondType, double *elec){
  if(bondType==12||bondType==13) return Success;
  if(distance12 > ELEC_DISTANCE_CUTOFF) return Success;
  if(fabs(pAtom1->charge)<1e-2 || fabs(pAtom2->charge)<1e-2) return Success;
  else if(distance12<0.8*(pAtom1->vdw_radius + pAtom2->vdw_radius)) distance12=0.8*(pAtom1->vdw_radius + pAtom2->vdw_radius);

  double energy = COULOMB_CONSTANT*pAtom1->charge*pAtom2->charge/distance12/distance12/40.0;
  double scale=0.0;
  if(bondType==14) scale=ENERGY_SCALE_FACTOR_BOND_14;
  else if(bondType==15) scale=ENERGY_SCALE_FACTOR_BOND_15;
  energy*=scale;
  *elec = energy;

  if(ENERGY_DEBUG_MODE_ELEC){
    if(fabs(energy) > 0.0 && distance12 < 4.0){
      printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, bondType: %2d, dist: %5.2f, elec: %5.2f\n", 
        AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1), 
        AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
        bondType, distance12, energy);
    }
  }
  return Success;
}

int LKDesolvationEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance,int bondType, double *energyP, double *energyH){
  if(bondType==12||bondType==13) return Success;
  if(AtomIsHydrogen(pAtom1) || AtomIsHydrogen(pAtom2)) return Success;
  if(distance>ENERGY_DISTANCE_CUTOFF) return Success;
  double volume1 = pAtom1->EEF1_volume;
  double volume2 = pAtom2->EEF1_volume;
  double dGFreeAtom1 = pAtom1->EEF1_freeDG;
  double dGFreeAtom2 = pAtom2->EEF1_freeDG;
  double coefficient = -0.089793561062582974; // 0.5/(pi^1.5)
  double r1=pAtom1->vdw_radius*RADIUS_SCALE_FOR_DESOLV;
  double r2=pAtom2->vdw_radius*RADIUS_SCALE_FOR_DESOLV;
  double r12 = r1+r2;

  distance = distance < r12 ? r12 : distance;
  double lamda1 = pAtom1->EEF1_lamda_ * distance * distance;
  double lamda2 = pAtom2->EEF1_lamda_ * distance * distance;
  double x1 = (distance - r1)/pAtom1->EEF1_lamda_;
  double x2 = (distance - r2)/pAtom2->EEF1_lamda_;

  double desolv12 = coefficient * volume2 * dGFreeAtom1 / lamda1;
  desolv12 *= exp( -1.0 * x1 * x1 );
  double desolv21 = coefficient * volume1 * dGFreeAtom2 / lamda2;
  desolv21 *= exp( -1.0 * x2 * x2 );
  if(pAtom1->polarity == Type_AtomPolarity_P || pAtom1->polarity == Type_AtomPolarity_C) *energyP += desolv12;
  else *energyH += desolv12;
  if(pAtom2->polarity == Type_AtomPolarity_P || pAtom2->polarity == Type_AtomPolarity_C) *energyP += desolv21;
  else *energyH += desolv21;

  if(ENERGY_DEBUG_MODE_DESOLV){
    printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, dist: %5.2f, desolv12: %7.3f, desolv21: %7.3f\n",
      AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1), 
      AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
      distance, desolv12, desolv21);
  }

  return Success;
}

///////////////////////////////////////////////////////////////////
//User can define new energy term here
///////////////////////////////////////////////////////////////////

int SSbondEnergyAtomAndAtom(Atom *pAtomS1,Atom *pAtomS2,Atom *pAtomCB1,Atom *pAtomCB2,Atom* pAtomCA1,Atom* pAtomCA2,double* sse){
  XYZ xyzC1S1=XYZDifference(&pAtomCB1->xyz,&pAtomS1->xyz);
  XYZ xyzS1S2=XYZDifference(&pAtomS1->xyz,&pAtomS2->xyz);
  XYZ xyzS2C2=XYZDifference(&pAtomS2->xyz,&pAtomCB2->xyz);
  double dSS=XYZDistance(&pAtomS1->xyz,&pAtomS2->xyz);
  double aC1S1S2=RadToDeg(PI-XYZAngle(&xyzS1S2,&xyzS2C2));
  double aC2S2S1=RadToDeg(PI-XYZAngle(&xyzC1S1,&xyzS1S2));
  double xC1S1S2C2=GetTorsionAngle(&pAtomCB1->xyz,&pAtomS1->xyz,&pAtomS2->xyz,&pAtomCB2->xyz);
  double xCA1CB1SG1SG2=GetTorsionAngle(&pAtomCA1->xyz,&pAtomCB1->xyz,&pAtomS1->xyz,&pAtomS2->xyz);
  double xCA2CB2SG2SG1=GetTorsionAngle(&pAtomCA2->xyz,&pAtomCB2->xyz,&pAtomS2->xyz,&pAtomS1->xyz);
  *sse=(0.8*(1-exp(-10.0*(dSS-SSBOND_DISTANCE)))*(1-exp(-10.0*(dSS-SSBOND_DISTANCE))))
    +0.005*(aC1S1S2-SSBOND_ANGLE)*(aC1S1S2-SSBOND_ANGLE)
    +0.005*(aC2S2S1-SSBOND_ANGLE)*(aC2S2S1-SSBOND_ANGLE)
    +cos(2.0*xC1S1S2C2)+1.0
    +1.25*sin(xCA1CB1SG1SG2+2.0*PI/3.0)-1.75
    +1.25*sin(xCA2CB2SG2SG1+2.0*PI/3.0)-1.75;
  if(*sse>0.0) *sse=0.0;
  return Success;
}



#define RESIDUE_PAIRWISE_ENERGY_COMPUTATION
int EVOEF_AminoAcidReferenceEnergy(char *AAname, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]){
  if(strcmp(AAname, "ALA")  == 0)      energyTerm[ 1] += 1.0;
  else if(strcmp(AAname, "CYS")  == 0) energyTerm[ 2] += 1.0;
  else if(strcmp(AAname, "ASP")  == 0) energyTerm[ 3] += 1.0;
  else if(strcmp(AAname, "GLU")  == 0) energyTerm[ 4] += 1.0;
  else if(strcmp(AAname, "PHE")  == 0) energyTerm[ 5] += 1.0;
  else if(strcmp(AAname, "GLY")  == 0) energyTerm[ 6] += 1.0;
  else if(strcmp(AAname, "HIS")  == 0) energyTerm[ 7] += 1.0;
  else if(strcmp(AAname, "HSE")  == 0) energyTerm[ 7] += 1.0;
  else if(strcmp(AAname, "HSD")  == 0) energyTerm[ 7] += 1.0;
  else if(strcmp(AAname, "ILE")  == 0) energyTerm[ 8] += 1.0;
  else if(strcmp(AAname, "LYS")  == 0) energyTerm[ 9] += 1.0;
  else if(strcmp(AAname, "LEU")  == 0) energyTerm[10] += 1.0;
  else if(strcmp(AAname, "MET")  == 0) energyTerm[11] += 1.0;
  else if(strcmp(AAname, "ASN")  == 0) energyTerm[12] += 1.0;
  else if(strcmp(AAname, "PRO")  == 0) energyTerm[13] += 1.0;
  else if(strcmp(AAname, "GLN")  == 0) energyTerm[14] += 1.0;
  else if(strcmp(AAname, "ARG")  == 0) energyTerm[15] += 1.0;
  else if(strcmp(AAname, "SER")  == 0) energyTerm[16] += 1.0;
  else if(strcmp(AAname, "THR")  == 0) energyTerm[17] += 1.0;
  else if(strcmp(AAname, "VAL")  == 0) energyTerm[18] += 1.0;
  else if(strcmp(AAname, "TRP")  == 0) energyTerm[19] += 1.0;
  else if(strcmp(AAname, "TYR")  == 0) energyTerm[20] += 1.0;
  return Success;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// the following are energy between residue and residue
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int EVOEF_EnergyResidueIntraEnergy(Residue* pThis,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<ResidueGetAtomCount(pThis); ++i){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    if(pAtom1->isXyzValid==FALSE)continue;
    for(int j=i+1; j<ResidueGetAtomCount(pThis);++j){
      Atom* pAtom2=ResidueGetAtom(pThis,j);
      if(pAtom2->isXyzValid==FALSE)continue;
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      //skip mainchain-mainchain
      if(pAtom2->isBBAtom==TRUE && pAtom1->isBBAtom==TRUE) continue;
      if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){
        if(strcmp(ResidueGetName(pThis),"ILE")==0 || strcmp(ResidueGetName(pThis),"MET")==0 ||
          strcmp(ResidueGetName(pThis),"GLN")==0||strcmp(ResidueGetName(pThis),"GLU")==0||
          strcmp(ResidueGetName(pThis),"LYS")==0||strcmp(ResidueGetName(pThis),"ARG")==0){
            int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtom1),AtomGetName(pAtom2),ResidueGetBonds(pThis));
            if(bondType==12||bondType==13) continue;
            double vdwAtt=0,vdwRep=0,desolvP=0,desolvH=0;
            VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
            VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
            LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
            energyTerms[21]+=vdwAtt;
            energyTerms[22]+=vdwRep;
            energyTerms[24]+=desolvP;
            energyTerms[25]+=desolvH;
        }
      }
      else if((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE) || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE)){
        if(strcmp(AtomGetName(pAtom1),"CB")==0 || strcmp(AtomGetName(pAtom2),"CB")==0) continue;
        int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtom1),AtomGetName(pAtom2),ResidueGetBonds(pThis));
        if(bondType==12||bondType==13) continue;
        double vdwAtt=0,vdwRep=0,desolvP=0,desolvH=0,ele=0;
        VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
        VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
        LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
        ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
        energyTerms[21]+=vdwAtt;
        energyTerms[22]+=vdwRep;
        energyTerms[23]+=ele;
        energyTerms[24]+=desolvP;
        energyTerms[25]+=desolvH;
        if(distance < HBOND_DISTANCE_CUTOFF_MAX){
          double hb_tot=0,hb_dist=0,hb_theta=0,hb_phi=0;
          if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
            HBondEnergyAtomAndAtom(pAtom1,pAtom2,ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom2)),
              distance,bondType,&hb_tot,&hb_dist,&hb_theta,&hb_phi);
          }
          else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
            HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom2)),ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
              distance,bondType,&hb_tot,&hb_dist,&hb_theta,&hb_phi);
          }
          energyTerms[26]+=hb_dist;
          energyTerms[27]+=hb_theta;
          energyTerms[28]+=hb_phi;
        }
      }
    }

  }  
  return Success;
}

int EVOEF_EnergyResidueAndNextResidue(Residue* pThis, Residue* pOther, double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<ResidueGetAtomCount(pThis); ++i){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    if(pAtom1->isXyzValid==FALSE)continue;
    for(int j=0; j<ResidueGetAtomCount(pOther); ++j){
      Atom* pAtom2=ResidueGetAtom(pOther,j);
      if(pAtom2->isXyzValid==FALSE)continue;
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      //skip the mainchain-mainchain
      if(pAtom2->isBBAtom==TRUE && pAtom1->isBBAtom==TRUE) continue;
      if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){
        int bondType=15;
        double vdwAtt=0,vdwRep=0,desolvP=0,desolvH=0,ele=0;
        VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
        VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
        LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
        ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
        energyTerms[31]+=vdwAtt;
        energyTerms[32]+=vdwRep;
        energyTerms[33]+=ele;
        energyTerms[34]+=desolvP;
        energyTerms[35]+=desolvH;
        if(distance < HBOND_DISTANCE_CUTOFF_MAX){
          double hb_tot=0,hbd=0,hbt=0,hbp=0;
          if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
            HBondEnergyAtomAndAtom(pAtom1,pAtom2,ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
              distance,bondType,&hb_tot,&hbd,&hbt,&hbp);
          }
          else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
            HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
              distance,bondType,&hb_tot,&hbd,&hbt,&hbp);
          }
          if(pThis->posInChain-pOther->posInChain<=2 && pThis->posInChain-pOther->posInChain>=-2){
            hbd *= HBOND_LOCAL_REDUCE;
            hbt *= HBOND_LOCAL_REDUCE;
            hbp *= HBOND_LOCAL_REDUCE;
          }
          if(pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE){
            energyTerms[41]+=hbd;
            energyTerms[42]+=hbt;
            energyTerms[43]+=hbp;
          }
          else if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
            energyTerms[47]+=hbd;
            energyTerms[48]+=hbt;
            energyTerms[49]+=hbp;
          }
          else{
            energyTerms[44]+=hbd;
            energyTerms[45]+=hbt;
            energyTerms[46]+=hbp;
          }
        }
      }
      else if((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE) || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE)){
        int bondType=ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1),AtomGetName(pAtom2),pOther->name);
        if(bondType==12||bondType==13)continue;
        double vdwAtt=0,vdwRep=0,desolvP=0,desolvH=0,ele=0;
        VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
        VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
        LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
        ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
        energyTerms[31]+=vdwAtt;
        energyTerms[32]+=vdwRep;
        energyTerms[33]+=ele;
        energyTerms[34]+=desolvP;
        energyTerms[35]+=desolvH;
        if(distance < HBOND_DISTANCE_CUTOFF_MAX){
          double hb_tot=0,hb_dist=0,hb_theta=0,hb_phi=0;
          if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
            HBondEnergyAtomAndAtom(pAtom1,pAtom2,ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
              distance,bondType,&hb_tot,&hb_dist,&hb_theta,&hb_phi);
          }
          else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
            HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
              distance,bondType,&hb_tot,&hb_dist,&hb_theta,&hb_phi);
          }
          if(pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE){
            energyTerms[41]+=hb_dist;
            energyTerms[42]+=hb_theta;
            energyTerms[43]+=hb_phi;
          }
          else if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
            energyTerms[47]+=hb_dist;
            energyTerms[48]+=hb_theta;
            energyTerms[49]+=hb_phi;
          }
          else{
            energyTerms[44]+=hb_dist;
            energyTerms[45]+=hb_theta;
            energyTerms[46]+=hb_phi;
          }
        }
      }
    }
  }
  return Success;
}

int EVOEF_EnergyResidueAndOtherResidueSameChain(Residue* pThis, Residue* pOther,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<ResidueGetAtomCount(pThis); ++i){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    if(pAtom1->isXyzValid==FALSE)continue;
    for(int j=0; j<ResidueGetAtomCount(pOther); ++j){
      Atom* pAtom2=ResidueGetAtom(pOther,j);
      if(pAtom2->isXyzValid==FALSE)continue;
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      int bondType=15;
      double vdwAtt=0,vdwRep=0,ele=0,desolvP=0,desolvH=0;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
      energyTerms[31]+=vdwAtt;
      energyTerms[32]+=vdwRep;
      energyTerms[33]+=ele;
      energyTerms[34]+=desolvP;
      energyTerms[35]+=desolvH;
      if(distance < HBOND_DISTANCE_CUTOFF_MAX){
        double hb_tot=0,hbd=0,hbt=0,hbp=0;
        if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hb_tot,&hbd,&hbt,&hbp);
        }
        else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hb_tot,&hbd,&hbt,&hbp);
        }
        if(pThis->posInChain-pOther->posInChain<=2 && pThis->posInChain-pOther->posInChain>=-2){
          hbd *= HBOND_LOCAL_REDUCE;
          hbt *= HBOND_LOCAL_REDUCE;
          hbp *= HBOND_LOCAL_REDUCE;
        }
        if(pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE){
          energyTerms[41]+=hbd;
          energyTerms[42]+=hbt;
          energyTerms[43]+=hbp;
        }
        else if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
          energyTerms[47]+=hbd;
          energyTerms[48]+=hbt;
          energyTerms[49]+=hbp;
        }
        else{
          energyTerms[44]+=hbd;
          energyTerms[45]+=hbt;
          energyTerms[46]+=hbp;
        }
      }
    }
  }

  //consider the SSBOND
  if(strcmp("CYS",ResidueGetName(pThis))==0 && strcmp("CYS",ResidueGetName(pOther))==0){
    Atom* pSG1=NULL,*pSG2=NULL,*pCB1=NULL,*pCB2=NULL,*pCA1=NULL,*pCA2=NULL;
    for(int i=0;i<ResidueGetAtomCount(pThis);i++){
      Atom* pAtom=ResidueGetAtom(pThis,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA1=pAtom;
    }
    for(int i=0;i<ResidueGetAtomCount(pOther);i++){
      Atom* pAtom=ResidueGetAtom(pOther,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA2=pAtom;
    }
    if(pSG1!=NULL && pCB1!=NULL &&pSG2!=NULL && pCB2!=NULL && pCA1!=NULL && pCA2!=NULL){
      double dist=XYZDistance(&pSG1->xyz,&pSG2->xyz);
      if(dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN){
        double ssbond=0;
        SSbondEnergyAtomAndAtom(pSG1,pSG2,pCB1,pCB2,pCA1,pCA2,&ssbond);
        energyTerms[36]+=ssbond;
      }
    }
  }

  return Success;
}


int EVOEF_EnergyResidueAndOtherResidueDifferentChain(Residue* pThis, Residue* pOther,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<ResidueGetAtomCount(pThis); ++i){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    if(pAtom1->isXyzValid==FALSE)continue;
    for(int j=0; j<ResidueGetAtomCount(pOther); ++j){
      Atom* pAtom2=ResidueGetAtom(pOther,j);
      if(pAtom2->isXyzValid==FALSE)continue;
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      int bondType=15;
      double vdwAtt=0,vdwRep=0,ele=0,desolvP=0,desolvH=0;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
      energyTerms[51]+=vdwAtt;
      energyTerms[52]+=vdwRep;
      energyTerms[53]+=ele;
      energyTerms[54]+=desolvP;
      energyTerms[55]+=desolvH;
      if(distance < HBOND_DISTANCE_CUTOFF_MAX){
        double hb_tot=0,hb_dist=0,hb_theta=0,hb_phi=0;
        if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hb_tot,&hb_dist,&hb_theta,&hb_phi);
        }
        else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hb_tot,&hb_dist,&hb_theta,&hb_phi);
        }
        if(pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE){
          energyTerms[61]+=hb_dist;
          energyTerms[62]+=hb_theta;
          energyTerms[63]+=hb_phi;
        }
        else if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
          energyTerms[67]+=hb_dist;
          energyTerms[68]+=hb_theta;
          energyTerms[69]+=hb_phi;
        }
        else{
          energyTerms[64]+=hb_dist;
          energyTerms[65]+=hb_theta;
          energyTerms[66]+=hb_phi;
        }
      }
    }
  }

  //consider the SSBOND
  if(strcmp("CYS",ResidueGetName(pThis))==0 && strcmp("CYS",ResidueGetName(pOther))==0){
    Atom* pSG1=NULL,*pSG2=NULL,*pCB1=NULL,*pCB2=NULL,*pCA1=NULL,*pCA2=NULL;
    for(int i=0;i<ResidueGetAtomCount(pThis);i++){
      Atom* pAtom=ResidueGetAtom(pThis,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA1=pAtom;
    }
    for(int i=0;i<ResidueGetAtomCount(pOther);i++){
      Atom* pAtom=ResidueGetAtom(pOther,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA2=pAtom;
    }
    if(pSG1!=NULL && pCB1!=NULL &&pSG2!=NULL && pCB2!=NULL &&pCA1!=NULL && pCA2!=NULL){
      double dist=XYZDistance(&pSG1->xyz,&pSG2->xyz);
      if(dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN){
        double ssbond=0;
        SSbondEnergyAtomAndAtom(pSG1,pSG2,pCB1,pCB2,pCA1,pCA2,&ssbond);
        energyTerms[56]+=ssbond;
      }
    }
  }

  return Success;
}


int EVOEF_EnergyResidueAndLigandResidue(Residue* pProtein, Residue* pLigand,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<ResidueGetAtomCount(pProtein); ++i){
    Atom* pAtom1=ResidueGetAtom(pProtein,i);
    if(pAtom1->isXyzValid==FALSE)continue;
    for(int j=0; j<ResidueGetAtomCount(pLigand); ++j){
      Atom* pAtom2=ResidueGetAtom(pLigand,j);
      if(pAtom2->isXyzValid==FALSE)continue;
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      int bondType=15;
      double vdwAtt=0,vdwRep=0,ele=0,desolvP=0,desolvH=0;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
      energyTerms[71]+=vdwAtt;
      energyTerms[72]+=vdwRep;
      energyTerms[73]+=ele;
      energyTerms[74]+=desolvP;
      energyTerms[75]+=desolvH;
      if(distance < HBOND_DISTANCE_CUTOFF_MAX){
        double hb_tot=0,hb_dist=0,hb_theta=0,hb_phi=0;
        if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,ResidueGetAtomByName(pProtein,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pLigand,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hb_tot,&hb_dist,&hb_theta,&hb_phi);
        }
        else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pLigand,AtomGetHbDorB(pAtom2)),ResidueGetAtomByName(pProtein,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hb_tot,&hb_dist,&hb_theta,&hb_phi);
        }
        if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
          energyTerms[84]+=hb_dist;
          energyTerms[85]+=hb_theta;
          energyTerms[86]+=hb_phi;
        }
        else{
          energyTerms[81]+=hb_dist;
          energyTerms[82]+=hb_theta;
          energyTerms[83]+=hb_phi;
        }
      }
    }
  }
  return Success;
}




#define ROTAMER_PAIRWISE_ENERGY_COMPUTATION
////////////////////////////////////////////////////////////////////////////////////////////
//these functions are used to calculate energy for rotamers
////////////////////////////////////////////////////////////////////////////////////////////

int EVOEF_EnergyProteinRotamerIntraEnergy(Rotamer* pThis,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
	for(int i=0; i<RotamerGetAtomCount(pThis); ++i){
		Atom* pAtom1=RotamerGetAtom(pThis,i);
		for(int j=i+1; j<RotamerGetAtomCount(pThis);++j){
			Atom* pAtom2=RotamerGetAtom(pThis,j);
			//skip mainchain-mainchain
			if(pAtom2->isBBAtom==TRUE && pAtom1->isBBAtom==TRUE) continue;
			double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
			if(distance>ENERGY_DISTANCE_CUTOFF) continue;
			if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){//two atoms in the sidechain
				if(strcmp(RotamerGetType(pThis),"ILE")==0 || strcmp(RotamerGetType(pThis),"MET")==0 ||
					strcmp(RotamerGetType(pThis),"GLN")==0||strcmp(RotamerGetType(pThis),"GLU")==0||
					strcmp(RotamerGetType(pThis),"LYS")==0||strcmp(RotamerGetType(pThis),"ARG")==0){
						int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtom1),AtomGetName(pAtom2),RotamerGetBonds(pThis));
						if(bondType==12||bondType==13) continue;
						double vdwAtt=0,vdwRep=0,desolvP=0,desolvH=0;
						VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
						VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
						LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
						energyTerms[21]+=vdwAtt;
						energyTerms[22]+=vdwRep;
						energyTerms[24]+=desolvP;
						energyTerms[25]+=desolvH;
				}
			}
			else if((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE) || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE)){
				//one atom in the mainchain, the other in the sidechain
				if(strcmp(AtomGetName(pAtom1),"CB")==0 || strcmp(AtomGetName(pAtom2),"CB")==0) continue;
				int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtom1),AtomGetName(pAtom2),RotamerGetBonds(pThis));
				if(bondType==12||bondType==13) continue;
				double vdwAtt=0,vdwRep=0,desolvP=0,desolvH=0,ele=0;
				VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
				VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
				LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
				ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
				energyTerms[21]+=vdwAtt;
				energyTerms[22]+=vdwRep;
				energyTerms[23]+=ele;
				energyTerms[24]+=desolvP;
				energyTerms[25]+=desolvH;
				if(distance < HBOND_DISTANCE_CUTOFF_MAX){
					double hbtot=0,hbd=0,hbt=0,hbp=0;
					if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
						HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom2)),
							distance,bondType,&hbtot,&hbd,&hbt,&hbp);
						//only hbscbb
						energyTerms[26]+=hbd;
					  energyTerms[27]+=hbt;
						energyTerms[28]+=hbp;
					}
					else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
						HBondEnergyAtomAndAtom(pAtom2,pAtom1,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
							distance,bondType,&hbtot,&hbd,&hbt,&hbp);
						//only hbscbb
						energyTerms[26]+=hbd;
					  energyTerms[27]+=hbt;
						energyTerms[28]+=hbp;
					}

				}
			}
		}
	}
	return Success;
}


int EVOEF_EnergyProteinRotamerAndRotamerSameChain(Rotamer* pThis, Rotamer* pOther, double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  //only consider sidechain-sidechain interactions
  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isBBAtom) continue;
    for(int j=0; j <RotamerGetAtomCount(pOther); j++){
      Atom* pAtom2 = RotamerGetAtom(pOther, j);
      if(pAtom2->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0, rep=0, ele=0, desP=0, desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[31] += att;
      energyTerms[32] += rep;
      energyTerms[33] += ele;
      energyTerms[34] += desP;
      energyTerms[35] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),RotamerGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          if(pThis->posInChain-pOther->posInChain<=2 && pThis->posInChain-pOther->posInChain>=-2){
            hbd *= HBOND_LOCAL_REDUCE;
            hbt *= HBOND_LOCAL_REDUCE;
            hbp *= HBOND_LOCAL_REDUCE;
          }
          //only hbscsc
          energyTerms[47] += hbd;
          energyTerms[48] += hbt;
          energyTerms[49] += hbp;
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,RotamerGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          if(pThis->posInChain-pOther->posInChain<=2 && pThis->posInChain-pOther->posInChain>=-2){
            hbd *= HBOND_LOCAL_REDUCE;
            hbt *= HBOND_LOCAL_REDUCE;
            hbp *= HBOND_LOCAL_REDUCE;
          }
          //only hbscsc
          energyTerms[47] += hbd;
          energyTerms[48] += hbt;
          energyTerms[49] += hbp;
        }
        
      }
    }
  }

  //consider the SSBOND
  if(strcmp("CYS",RotamerGetType(pThis))==0 && strcmp("CYS",RotamerGetType(pOther))==0){
    Atom* pSG1=NULL,*pSG2=NULL,*pCB1=NULL,*pCB2=NULL,*pCA1=NULL,*pCA2=NULL;
    for(int i=0;i<RotamerGetAtomCount(pThis);i++){
      Atom* pAtom=RotamerGetAtom(pThis,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA1=pAtom;
    }
    for(int i=0;i<RotamerGetAtomCount(pOther);i++){
      Atom* pAtom=RotamerGetAtom(pOther,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA2=pAtom;
    }
    if(pSG1!=NULL && pCB1!=NULL &&pSG2!=NULL && pCB2!=NULL && pCA1!=NULL && pCA2!=NULL){
      double dist=XYZDistance(&pSG1->xyz,&pSG2->xyz);
      if(dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN){
        double ssbond=0;
        SSbondEnergyAtomAndAtom(pSG1,pSG2,pCB1,pCB2,pCA1,pCA2,&ssbond);
        energyTerms[36]+=ssbond;
      }
    }
  }


  return Success;
}

int EVOEF_EnergyProteinRotamerAndRotamerDifferentChain(Rotamer* pThis, Rotamer* pOther, double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isBBAtom) continue;
    for(int j=0; j <RotamerGetAtomCount(pOther); j++){
      Atom* pAtom2 = RotamerGetAtom(pOther, j);
      if(pAtom2->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0, rep=0, ele=0, desP=0, desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[51] += att;
      energyTerms[52] += rep;
      energyTerms[53] += ele;
      energyTerms[54] += desP;
      energyTerms[55] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),RotamerGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          //only hbscsc
          energyTerms[67] += hbd;
          energyTerms[68] += hbt;
          energyTerms[69] += hbp;
        }
        else 
        if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,RotamerGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          //only hbscsc
          energyTerms[67] += hbd;
          energyTerms[68] += hbt;
          energyTerms[69] += hbp;
        }
      }
    }
  }

  //consider the SSBOND
  if(strcmp("CYS",RotamerGetType(pThis))==0 && strcmp("CYS",RotamerGetType(pOther))==0){
    Atom* pSG1=NULL,*pSG2=NULL,*pCB1=NULL,*pCB2=NULL,*pCA1=NULL,*pCA2=NULL;
    for(int i=0;i<RotamerGetAtomCount(pThis);i++){
      Atom* pAtom=RotamerGetAtom(pThis,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA1=pAtom;
    }
    for(int i=0;i<RotamerGetAtomCount(pOther);i++){
      Atom* pAtom=RotamerGetAtom(pOther,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA2=pAtom;
    }
    if(pSG1!=NULL && pCB1!=NULL &&pSG2!=NULL && pCB2!=NULL && pCA1!=NULL && pCA2!=NULL){
      double dist=XYZDistance(&pSG1->xyz,&pSG2->xyz);
      if(dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN){
        double ssbond=0;
        SSbondEnergyAtomAndAtom(pSG1,pSG2,pCB1,pCB2,pCA1,pCA2,&ssbond);
        energyTerms[56]+=ssbond;
      }
    }
  }


  return Success;
}



int EVOEF_EnergyLigandRotamerAndRotamer(Rotamer* pThis, Rotamer* pOther, double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isBBAtom) continue;
    for(int j=0; j <RotamerGetAtomCount(pOther); j++){
      Atom* pAtom2 = RotamerGetAtom(pOther, j);
      if(pAtom2->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0, rep=0, ele=0, desP=0, desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[71] += att;
      energyTerms[72] += rep;
      energyTerms[73] += ele;
      energyTerms[74] += desP;
      energyTerms[75] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),RotamerGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          //only hbscsc
          energyTerms[84] += hbd;
          energyTerms[85] += hbt;
          energyTerms[86] += hbp;
        }
        else{
          if(pAtom1->isHBatomA && pAtom2->isHBatomH){
            HBondEnergyAtomAndAtom(pAtom2,pAtom1,RotamerGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
              distance,15,&hbtot,&hbd,&hbt,&hbp);
            //only hbscsc
            energyTerms[84] += hbd;
            energyTerms[85] += hbt;
            energyTerms[86] += hbp;
          }
        }
      }
    }
  }

  return Success;
}



int EVOEF_EnergyProteinRotamerAndFixedResidueSameChain(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  int neighborCheck = 0;
  if(RotamerGetPosInChain(pThis)+1==ResidueGetPosInChain(pOther)) neighborCheck=12;
  else if(RotamerGetPosInChain(pThis)-1==ResidueGetPosInChain(pOther)) neighborCheck=21;

  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isBBAtom) continue;
    for(int j=0; j <ResidueGetAtomCount(pOther); j++){
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      int bondType=15;
      if(neighborCheck==12) bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1),AtomGetName(pAtom2),ResidueGetName(pOther));
      else if(neighborCheck==21) bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom2),AtomGetName(pAtom1),RotamerGetType(pThis));
      double att=0, rep=0, ele=0, desP=0, desH=0;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[31] += att;
      energyTerms[32] += rep;
      energyTerms[33] += ele;
      energyTerms[34] += desP;
      energyTerms[35] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          if(pThis->posInChain-pOther->posInChain<=2 && pThis->posInChain-pOther->posInChain>=-2){
            hbd *= HBOND_LOCAL_REDUCE;
            hbt *= HBOND_LOCAL_REDUCE;
            hbp *= HBOND_LOCAL_REDUCE;
          }
          //only hbscbb or hbscsc
          if(pAtom2->isBBAtom == TRUE){
            energyTerms[44] += hbd;
            energyTerms[45] += hbt;
            energyTerms[46] += hbp;
          }
          else{
            energyTerms[47] += hbd;
            energyTerms[48] += hbt;
            energyTerms[49] += hbp;
          }
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          if(pThis->posInChain-pOther->posInChain<=2 && pThis->posInChain-pOther->posInChain>=-2){
            hbd *= HBOND_LOCAL_REDUCE;
            hbt *= HBOND_LOCAL_REDUCE;
            hbp *= HBOND_LOCAL_REDUCE;
          }
          //only hbscbb or hbscsc
          if(pAtom2->isBBAtom == TRUE){
            energyTerms[44] += hbd;
            energyTerms[45] += hbt;
            energyTerms[46] += hbp;
          }
          else{
            energyTerms[47] += hbd;
            energyTerms[48] += hbt;
            energyTerms[49] += hbp;
          }
        }
        
      }
    }
  }

  //consider the SSBOND
  if(strcmp("CYS",RotamerGetType(pThis))==0 && strcmp("CYS",ResidueGetName(pOther))==0){
    Atom* pSG1=NULL,*pSG2=NULL,*pCB1=NULL,*pCB2=NULL,*pCA1=NULL,*pCA2=NULL;
    for(int i=0;i<RotamerGetAtomCount(pThis);i++){
      Atom* pAtom=RotamerGetAtom(pThis,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA1=pAtom;
    }
    for(int i=0;i<ResidueGetAtomCount(pOther);i++){
      Atom* pAtom=ResidueGetAtom(pOther,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA2=pAtom;
    }
    if(pSG1!=NULL && pCB1!=NULL &&pSG2!=NULL && pCB2!=NULL && pCA1!=NULL && pCA2!=NULL){
      double dist=XYZDistance(&pSG1->xyz,&pSG2->xyz);
      if(dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN){
        double ssbond=0;
        SSbondEnergyAtomAndAtom(pSG1,pSG2,pCB1,pCB2,pCA1,pCA2,&ssbond);
        energyTerms[36]+=ssbond;
      }
    }
  }


  return Success;
}

int EVOEF_EnergyProteinRotamerAndFixedResidueDifferentChain(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isBBAtom) continue;
    for(int j=0; j <ResidueGetAtomCount(pOther); j++){
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      //if(pAtom2->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0, rep=0, ele=0, desP=0, desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[51] += att;
      energyTerms[52] += rep;
      energyTerms[53] += ele;
      energyTerms[54] += desP;
      energyTerms[55] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          //only hbscbb or hbscsc
          if(pAtom2->isBBAtom == TRUE){
            energyTerms[64] += hbd;
            energyTerms[65] += hbt;
            energyTerms[66] += hbp;
          }
          else{
            energyTerms[67] += hbd;
            energyTerms[68] += hbt;
            energyTerms[69] += hbp;
          }
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          //only hbscbb or hbscsc
          if(pAtom2->isBBAtom == TRUE){
            energyTerms[64] += hbd;
            energyTerms[65] += hbt;
            energyTerms[66] += hbp;
          }
          else{
            energyTerms[67] += hbd;
            energyTerms[68] += hbt;
            energyTerms[69] += hbp;
          }
        }
      }
    }
  }

  //consider the SSBOND
  if(strcmp("CYS",RotamerGetType(pThis))==0 && strcmp("CYS",ResidueGetName(pOther))==0){
    Atom* pSG1=NULL,*pSG2=NULL,*pCB1=NULL,*pCB2=NULL,*pCA1=NULL,*pCA2=NULL;
    for(int i=0;i<RotamerGetAtomCount(pThis);i++){
      Atom* pAtom=RotamerGetAtom(pThis,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA1=pAtom;
    }
    for(int i=0;i<ResidueGetAtomCount(pOther);i++){
      Atom* pAtom=ResidueGetAtom(pOther,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA2=pAtom;
    }
    if(pSG1!=NULL && pCB1!=NULL &&pSG2!=NULL && pCB2!=NULL && pCA1!=NULL && pCA2!=NULL){
      double dist=XYZDistance(&pSG1->xyz,&pSG2->xyz);
      if(dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN){
        double ssbond=0;
        SSbondEnergyAtomAndAtom(pSG1,pSG2,pCB1,pCB2,pCA1,pCA2,&ssbond);
        energyTerms[56]+=ssbond;
      }
    }
  }


  return Success;
}


int EVOEF_EnergyRotamerAndFixedLigandResidue(Rotamer* pThis, Residue* pLigand, double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isBBAtom) continue;
    for(int j=0; j <ResidueGetAtomCount(pLigand); j++){
      Atom* pAtom2 = ResidueGetAtom(pLigand, j);
      //if(pAtom2->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0, rep=0, ele=0, desP=0, desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[71] += att;
      energyTerms[72] += rep;
      energyTerms[73] += ele;
      energyTerms[74] += desP;
      energyTerms[75] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pLigand,AtomGetHbDorB(pAtom2)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          //only hbscsc
          energyTerms[84] += hbd;
          energyTerms[85] += hbt;
          energyTerms[86] += hbp;
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pLigand,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          //only hbscsc
          energyTerms[84] += hbd;
          energyTerms[85] += hbt;
          energyTerms[86] += hbp;
        }
      }
    }
  }

  return Success;
}


int EVOEF_EnergyLigandRotamerAndFixedResidue(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isBBAtom) continue;
    for(int j=0; j <ResidueGetAtomCount(pOther); j++){
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      //if(pAtom2->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0, rep=0, ele=0, desP=0, desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[71] += att;
      energyTerms[72] += rep;
      energyTerms[73] += ele;
      energyTerms[74] += desP;
      energyTerms[75] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          //only hbscbb or hbscsc
          if(pAtom2->isBBAtom == TRUE){
            energyTerms[81] += hbd;
            energyTerms[82] += hbt;
            energyTerms[83] += hbp;
          }
          else{
            energyTerms[84] += hbd;
            energyTerms[85] += hbt;
            energyTerms[86] += hbp;
          }
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          //only hbscbb or hbscsc
          if(pAtom2->isBBAtom == TRUE){
            energyTerms[81] += hbd;
            energyTerms[82] += hbt;
            energyTerms[83] += hbp;
          }
          else{
            energyTerms[84] += hbd;
            energyTerms[85] += hbt;
            energyTerms[86] += hbp;
          }
        }
      }
    }
  }

  return Success;
}



int EVOEF_EnergyProteinRotamerAndDesignResidueSameChain(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  int neighborCheck = 0;
  if(RotamerGetPosInChain(pThis)+1==ResidueGetPosInChain(pOther)) neighborCheck=12;
  else if(RotamerGetPosInChain(pThis)-1==ResidueGetPosInChain(pOther)) neighborCheck=21;

  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isBBAtom) continue;
    for(int j=0; j <ResidueGetAtomCount(pOther); j++){
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      //for designed residues,only consider its mainchain atoms
      if(pAtom2->isBBAtom==FALSE) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      int bondType=15;
      if(neighborCheck==12) bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1),AtomGetName(pAtom2),ResidueGetName(pOther));
      else if(neighborCheck==21) bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom2),AtomGetName(pAtom1),RotamerGetType(pThis));
      double att=0, rep=0, ele=0, desP=0, desH=0;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[31] += att;
      energyTerms[32] += rep;
      energyTerms[33] += ele;
      energyTerms[34] += desP;
      energyTerms[35] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          if(pThis->posInChain-pOther->posInChain<=2 && pThis->posInChain-pOther->posInChain>=-2){
            hbd *= HBOND_LOCAL_REDUCE;
            hbt *= HBOND_LOCAL_REDUCE;
            hbp *= HBOND_LOCAL_REDUCE;
          }
          //only hbscbb
          energyTerms[44] += hbd;
          energyTerms[45] += hbt;
          energyTerms[46] += hbp;
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          if(pThis->posInChain-pOther->posInChain<=2 && pThis->posInChain-pOther->posInChain>=-2){
            hbd *= HBOND_LOCAL_REDUCE;
            hbt *= HBOND_LOCAL_REDUCE;
            hbp *= HBOND_LOCAL_REDUCE;
          }
          //only hbscbb
          energyTerms[44] += hbd;
          energyTerms[45] += hbt;
          energyTerms[46] += hbp;
        }
       
      }
    }
  }

  return Success;
}



int EVOEF_EnergyProteinRotamerAndDesignResidueDifferentChain(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isBBAtom) continue;
    for(int j=0; j <ResidueGetAtomCount(pOther); j++){
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if(pAtom2->isBBAtom==FALSE) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0, rep=0, ele=0, desP=0, desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[51] += att;
      energyTerms[52] += rep;
      energyTerms[53] += ele;
      energyTerms[54] += desP;
      energyTerms[55] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          //only hbscbb
          energyTerms[64] += hbd;
          energyTerms[65] += hbt;
          energyTerms[66] += hbp;
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          //only hbscbb
          energyTerms[64] += hbd;
          energyTerms[65] += hbt;
          energyTerms[66] += hbp;
        }
        
      }
    }
  }

  return Success;
}

int EVOEF_EnergyLigandRotamerAndDesignResidue(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isBBAtom) continue;
    for(int j=0; j <ResidueGetAtomCount(pOther); j++){
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if(pAtom2->isBBAtom==FALSE) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0, rep=0, ele=0, desP=0, desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[71] += att;
      energyTerms[72] += rep;
      energyTerms[73] += ele;
      energyTerms[74] += desP;
      energyTerms[75] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          //only hbscbb
          energyTerms[84] += hbd;
          energyTerms[85] += hbt;
          energyTerms[86] += hbp;
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,15,&hbtot,&hbd,&hbt,&hbp);
          //only hbscbb
          energyTerms[84] += hbd;
          energyTerms[85] += hbt;
          energyTerms[86] += hbp;
        }

      }
    }
  }

  return Success;
}


#define DEAL_WITH_STATISTICAL_ENERGY_TERMS
int RamaTableReadFromFile(RamaTable* pRama,char* ramafile){
  FILE* fin=fopen(ramafile,"r");
  if(fin==NULL){
    return IOError;
  }
  char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(fgets(buffer,MAX_LENGTH_ONE_LINE_IN_FILE,fin)){
    if(buffer[0]==' '|| buffer[0]=='#') continue;
    int phipsi[2];
    double aae[20];
    sscanf(buffer,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
      &phipsi[0],&phipsi[1],&aae[0],&aae[1],&aae[2],&aae[3],&aae[4],&aae[5],&aae[6],&aae[7],&aae[8],&aae[9],
      &aae[10],&aae[11],&aae[12],&aae[13],&aae[14],&aae[15],&aae[16],&aae[17],&aae[18],&aae[19]);
    int phiindex=(phipsi[0]+180)/10;
    int psiindex=(phipsi[1]+180)/10;
    for(int j=0;j<20;j++){
      pRama->ramatable[phiindex][psiindex][j]=aae[j];
    }
  }
  fclose(fin);
  return Success;
}

int AApropensityTableReadFromFile(AAppTable* pAAppTable,char* aappfile){
  FILE* fin=fopen(aappfile,"r");
  if(fin==NULL){
    return IOError;
  }
  char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(fgets(buffer,MAX_LENGTH_ONE_LINE_IN_FILE,fin)){
    if(buffer[0]==' '|| buffer[0]=='#') continue;
    int phipsi[2];
    double aae[20];
    sscanf(buffer,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
      &phipsi[0],&phipsi[1],&aae[0],&aae[1],&aae[2],&aae[3],&aae[4],&aae[5],&aae[6],&aae[7],&aae[8],&aae[9],
      &aae[10],&aae[11],&aae[12],&aae[13],&aae[14],&aae[15],&aae[16],&aae[17],&aae[18],&aae[19]);
    int phiindex=(phipsi[0]+180)/10;
    int psiindex=(phipsi[1]+180)/10;
    for(int j=0;j<20;j++){
      pAAppTable->aapptable[phiindex][psiindex][j]=aae[j];
    }
  }
  fclose(fin);
  return Success;
}

int AminoAcidPropensityAndRamachandranEnergy(Residue* pThis,AAppTable* pAAppTable,RamaTable* pRama,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  char aa1=ThreeLetterAAToOneLetterAA(ResidueGetName(pThis));
  int aaindex=AminoAcidGetIndex(aa1);
  int phiindex=((int)pThis->phipsi[0]+180)/10;
  int psiindex=((int)pThis->phipsi[1]+180)/10;
  energyTerms[91]+=pAAppTable->aapptable[phiindex][psiindex][aaindex];
  energyTerms[92]+=pRama->ramatable[phiindex][psiindex][aaindex];
  energyTerms[93]+=0;
  return Success;
}


int EVOEF_RotamerPropensityAndRamachandranEnergy(Rotamer* pThis,Residue* pResidue,AAppTable* pAAppTable,RamaTable* pRama,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  char aa1=ThreeLetterAAToOneLetterAA(RotamerGetType(pThis));
  int aaindex=AminoAcidGetIndex(aa1);
  int phiindex=((int)pResidue->phipsi[0]+180)/10;
  int psiindex=((int)pResidue->phipsi[1]+180)/10;
  energyTerms[91]+=pAAppTable->aapptable[phiindex][psiindex][aaindex];
  energyTerms[92]+=pRama->ramatable[phiindex][psiindex][aaindex];
  energyTerms[93]+=pThis->dunbrack;
  return Success;
}


int AminoAcidDunbrackEnergy(Residue* pThis,BBdepRotamerLib* pBBdepRotLib,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  if(strcmp(ResidueGetName(pThis),"ALA")==0 || strcmp(ResidueGetName(pThis),"GLY")==0){
    energyTerms[93]+=0;
    return Success;
  }
  double DUNBRACK_ENERGY=0.0;
  int binindex=((int)(pThis->phipsi[0]+180)/10)*36+(int)(pThis->phipsi[1]+180)/10;
  RotLibPhiPsi* pRotLibPhiPsi=&pBBdepRotLib->rotlibphipsis[binindex];
  int rotTypeIndex=-1;
  StringArrayFind(&pRotLibPhiPsi->rotTypes,ResidueGetName(pThis),&rotTypeIndex);
  DoubleArray* pTorsionsArrayForTypeI = pRotLibPhiPsi->torsions[rotTypeIndex];
  DoubleArray* pDeviationsArrayForTypeI = pRotLibPhiPsi->deviations[rotTypeIndex];
  int matchIndex=-1;
  for(int i=0;i<IntArrayGet(&pRotLibPhiPsi->rotamerCounts,rotTypeIndex);i++){
    DoubleArray* pTorsions=&pTorsionsArrayForTypeI[i];
    DoubleArray* pDeviations=&pDeviationsArrayForTypeI[i];
    BOOL match=TRUE;
    for(int j=0;j<DoubleArrayGetLength(&pThis->xtorsions);j++){
      double min=DoubleArrayGet(pTorsions,j)-DegToRad(30);
      double max=DoubleArrayGet(pTorsions,j)+DegToRad(30);
      double torsion=DoubleArrayGet(&pThis->xtorsions,j);
      double torsionm2pi=torsion-2*PI;
      double torsionp2pi=torsion+2*PI;
      double torsion2=torsion;
      if((strcmp(ResidueGetName(pThis),"PHE")==0 && j==1)||
        (strcmp(ResidueGetName(pThis),"TYR")==0 && j==1)||
        (strcmp(ResidueGetName(pThis),"ASP")==0 && j==1)||
        strcmp(ResidueGetName(pThis),"GLU")==0 && j==2){
          torsion2=torsion+PI;
          torsion2=torsion>0?torsion-PI:torsion2;
      }
      double torsion2m2pi=torsion2-2*PI;
      double torsion2p2pi=torsion2+2*PI;
      if(!(
        (torsion    <=max && torsion>=min) ||
        (torsionm2pi<=max && torsionm2pi>=min) ||
        (torsionp2pi<=max && torsionp2pi>=min) ||
        (torsion2    <=max && torsion2>=min) ||
        (torsion2m2pi<=max && torsion2m2pi>=min) ||
        (torsion2p2pi<=max && torsion2p2pi>=min)
        )){
          match=FALSE;
          break;
      }
    }
    if(match==TRUE){
      matchIndex=i;
      break;
    }
  }
  double delta_prob=1e-7;
  if(matchIndex!=-1){
    DUNBRACK_ENERGY=-1.0*log(DoubleArrayGet(&pRotLibPhiPsi->probability[rotTypeIndex],matchIndex)+delta_prob);
  }
  else{
    //printf("cannot find similar rotamer for residue %s %d %s in bbdep rotlib\n",ResidueGetChainName(pDestResidue),ResidueGetPosInChain(pDestResidue),ResidueGetName(pDestResidue));
    DUNBRACK_ENERGY=-1.0*log(DoubleArrayGet(&pRotLibPhiPsi->probability[rotTypeIndex],DoubleArrayGetLength(&pRotLibPhiPsi->probability[rotTypeIndex])-1)+delta_prob);
  }
  energyTerms[93]+=DUNBRACK_ENERGY;
  return Success;
}


int EVOEF_RotamerDunbrackEnergy(Rotamer* pThis,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  energyTerms[93]+=pThis->dunbrack;
  return Success;
}

