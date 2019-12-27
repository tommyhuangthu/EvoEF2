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

#include "EEF1Ligand.h"
#include <string.h>
#include <stdio.h>



int AtomAssignParameter(Atom* pAtom,char* atomType,double dgfree,double lambda,double volume,double epsilon,double radius){
  strcpy(pAtom->type,atomType);
  pAtom->EEF1_freeDG=dgfree;
  pAtom->EEF1_lamda_=lambda;
  pAtom->EEF1_volume=volume;
  pAtom->vdw_epsilon=epsilon;
  pAtom->vdw_radius=radius;
  //pAtom->charge=charge;
  return Success;
}

int AssignAtomParameterByEEF1Type(Atom* pAtom,int eef1type){
  switch(eef1type){
    case EEF1_AtomType_CO:
      AtomAssignParameter(pAtom,"CO",EEF1_dg_free_CO,EEF1_lamda_other,EEF1_volume_CO,epsilon_CO,radius_CO);
      break;
    case EEF1_AtomType_COO:
      AtomAssignParameter(pAtom,"COO",EEF1_dg_free_COO,EEF1_lamda_other,EEF1_volume_COO,epsilon_COO,radius_COO);
      break;
    case EEF1_AtomType_CR:
      AtomAssignParameter(pAtom,"CR",EEF1_dg_free_CR,EEF1_lamda_other,EEF1_volume_CR,epsilon_CR,radius_CR);
      break;
    case EEF1_AtomType_CH1E:
      AtomAssignParameter(pAtom,"CH1E",EEF1_dg_free_CH1E,EEF1_lamda_other,EEF1_volume_CH1E,epsilon_CH1E,radius_CH1E);
      break;
    case EEF1_AtomType_CH2E:
      AtomAssignParameter(pAtom,"CH2E",EEF1_dg_free_CH2E,EEF1_lamda_other,EEF1_volume_CH2E,epsilon_CH2E,radius_CH2E);
      break;
    case EEF1_AtomType_CH3E:
      AtomAssignParameter(pAtom,"CH3E",EEF1_dg_free_CH3E,EEF1_lamda_other,EEF1_volume_CH3E,epsilon_CH3E,radius_CH3E);
      break;
    case EEF1_AtomType_CR1E:
      AtomAssignParameter(pAtom,"CR1E",EEF1_dg_free_CR1E,EEF1_lamda_other,EEF1_volume_CR1E,epsilon_CR1E,radius_CR1E);
      break;
    case EEF1_AtomType_Carg:
      AtomAssignParameter(pAtom,"Carg",EEF1_dg_free_Carg,EEF1_lamda_other,EEF1_volume_Carg,epsilon_Carg,radius_Carg);
      break;
    case EEF1_AtomType_Hpol:
      AtomAssignParameter(pAtom,"Hpol",EEF1_dg_free_Hpol,EEF1_lamda_other,EEF1_volume_Hpol,epsilon_Hpol,radius_Hpol);
      break;
    case EEF1_AtomType_Hnpl:
      AtomAssignParameter(pAtom,"Hnpl",EEF1_dg_free_Hnpl,EEF1_lamda_other,EEF1_volume_Hnpl,epsilon_Hnpl,radius_Hnpl);
      break;
    case EEF1_AtomType_OH1:
      AtomAssignParameter(pAtom,"OH1",EEF1_dg_free_OH1,EEF1_lamda_other,EEF1_volume_OH1,epsilon_OH1,radius_OH1);
      break;
    case EEF1_AtomType_OH2:
      AtomAssignParameter(pAtom,"OH2",EEF1_dg_free_OH2,EEF1_lamda_other,EEF1_volume_OH2,epsilon_OH2,radius_OH2);
      break;
    case EEF1_AtomType_Oest:
      AtomAssignParameter(pAtom,"Oest",EEF1_dg_free_Oest,EEF1_lamda_other,EEF1_volume_Oest,epsilon_Oest,radius_Oest);
      break;
    case EEF1_AtomType_OC:
      AtomAssignParameter(pAtom,"OC",EEF1_dg_free_OC,EEF1_lamda_other,EEF1_volume_OC,epsilon_OC,radius_OC);
      break;
    case EEF1_AtomType_OOC:
      AtomAssignParameter(pAtom,"OOC",EEF1_dg_free_OOC,EEF1_lamda_other,EEF1_volume_OOC,epsilon_OOC,radius_OOC);
      break;
    case EEF1_AtomType_NH1:
      AtomAssignParameter(pAtom,"NH1",EEF1_dg_free_NH1,EEF1_lamda_other,EEF1_volume_NH1,epsilon_NH1,radius_NH1);
      break;
    case EEF1_AtomType_NH2:
      AtomAssignParameter(pAtom,"NH2",EEF1_dg_free_NH2,EEF1_lamda_other,EEF1_volume_NH2,epsilon_NH2,radius_NH2);
      break;
    case EEF1_AtomType_NH3:
      AtomAssignParameter(pAtom,"NH3",EEF1_dg_free_NH3,EEF1_lamda_other,EEF1_volume_NH3,epsilon_NH3,radius_NH3);
      break;
    case EEF1_AtomType_NR:
      AtomAssignParameter(pAtom,"NR",EEF1_dg_free_NR,EEF1_lamda_other,EEF1_volume_NR,epsilon_NR,radius_NR);
      break;
    case EEF1_AtomType_Narg:
      AtomAssignParameter(pAtom,"Narg",EEF1_dg_free_Narg,EEF1_lamda_other,EEF1_volume_Narg,epsilon_Narg,radius_Narg);
      break;
    case EEF1_AtomType_Npro:
      AtomAssignParameter(pAtom,"Npro",EEF1_dg_free_Npro,EEF1_lamda_other,EEF1_volume_Npro,epsilon_Npro,radius_Npro);
      break;
    case EEF1_AtomType_P:
      AtomAssignParameter(pAtom,"P",EEF1_dg_free_P,EEF1_lamda_other,EEF1_volume_P,epsilon_P,radius_P);
      break;
    case EEF1_AtomType_S:
      AtomAssignParameter(pAtom,"S",EEF1_dg_free_S,EEF1_lamda_other,EEF1_volume_S,epsilon_S,radius_S);
      break;
    case EEF1_AtomType_SH1E:
      AtomAssignParameter(pAtom,"SH1E",EEF1_dg_free_SH1E,EEF1_lamda_other,EEF1_volume_SH1E,epsilon_SH1E,radius_SH1E);
      break;
    case EEF1_AtomType_B:
      AtomAssignParameter(pAtom,"B",EEF1_dg_free_B,EEF1_lamda_other,EEF1_volume_B,epsilon_B,radius_B);
      break;
    case EEF1_AtomType_F:
      AtomAssignParameter(pAtom,"F",EEF1_dg_free_F,EEF1_lamda_other,EEF1_volume_F,epsilon_F,radius_F);
      break;
    case EEF1_AtomType_Cl:
      AtomAssignParameter(pAtom,"Cl",EEF1_dg_free_Cl,EEF1_lamda_other,EEF1_volume_Cl,epsilon_Cl,radius_Cl);
      break;
    case EEF1_AtomType_Br:
      AtomAssignParameter(pAtom,"Br",EEF1_dg_free_Br,EEF1_lamda_other,EEF1_volume_Br,epsilon_Br,radius_Br);
      break;
    case EEF1_AtomType_I:
      AtomAssignParameter(pAtom,"I",EEF1_dg_free_I,EEF1_lamda_other,EEF1_volume_I,epsilon_I,radius_I);
      break;
    case EEF1_AtomType_Fe:
      AtomAssignParameter(pAtom,"Fe",EEF1_dg_free_Fe,EEF1_lamda_other,EEF1_volume_Fe,epsilon_Fe,radius_Fe);
      break;
    case EEF1_AtomType_Zn:
      AtomAssignParameter(pAtom,"Zn",EEF1_dg_free_Zn,EEF1_lamda_other,EEF1_volume_Zn,epsilon_Zn,radius_Zn);
      break;
    case EEF1_AtomType_Al:
      AtomAssignParameter(pAtom,"Al",EEF1_dg_free_Al,EEF1_lamda_other,EEF1_volume_Al,epsilon_Al,radius_Al);
      break;
    case EEF1_AtomType_Mg:
      AtomAssignParameter(pAtom,"Mg",EEF1_dg_free_Mg,EEF1_lamda_other,EEF1_volume_Mg,epsilon_Mg,radius_Mg);
      break;
    case EEF1_AtomType_Ca:
      AtomAssignParameter(pAtom,"Ca",EEF1_dg_free_Ca,EEF1_lamda_other,EEF1_volume_Ca,epsilon_Ca,radius_Ca);
      break;
    case EEF1_AtomType_Na:
      AtomAssignParameter(pAtom,"Na",EEF1_dg_free_Na,EEF1_lamda_other,EEF1_volume_Na,epsilon_Na,radius_Na);
      break;
    case EEF1_AtomType_K:
      AtomAssignParameter(pAtom,"K",EEF1_dg_free_K,EEF1_lamda_other,EEF1_volume_K,epsilon_K,radius_K);
      break;
    case EEF1_AtomType_Other:
      AtomAssignParameter(pAtom,"X",EEF1_dg_free_Other,EEF1_lamda_other,EEF1_volume_Other,epsilon_Other,radius_Other);
      break;
  }

  return Success;
}