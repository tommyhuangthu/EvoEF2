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

#ifndef EEF1_H
#define EEF1_H

#include "Atom.h"

#define  EEF1_volume_CO        14.7
#define  EEF1_volume_COO       14.7
#define  EEF1_volume_CR         8.3
#define  EEF1_volume_CH1E      23.7
#define  EEF1_volume_CH2E      22.4
#define  EEF1_volume_CH3E      30.0
#define  EEF1_volume_CR1E      18.4
#define  EEF1_volume_Carg       8.3
#define  EEF1_volume_Hpol       0.0
#define  EEF1_volume_Hnpl       0.0
#define  EEF1_volume_OH1       10.8
#define  EEF1_volume_OH2       10.8
#define  EEF1_volume_Oest      10.8
#define  EEF1_volume_OC        10.8
#define  EEF1_volume_OOC       10.8
#define  EEF1_volume_NH1        4.4
#define  EEF1_volume_NR         4.4
#define  EEF1_volume_NH2       11.2
#define  EEF1_volume_NH3       11.2
#define  EEF1_volume_Narg      11.2
#define  EEF1_volume_Npro       0.0
#define  EEF1_volume_S         14.7
#define  EEF1_volume_SH1E      21.4
#define  EEF1_volume_P         14.7
#define  EEF1_volume_B          8.3
#define  EEF1_volume_F         10.8
#define  EEF1_volume_Cl        14.7
#define  EEF1_volume_Br        14.7
#define  EEF1_volume_I         14.7
#define  EEF1_volume_Fe         4.4
#define  EEF1_volume_Zn         4.4
#define  EEF1_volume_Al         4.4
#define  EEF1_volume_Mg         4.4
#define  EEF1_volume_Ca         4.4
#define  EEF1_volume_Na         4.4
#define  EEF1_volume_K          4.4
#define  EEF1_volume_Other      0.0

#define  EEF1_dg_free_CO        0.00
#define  EEF1_dg_free_COO       0.00
#define  EEF1_dg_free_CR        0.08
#define  EEF1_dg_free_CH1E      0.25
#define  EEF1_dg_free_CH2E      0.52
#define  EEF1_dg_free_CH3E      1.50
#define  EEF1_dg_free_CR1E      0.80
#define  EEF1_dg_free_Carg     -1.40
#define  EEF1_dg_free_Hpol      0.00
#define  EEF1_dg_free_Hnpl      0.00
#define  EEF1_dg_free_OH1      -6.70
#define  EEF1_dg_free_OH2      -6.70
#define  EEF1_dg_free_Oest     -6.70
#define  EEF1_dg_free_OC       -5.85
#define  EEF1_dg_free_OOC     -10.00
#define  EEF1_dg_free_NH1      -8.90
#define  EEF1_dg_free_NR       -4.00
#define  EEF1_dg_free_NH2      -7.80
#define  EEF1_dg_free_NH3     -20.00
#define  EEF1_dg_free_Narg    -10.00
#define  EEF1_dg_free_Npro     -1.55
#define  EEF1_dg_free_S         0.52//same as CH2E
#define  EEF1_dg_free_SH1E      1.50//same as CH3E
#define  EEF1_dg_free_P         0.52//same as S
#define  EEF1_dg_free_B         0.00
#define  EEF1_dg_free_F        -6.70//same as OH1
#define  EEF1_dg_free_Cl        0.52//same as S
#define  EEF1_dg_free_Br        0.52//same as S
#define  EEF1_dg_free_I         0.52//same as S
#define  EEF1_dg_free_Fe        0.00
#define  EEF1_dg_free_Zn        0.00
#define  EEF1_dg_free_Al        0.00
#define  EEF1_dg_free_Mg        0.00
#define  EEF1_dg_free_Ca        0.00
#define  EEF1_dg_free_Na        0.00
#define  EEF1_dg_free_K         0.00
#define  EEF1_dg_free_Other     0.00


#define  EEF1_lamda_charged     6.0
#define  EEF1_lamda_other       3.5

//from CHARMM19
#define  epsilon_CO             0.1200
#define  epsilon_COO            0.1200
#define  epsilon_CR             0.1200
#define  epsilon_CH1E           0.0486
#define  epsilon_CH2E           0.1142
#define  epsilon_CH3E           0.1811
#define  epsilon_CR1E           0.1200
#define  epsilon_Carg           0.1200
#define  epsilon_Hpol           0.0498
#define  epsilon_Hnpl           0.0450
#define  epsilon_OH1            0.1591
#define  epsilon_OH2            0.1591
#define  epsilon_Oest           0.1591
#define  epsilon_OC             0.6469
#define  epsilon_OOC            0.6469
#define  epsilon_NH1            0.2384
#define  epsilon_NH2            0.2384
#define  epsilon_NH3            0.2384
#define  epsilon_NR             0.2384
#define  epsilon_Narg           0.2384
#define  epsilon_Npro           0.2384
#define  epsilon_P              0.0430
#define  epsilon_S              0.0430
#define  epsilon_SH1E           0.0430
#define  epsilon_B              0.0486
#define  epsilon_F              0.1591
#define  epsilon_Cl             0.0430
#define  epsilon_Br             0.0430
#define  epsilon_I              0.0430
#define  epsilon_Fe             0.0000
#define  epsilon_Zn             0.0000
#define  epsilon_Al             0.0000
#define  epsilon_Mg             0.0000
#define  epsilon_Ca             0.0000
#define  epsilon_Na             0.0000
#define  epsilon_K              0.0000
#define  epsilon_Other          0.0000

#define  radius_CO              2.1000
#define  radius_COO             2.1000
#define  radius_CR              2.1000
#define  radius_CH1E            2.3650
#define  radius_CH2E            2.2350
#define  radius_CH3E            2.1650
#define  radius_CR1E            2.1000
#define  radius_Carg            2.1000
#define  radius_Hpol            0.8000
#define  radius_Hnpl            1.4680
#define  radius_OH1             1.6000
#define  radius_OH2             1.6000
#define  radius_Oest            1.6000
#define  radius_OC              1.6000
#define  radius_OOC             1.6000
#define  radius_NH1             1.6000
#define  radius_NH2             1.6000
#define  radius_NH3             1.6000
#define  radius_NR              1.6000
#define  radius_Narg            1.6000
#define  radius_Npro            1.6000
#define  radius_P               1.8900
#define  radius_S               1.8900
#define  radius_SH1E            1.8900
#define  radius_B               2.1000
#define  radius_F               1.6000
#define  radius_Cl              1.8900
#define  radius_Br              1.8900
#define  radius_I               1.8900
#define  radius_Fe              0.6500
#define  radius_Zn              0.6500
#define  radius_Al              0.6500
#define  radius_Mg              0.6500
#define  radius_Ca              0.6500
#define  radius_Na              0.6500
#define  radius_K               0.6500
#define  radius_Other           0.6500


typedef enum _EEF1_AtomType{
  EEF1_AtomType_CO, 
  EEF1_AtomType_COO, 
  EEF1_AtomType_CR, 
  EEF1_AtomType_CH1E, 
  EEF1_AtomType_CH2E, 
  EEF1_AtomType_CH3E, 
  EEF1_AtomType_CR1E,
  EEF1_AtomType_Carg,
  EEF1_AtomType_Hpol, 
  EEF1_AtomType_Hnpl, 
  EEF1_AtomType_OH1, 
  EEF1_AtomType_OH2,
  EEF1_AtomType_Oest, 
  EEF1_AtomType_OC, 
  EEF1_AtomType_OOC, 
  EEF1_AtomType_NH1, 
  EEF1_AtomType_NH2, 
  EEF1_AtomType_NH3, 
  EEF1_AtomType_NR, 
  EEF1_AtomType_Narg, 
  EEF1_AtomType_Npro, 
  EEF1_AtomType_P, 
  EEF1_AtomType_S, 
  EEF1_AtomType_SH1E,
  EEF1_AtomType_B,
  EEF1_AtomType_F,
  EEF1_AtomType_Cl,
  EEF1_AtomType_Br,
  EEF1_AtomType_I,
  EEF1_AtomType_Fe,
  EEF1_AtomType_Zn,
  EEF1_AtomType_Al,
  EEF1_AtomType_Mg,
  EEF1_AtomType_Ca,
  EEF1_AtomType_Na,
  EEF1_AtomType_K,
  EEF1_AtomType_Other,
} EEF1_AtomType;

int AtomAssignParameter(Atom* pAtom,char* atomType,double dgfree,double lambda,double volume,double epsilon,double radius);
int AssignAtomParameterByEEF1Type(Atom* pAtom,int eef1type);

#endif // EEF1_H