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

#ifndef ENERGY_FUNCTION_H
#define ENERGY_FUNCTION_H
#include "Residue.h"
#include "Rotamer.h"


#define ENERGY_DEBUG_MODE_VDW_ATT           0
#define ENERGY_DEBUG_MODE_VDW_REP           0
#define ENERGY_DEBUG_MODE_HBOND             0
#define ENERGY_DEBUG_MODE_ELEC              0
#define ENERGY_DEBUG_MODE_DESOLV            0


// a maximum of energy weights
#define MAX_EVOEF_ENERGY_TERM_NUM          100

#define ENERGY_DISTANCE_CUTOFF             6.0
// scale for 1-2, 1-3, 1-4 and >=1-5 interactions
#define ENERGY_SCALE_FACTOR_BOND_123       0.0
#define ENERGY_SCALE_FACTOR_BOND_14        0.2
#define ENERGY_SCALE_FACTOR_BOND_15        1.0
// VDW interactions
#define RADIUS_SCALE_FOR_VDW              0.95

//HBond interactions
#define HBOND_DISTANCE_CUTOFF_MAX          3.0
#define HBOND_WELL_DEPTH                   1.0
#define HBOND_OPTIMAL_DISTANCE             1.9
//reduce local hbonds if |i-j|<=2
#define HBOND_LOCAL_REDUCE                 0.5


// Electrostatics interactions
#define ELEC_DISTANCE_CUTOFF               6.0
#define COULOMB_CONSTANT                 332.0
#define DIELECTRIC_CONST_PROTEIN           8.0
#define DIELECTRIC_CONSTANT_WATER         80.0
#define DIELECTRIC_CONST_PROTEIN_AVE      20.0
#define IONIC_STRENGTH                     0.05 // (unit: M)
#define PROTEIN_DESIGN_TEMPERATURE         298  // (unit: K)
// LK solvation interaction
#define LK_SOLV_DISTANCE_CUTOFF            6.0
#define RADIUS_SCALE_FOR_DESOLV            1.00


//SSBOND parameter
#define SSBOND_DISTANCE                    2.03
#define SSBOND_ANGLE                     105.0
#define SSBOND_TORSION                    90.0
#define SSBOND_CUTOFF_MAX                  2.15
#define SSBOND_CUTOFF_MIN                  1.95


int EnergyTermInitialize(double *energyTerms);
int EnergyTermWeighting(double *energyTerms);
int EnergyWeightRead(char* weightfile);
int EnergyWeightWrite(char* weightfile);

double CalcResidueBuriedRatio(Residue* pResi1);
double CalcAverageBuriedRatio(double ratio1, double ratio2);

BOOL ResidueIntraBond12Check(char *atom1,char *atom2, BondSet* pBondSet);
BOOL ResidueIntraBond13Check(char *atom1,char *atom2, BondSet* pBondSet);
BOOL ResidueIntraBond14Check(char *atom1,char *atom2, BondSet* pBondSet);
int ResidueIntraBondConnectionCheck(char *atom1,char *atom2, BondSet* pBondSet);
int ResidueAndNextResidueInterBondConnectionCheck_charmm22(char *atomOnPreResi, char *atomOnNextResi, Residue *pPreResi, Residue *pNextResi);
int ResidueAndNextResidueInterBondConnectionCheck_charmm19(char *atomOnPreResi, char *atomOnNextResi, char *nextResiName);

//physics- and knowledge-based energy terms
int VdwAttEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance, int bondType,double *vdwAtt);
int VdwRepEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance, int bondType, double *vdwRep);
int HBondEnergyAtomAndAtom(Atom *atomH, Atom *atomA, Atom*atomD, Atom *atomB,double distanceHA, int bondType, double *etotal, double *edist, double *etheta, double *ephi);
int ElecEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance12,int bondType, double *elec);
int LKDesolvationEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance,int bondType, double *energyP, double *energyH);
int SSbondEnergyAtomAndAtom(Atom *pAtomS1,Atom *pAtomS2,Atom *pAtomCB1,Atom *pAtomCB2,Atom* pAtomCA1,Atom* pAtomCA2,double* sse);


//reference energy to balance the composition of designed sequence
//make the designed sequences more native-like
int EVOEF_AminoAcidReferenceEnergy(char *AAname, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);


//residue-residue energy, including backbone and sidechain
int EVOEF_EnergyResidueIntraEnergy(Residue* pThis, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_EnergyResidueAndNextResidue(Residue* pThis, Residue* pOther, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_EnergyResidueAndOtherResidueSameChain(Residue* pThis, Residue* pOther, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_EnergyResidueAndOtherResidueDifferentChain(Residue* pThis, Residue* pOther, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_EnergyResidueAndLigandResidue(Residue* pProtein, Residue* pLigand,double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);


//protein rotamer-specific energy, similar to residue-specific energy
int EVOEF_EnergyProteinRotamerIntraEnergy(Rotamer* pThis, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_EnergyProteinRotamerAndRotamerSameChain(Rotamer* pThis, Rotamer* pOther, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_EnergyProteinRotamerAndRotamerDifferentChain(Rotamer* pThis, Rotamer* pOther, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_EnergyProteinRotamerAndFixedResidueSameChain(Rotamer* pThis, Residue* pOther, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_EnergyProteinRotamerAndFixedResidueDifferentChain(Rotamer* pThis, Residue* pOther, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_EnergyProteinRotamerAndDesignResidueSameChain(Rotamer* pThis, Residue* pOther, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_EnergyProteinRotamerAndDesignResidueDifferentChain(Rotamer* pThis, Residue* pOther, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);

//ligand related energy
int EVOEF_EnergyRotamerAndFixedLigandResidue(Rotamer* pThis, Residue* pLigand, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_EnergyLigandRotamerAndFixedResidue(Rotamer* pThis, Residue* pOther, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_EnergyLigandRotamerAndDesignResidue(Rotamer* pThis, Residue* pOther, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_EnergyLigandRotamerAndRotamer(Rotamer* pThis, Rotamer* pOther, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]);


typedef struct _RamaTable{
  double ramatable[36][36][20];
}RamaTable;

int RamaTableReadFromFile(RamaTable* pRama,char* ramafile);

typedef struct _AAppTable{
  double aapptable[36][36][20];
}AAppTable;

int AApropensityTableReadFromFile(AAppTable* pAAppTable,char* aappfile);
int EVOEF_RotamerPropensityAndRamachandranEnergy(Rotamer* pThis,Residue* pResidue,AAppTable* pAAppTable,RamaTable* pRama,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_RotamerDunbrackEnergy(Rotamer* pThis,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]);
int AminoAcidPropensityAndRamachandranEnergy(Residue* pThis,AAppTable* pAAppTable,RamaTable* pRama,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]);
int AminoAcidDunbrackEnergy(Residue* pThis,BBdepRotamerLib* pBBdepRotLib,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]);

#endif // ENERGY_FUNCTION_H


