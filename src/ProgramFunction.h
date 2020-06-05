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

#ifndef PROGRAM_FUNCTION_H
#define PROGRAM_FUNCTION_H

#include "Structure.h"
#include "EnergyFunction.h"


int EVOEF_help();
int EVOEF_version();
int EVOEF_interface();
int EVOEF_ComputeChainStability(Structure *pStructure, int chainIndex, double *energyTerms);
int EVOEF_ComputeStability(Structure *pStructure,AAppTable* pAAppTable,RamaTable* pRama,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_ComputeStabilityWithBBdepRotLib(Structure *pStructure,AAppTable* pAAppTable,RamaTable* pRama,BBdepRotamerLib* pRotLib,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]);
int EVOEF_ComputeStabilityForSelectedChains(Structure *pStructure, double *energyTerms,char selechains[]);
int EVOEF_ComputeBinding(Structure *pStructure);
int EVOEF_ComputeBindingWithSplitting(Structure *pStructure,char split1[], char split2[]);
int EVOEF_ComputeBindingWithSplittingNew(Structure *pStructure,char split1[], char split2[]);
int EVOEF_BuildMutant(Structure* pStructure, char* mutantfile, BBindRotamerLib* rotlib, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid);
int EVOEF_BuildMutantWithBBdepRotLib(Structure* pStructure, char* mutantfile, BBdepRotamerLib* pBBdepRotLib, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid);
int EVOEF_RepairStructure(Structure* pStructure, BBindRotamerLib* rotlib, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid);
int EVOEF_WriteStructureToFile(Structure* pStructure, char* pdbfile);
int EVOEF_AddHydrogen(Structure* pStructure, char* pdbid);
int EVOEF_OptimizeHydrogen(Structure* pStructure, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid);
int EVOEF_StructureComputeResidueInteractionWithFixedSurroundingResidues(Structure *pStructure, int chainIndex, int residueIndex);


int EVOEF_ComputeRotamersEnergy(Structure* pStructure,BBindRotamerLib* rotlib,AAppTable* pAAppTable,RamaTable* pRama,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid);
int EVOEF_ComputeWildtypeRotamersEnergy(Structure* pStructure,BBindRotamerLib* rotlib,AAppTable* pAAppTable,RamaTable* pRama,AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid);

int EVOEF_StructureFindInterfaceResidues(Structure *pStructure,double cutoff,char* outputfile);
int StructureCalcPhiPsi(Structure* pStructure);
int EVOEF_StructureShowPhiPsi(Structure* pStructure,char* phipsifile);
int EVOEF_ComputeRotamersEnergyByBBdepRotLib(Structure* pStructure,BBdepRotamerLib* rotlib,AAppTable* pAAppTable,RamaTable* pRama,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid);
int EVOEF_ComputeWildtypeRotamersEnergyByBBdepRotLib(Structure* pStructure,BBdepRotamerLib* rotlib,AAppTable* pAAppTable,RamaTable* pRama,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid);

int EVOEF_CheckRotamerInBBindRotLib(Structure* pStructure,BBindRotamerLib* pRotLib,ResiTopoSet* pTopos,double cutoff,char* pdbid);
int EVOEF_CheckRotamerInBBdepRotLib(Structure* pStructure,BBdepRotamerLib* pRotLib,ResiTopoSet* pTopos,double cutoff,char* pdbid);
int EVOEF_RepairStructureWithBBdepRotLib(Structure* pStructure, BBdepRotamerLib* pBBdepRotLib, AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid);
int EVOEF_GetResiMinRmsdRotFromLab(Structure* pStructure,char* pdbid);



#endif //PROGRAM_FUNCTION_H