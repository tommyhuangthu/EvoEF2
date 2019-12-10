///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#ifndef ROTAMER_OPTIMIZER_H
#define ROTAMER_OPTIMIZER_H

#include "Structure.h"
#include "EnergyFunction.h"


int ProteinSiteOptimizeRotamer(Structure *pStructure, int chainIndex, int resiIndex);
int ProteinSiteOptimizeRotamerHBondEnergy(Structure *pStructure, int chainIndex, int resiIndex);
int ProteinSiteOptimizeRotamerLocally(Structure *pStructure, int chainIndex, int resiIndex, double rmsdcutoff);
int ProteinSiteCalcRotamersEnergy(Structure *pStructure,AAppTable* pAAppTable,RamaTable* pRama,int chainIndex,int resiIndex,FILE* fp);
int ProteinSiteOptimizeRotamerWithBBdepRotLib(Structure *pStructure, int chainIndex, int resiIndex,BBdepRotamerLib *pBBdepRotLib);
int ProteinSiteOptimizeRotamerLocallyWithBBdepRotLib(Structure *pStructure, int chainIndex, int resiIndex,double rmsdcutoff,BBdepRotamerLib *pBBdepRotLib);


#endif