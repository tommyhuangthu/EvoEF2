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