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

#ifndef ROTAMER_BUILDER_H
#define ROTAMER_BUILDER_H

#include "Structure.h"

#define NATIVE_ROTAMER_DUNBRACK 0.0//kcal

#define EXPANDED_ROT_SER  5
#define EXPANDED_ROT_THR  5
#define EXPANDED_ROT_TYR  1


int ProteinSiteBuildAllRotamers(Structure* pThis, int chainIndex, int resiIndex, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int ProteinSiteBuildMutatedRotamers(Structure* pThis, int chainIndex, int resiIndex, BBindRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos, StringArray *pDesignTypes, StringArray *pPatchTypes);
int ProteinSiteBuildWildtypeRotamers(Structure* pThis, int chainIndex, int resiIndex, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int ProteinSiteAddCrystalRotamer(Structure* pStructure, int chainIndex, int resiIndex, ResiTopoSet *pResiTopos);
int ProteinSiteExpandHydroxylRotamers(Structure *pStructure, int chainIndex, int resiIndex, ResiTopoSet *pTopos);
int ProteinSiteBuildFlippedCrystalRotamer(Structure* pStructure, int chainIndex, int resiIndex, ResiTopoSet *pResiTopos);
int ProteinSiteWriteRotamers(Structure *pStructure, int chainIndex, int resiIndex, char *rotamerFilePath);




int StructureGenerateSpecifiedProteinRotamers(Structure* pThis, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* designsitefile);
int StructureGenerateBindingSiteRotamers(Structure* pThis, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, double range);
int StructureGenerateAllRotamers(Structure* pThis,BBindRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* deschns);
int StructureGenerateWildtypeRotamers(Structure* pThis,BBindRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* DESCHNS);


//deal with small molecule
//int StructureGetTruncatedBackbone(Structure* pThis,Residue* pSmallMol, double activeSiteRange, BOOL withHydrogen,AtomArray* pBackboneAtoms);
//int StructureDeployCataConsSitePair(Structure* pThis, CataConsSitePair* pCataConsSitePair);
//int StructurePlaceSmallMol(Structure* pThis, PlacingRule* pPlacingRule, CataConsSitePairArray* pCataConsCollection, int relatedProteinSiteCount,DesignSite** relatedProteinSites, RotamerSet* pSmallMolRotSet);
//int StructureGenerateSmallMolRotamers(Structure* pThis,char* cataConsFileName, char* placingRuleFileName);
//int StructureReadSmallMolRotamers(Structure* pThis,ResiTopoSet* resiTopos,char* smallMolFileName);
//int StructureWriteSmallMolRotamers(Structure* pThis,char* smallMolFile);
//int StructureSmallmolOrientationScreen(Structure* pStructure, ResiTopoSet* pResiTopo, char* oriFileName, char* newFileName, char* screenRuleFileName);


int ProteinSiteBuildAllRotamersByBBdepRotLib(Structure* pThis,int chainIndex,int resiIndex,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos);
int ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(Structure* pThis, int chainIndex, int resiIndex,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos, Type_ResidueDesignType type);
int ProteinSiteBuildMutatedRotamersByBBdepRotLib(Structure* pThis, int chainIndex, int resiIndex, BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos, StringArray *pDesignTypes, StringArray *pPatchTypes);
int ProteinSiteBuildWildtypeRotamersByBBdepRotLib(Structure* pThis,int chainIndex,int resiIndex,BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int StructureGenerateSpecifiedProteinRotamersByBBdepRotLib(Structure* pThis,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* designsitefile);
int StructureGenerateBindingSiteRotamersByBBdepRotLib(Structure* pThis,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,double range);
int StructureGenerateAllRotamersByBBdepRotLib(Structure* pThis,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* DESCHNS);
int ProteinSiteAddCrystalRotamerByBBdepRotLib(Structure* pThis,int chainIndex,int resiIndex,ResiTopoSet *pResiTopos,BBdepRotamerLib* pBBdepRotLib);
int StructureGenerateWildtypeRotamersByBBdepRotLib(Structure* pThis,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* DESCHNS);


BOOL ProteinSiteCheckCrystalRotamerInBBdepRotLib(Structure* pThis,int chainIndex,int resiIndex,ResiTopoSet *pResiTopos,BBdepRotamerLib* pBBdepRotLib,double torsionStd);
BOOL ProteinSiteCheckCrystalRotamerInBBindRotLib(Structure* pThis,int chainIndex,int resiIndex,ResiTopoSet *pResiTopos,BBindRotamerLib* pBBindRotLib,double torsionStd);

int StructureGenerateCatalyticSiteProteinRotamersByBBdepRotLib(Structure* pThis,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* catasitefile);
int StructureGenerateFirstShellProteinRotamersByBBdepRotLib(Structure* pThis,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* shell1sitefile,double shell1cutoff);
int StructureGenerateSecondShellProteinRotamersByBBdepRotLib(Structure* pThis,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* shell2sitefile,double shell2cutoff);


#endif