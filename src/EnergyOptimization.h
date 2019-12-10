///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#ifndef ENERGY_OPTIMIZATION_H
#define ENERGY_OPTIMIZATION_H

#include "EnergyOptimization.h"
#include "DEE.h"
#include "EnergyMatrix.h"
#include "SmallMol.h"
#include "Sequence.h"

//Common parameters
#define STEP_SCALE   100
#define	MUT_NUM      1

//Monte Carlo
#define	MC_TMIN      0.05
#define	MC_STEP      5000000

//Simulated Annealing
#define SA_CYCLE     3
#define	SA_TMAX      5.0
#define	SA_TMIN	     0.03
#define	SA_TDEC      0.8
#define	SA_STEP      20000

//Reannealing
#define RA_TMAX      10.0
#define RA_TMIN      0.03
#define RA_CYCLE     100 //should not be too large, [100,200]?


int SitePairConsDeploy(CataConsSitePairArray *pSitePairArray, Structure *pStructure);
int EnergyMatrixUpdateForCataConsNew(EnergyMatrix *pMatrix, RotamerList* pList, CataConsSitePair *pSitePair, Structure *pStructure);
int EnergyMatrixUpdateForCataConsArrayNew(EnergyMatrix *pMatrix, RotamerList* pList, CataConsSitePairArray *pSitePairArray, Structure *pStructure);
int DesignSiteGetRotamerTypeAndCount(RotamerList* pList, Structure* pStructure, StringArray** ppRotamerType, IntArray** ppRotamerCount);

int SequenceRandomSingleSiteIndex(int *mutSiteIndex, int designSiteCount);
int SequenceRandomSingleRotamerIndex(Structure* pStructure, RotamerList* pList, int mutSiteIndex, int *mutRotIndex);
int SequenceRandomSingleRotamerIndexNew(Structure* pStructure,RotamerList* pList,StringArray** ppRotamerType,IntArray** ppRotamerCount,int siteIndex,int* rotIndex);
int SequenceChangeSingleSite(Sequence* pThis, int mutationSiteIndex, int mutationRotamerIndex);

int StructureCalcSequenceTemplateEnergy(Structure* pStructure, Sequence* pSequence, double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]);
int StructureCalcSequenceEnergy(Structure* pStructure,Sequence* pSequence,char* desChnName,char* programpath,char* prf,char* ss,char* sa,char* phipsi);
int StructureCalcEnergyChangeUponSingleMutation(Structure* pStructure,Sequence* pSequence,int mutSiteIndex,int mutRotIndex,double* dtot,double* dphy,double *devo,char* desChnName,char* programpath,char* prf,char* ss, char* sa,char* phipsi);


int SequenceGenerateRandomSeed(Sequence* pThis, RotamerList* pList);
int SequenceGenerateInitialSequenceSeed(Structure* pStructure, RotamerList* pList,Sequence* pThis);
int SequenceGenerateCrystalSequenceSeed(Structure* pStructure, RotamerList* pList,Sequence* pSequence);
int MetropolisCriteria(Sequence* pOld,Sequence* pBest,Structure* pStructure,RotamerList* pList,int *seqIndex,double temp,int stepCount, FILE* fp,FILE *fp2,char* deschn,char* programpath, char* prf, char* ss, char* sa, char* phipsi);
int MetropolisCriteriaNew(Sequence* pOld,Sequence* pBest,Structure* pStructure,RotamerList* pList,StringArray** ppRotamerType,IntArray** ppRotamerCount,int *seqIndex,double temp,int stepCount, FILE* fp,FILE *fp2,char* deschn,char* programpath, char* prf, char* ss, char* sa, char* phipsi);

int MonteCarloOptimization(Structure* pStructure,RotamerList* pList,char* deschn,char* programpath,char* prf, char* ss, char* sa, char* phipsi,char* rotindexfile,char* fastafile);
int SimulatedAnnealingOptimization(Structure* pStructure,RotamerList* pList,char* deschn,char* programpath,char* prf, char* ss, char* sa, char* phipsi,char* design_rotindex_file,char* fastafile);
int ReannealingOptimization(Structure* pStructure,RotamerList* pList,char* deschn,char* programpath,char* prf, char* ss, char* sa, char* phipsi,char* design_rotindex_file,char* fastafile);

int SimulatedAnnealingOptimizationForSCP(Structure* pStructure,RotamerList* pList,char* deschn,char* programpath,char* prf,char* ss,char* sa,char* phipsi,char* design_rotindex_file,char* fastafile);
int MetropolisCriteriaForSCP(Sequence* pOld,Sequence* pBest,Structure* pStructure,RotamerList* pList,StringArray** ppRotamerType,IntArray** ppRotamerCount,int *seqIndex,double temp,int stepCount, FILE* fp,FILE *fp2,char* deschn,char* programpath, char* prf, char* ss, char* sa, char* phipsi);


#endif
