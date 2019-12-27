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

#ifndef SMALL_MOL_H
#define SMALL_MOL_H

#include "DesignSite.h"
#include "Rotamer.h"

// a disturtance for the range check in the catalytic constraint checking functions;
// for example, if the upper bound of a hydrogen bond is 3.0 angstroms, the constraint
// will be sitisfied as long as the distance is less than 3.0+DISTURBANCE_IN_RANGE_CHECK;
#define DISTURBANCE_IN_RANGE_CHECK 1e-3 

#define MAX_COUNT_CHECK_MULTI_CONS 50

// an atom name prefixed with symbol '#' is a small molecule atom in the placing rule file;
#define SMALLMOL_ATOM_SYMBOL       '#'
// an atom name prefixed with symbol '#' is an atom on the second catalytic site in the catalytic constraint file;
#define SECOND_CATACON_SITE_SYMBOL '#'
// an atom name prefixed with symbol '+' is a pseudo atom in the placing rule file and the catalytic constraint file;
#define PSEUDO_ATOM_SYMBOL         '+'

//////////////////////////////////////////////////////////////
//The Catalytic Constraints Hierarchy
//////////////////////////////////////////////////////////////
typedef enum _Type_CataCons{
    Type_CataCons_NewPseudoAtom,
    Type_CataCons_Distance,
    Type_CataCons_Angle,
    Type_CataCons_Torsion
} Type_CataCons;
typedef enum _Type_CataConsAtomLocation{
    Type_CataConsAtomLocation_FirstSite,
    Type_CataConsAtomLocation_SecondSite,
    Type_CataConsAtomLocation_Pseudo,
    Type_CataConsAtomLocation_Undefined
} Type_CataConsAtomLocation;


typedef struct _CataConsItem{
	Type_CataCons type;
	char atomName[4][MAX_LENGTH_ATOM_NAME+1];
	Type_CataConsAtomLocation atomLocation[4];
	int  atomIndex[4];
	double min;
	double max;
} CataConsItem;
int    CataConsItemShow(CataConsItem* pThis);

typedef struct _CataConsGroup{
	char firstSiteRotamerType[MAX_LENGTH_RESIDUE_NAME+1];
	char secondSiteRotamerType[MAX_LENGTH_RESIDUE_NAME+1];

	int consItemsCount;
	CataConsItem* consItems;

	int pseudoAtomsCount;
	XYZ* pseudoAtoms;
} CataConsGroup;
int    CataConsGroupDeploy(CataConsGroup* pThis,Rotamer* pFirstSiteRotamer,Rotamer* pSecondSiteRotamer);

typedef struct _CataConsSitePair{

	//the following information is determined when the catalytic constraint file is read.
  char firstSiteChainName[MAX_LENGTH_CHAIN_NAME+1];
  int  firstSitePosInChain;
  char firstSiteResidueName[MAX_LENGTH_RESIDUE_NAME+1]; 
	char secondSiteChainName[MAX_LENGTH_CHAIN_NAME+1];
	int  secondSitePosInChain;
	char secondSiteResidueName[MAX_LENGTH_RESIDUE_NAME+1];

	//Memories will be allocated for 'groups' when the file have been read;
	//but the information is complete only when the CataConsSitePair is deployed;
	int groupsCount;
    CataConsGroup* groups;

	//a BOOL flag is used to record if the site pair has been deployed;
	//it is initialized as FALSE, and set as TRUE after the site pair is deployed;
	//the site pair is not allowed to use before the flag is deployed;
	//besides, the flag must be deployed only once.
    BOOL deployedFlag;

} CataConsSitePair;
int  CataConsSitePairCreate(CataConsSitePair* pThis,FileReader* file);
void CataConsSitePairDestroy(CataConsSitePair* pThis);
BOOL CataConsSitePairGetDeployedFlag(CataConsSitePair* pThis);
int  CataConsSitePairDeploy(CataConsSitePair* pThis,RotamerSet* pFirstSiteRotSet,RotamerSet* pSecondSiteRotSet);
int  CataConsSitePairShow(CataConsSitePair* pThis);


typedef struct _CataConsSitePairArray{
    int count;
    CataConsSitePair* sites;
} CataConsSitePairArray;
int           CataConsSitePairArrayCreate(CataConsSitePairArray* pThis,char* cataConsFile);
void          CataConsSitePairArrayDestroy(CataConsSitePairArray* pThis);
int           CataConsSitePairArrayGetCount(CataConsSitePairArray* pThis);
CataConsSitePair* CataConsSitePairArrayGet(CataConsSitePairArray* pThis,int index);
CataConsSitePair* CataConsSitePairArrayFind(CataConsSitePairArray* pThis,
									char* firstSiteChainName,
									int firstSitePosInChain,
									char* firstSiteResidueName,
									char* secondSiteChainName,
									int secondSitePosInChain,
									char* secondSiteResidueName);

int           CataConsSitePairArrayShow(CataConsSitePairArray* pThis);
int           CataConsSitePairArrayTester(char* filename);

//Catalytic Constraint Checking Operations
BOOL CataConsItemCheck(CataConsItem* pThis,XYZArray* pOnFirstSite, XYZArray* pOnSecondSite,XYZ* pseudoAtoms);
BOOL CataConsGroupCheck(CataConsGroup* pThis,Rotamer* pOnFirstSite,Rotamer* pOnSecondSite);
BOOL CataConsSitePairCheck(CataConsSitePair* pThis,Rotamer* pOnFirstSite,Rotamer* pOnSecondSite);

//Various Placing Actions
typedef enum _Type_PlacingAction{
    Type_PlacingAction_Load,
    Type_PlacingAction_Evaluate,
    Type_PlacingAction_Variate,
    Type_PlacingAction_Calc,
    Type_PlacingAction_CheckCataCons,
	Type_PlacingAction_CheckMultiCons,
    Type_PlacingAction_CheckVDW_Backbone,
	Type_PlacingAction_CheckVDW_Internal,
    Type_PlacingAction_CheckRMSD
} Type_PlacingAction;
typedef enum _Type_PlacingActionEvalParam{
    Type_PlacingActionEvaluate_Distance,
    Type_PlacingActionEvaluate_HalfDistance,
    Type_PlacingActionEvaluate_Angle,
    Type_PlacingActionEvaluate_Torsion
} Type_PlacingActionEvalParam;

typedef struct _CheckMultiConsStep{
	int designSite;
	int countOfConsToCheckAtThisStep;
	int consToCheckAtThisStep[MAX_COUNT_CHECK_MULTI_CONS];
} CheckMultiConsStep;

typedef struct _PlacingAction{
    Type_PlacingAction actionType;

    //Fields if actionType is "LoadAtom"
    IntArray load_atoms;


    //Fields if actionType is "VariateParam" 
    int variate_param;
    double variate_from;
    double variate_to;
    double variate_increment;

    //Fields if actionType is "CalcAtom" 
    IntArray calc_atoms;
    IntArray calc_params;

    //Fields if actionType is "CheckCataCons" 
    char checkCataCons_firstSiteChainName[MAX_LENGTH_CHAIN_NAME+1];
    char checkCataCons_firstSiteResidueName[MAX_LENGTH_RESIDUE_NAME+1];
    int  checkCataCons_firstSitePosInChain;
	char checkCataCons_secondSiteChainName[MAX_LENGTH_CHAIN_NAME+1];
	char checkCataCons_secondSiteResidueName[MAX_LENGTH_RESIDUE_NAME+1];
	int  checkCataCons_secondSitePosInChain;
	CataConsSitePair* checkCataCons_pCataCon;
	int checkCataCons_pSite[2];

	//Fields if actionType is "CheckMultiCons"
	int checkMultiCons_cataConsCount;
	StringArray checkMultiCons_firstSiteChainNames;
	StringArray checkMultiCons_firstSiteResidueNames;
	IntArray    checkMultiCons_firstSitePosInChains;
	StringArray checkMultiCons_secondSiteChainNames;
	StringArray checkMultiCons_secondSiteResidueNames;
	IntArray    checkMultiCons_secondSitePosInChains;
	CataConsSitePair* checkMultiCons_pCataCons[MAX_COUNT_CHECK_MULTI_CONS];
	int checkMultiCons_siteIndexes[MAX_COUNT_CHECK_MULTI_CONS][2];
	int checkMultiCons_rotamerCountOnEachSite[MAX_COUNT_CHECK_MULTI_CONS][2];
	BOOL** checkMultiCons_predeterminedConsRelations[MAX_COUNT_CHECK_MULTI_CONS];
	int checkMultiCons_stepCount;
	CheckMultiConsStep checkMultiCons_steps[MAX_COUNT_CHECK_MULTI_CONS*2];
	


    //Fields if actionType is "CheckVDW_backbone"
    double checkVDW_backbone_maxAllowed;
    BOOL checkVDW_backbone_withHydrogen;
    double checkVDW_backbone_activeRange;
    double checkVDW_backbone_totalPotential;
	IntArray checkVDW_backbone_smallmolAtomHasXyz;

	//Fields if actionType is "CheckVDW_internal"
	double checkVDW_internal_maxAllowed;
	BOOL checkVDW_internal_withHydrogen;
	double checkVDW_internal_totalPotential;
	int checkVDW_internal_smallMolAtomCount;
	BOOL** checkVDW_internal_smallMolAtom13bondedMatrix;
	IntArray checkVDW_internal_smallmolAtomHasXyz;

    //Fields if actionType is CheckRMSD 
	IntArray checkRMSD_smallmolAtomHasXyz;
    double checkRMSD_minDifference;
    BOOL checkRMSD_withHydrogen;

    //Fields if actionType is "EvaluateParam" 
    int evaluate_param;
    Type_PlacingActionEvalParam evaluate_type;
    IntArray evaluate_atoms;
   
} PlacingAction;
int    PlacingActionCreate(PlacingAction* pThis,Type_PlacingAction type);
int    PlacingActionDestroy(PlacingAction* pThis);
BOOL   PlacingActionValidParamName(char* paramName);
int    PlacingActionRead_CALC(PlacingAction* pThis,StringArray* pAtomNames,XYZArray* pAtomXYZs,
							  StringArray* pParamNames,DoubleArray* pParams,StringArray* pContent);
//The following functions begin with "PlacingActionRead" are all sub-routines of "PlacingRuleReadFile()" 
int    PlacingActionRead_CHECK_CATA_CONS(PlacingAction* pThis,StringArray* pContent);
int    PlacingActionRead_CHECK_MULTI_CONS(PlacingAction* pThis,StringArray* pContent);
int    PlacingActionRead_CHECK_RMSD(PlacingAction* pThis,StringArray* pContent);
int    PlacingActionRead_CHECK_VDW_BACKBONE(PlacingAction* pThis,StringArray* pContent);
int    PlacingActionRead_CHECK_VDW_INTERNAL(PlacingAction* pThis,StringArray* pContent);
int    PlacingActionRead_EVALUATE(PlacingAction* pThis,StringArray* pAtomNames,
								  StringArray* pParamNames,DoubleArray* pParams,
								  StringArray* pContent);
int    PlacingActionRead_LOAD(PlacingAction* pThis,StringArray* pAtomNames,
							  XYZArray* pAtomXYZs,StringArray* pContent);
int    PlacingActionRead_VARIATE(PlacingAction* pThis,StringArray* pParamNames,
								 DoubleArray* pParams,StringArray* pContent);



typedef struct _PlacingRule{

    char chainName[MAX_LENGTH_CHAIN_NAME+1];
    char residueName[MAX_LENGTH_RESIDUE_NAME+1];
    int  posInChain;
    char rotamerType[MAX_LENGTH_RESIDUE_NAME+1];

    StringArray atomNames;
    XYZArray    atomXYZs;

    StringArray paramNames;
    DoubleArray params;

    PlacingAction* actions;
    int actionCount;

	// these two fields have no meaning before the rule is deployed;
    IntArray    atomPosOnSmallMol;
    XYZArray    smallMolAtomXYZs;

	Residue* pSmallMol;
    AtomArray* pTruncatedBackbone;

	// a flag for recording if the rule has been deployed.
    BOOL deployedFlag;

	double vdwBackbone;
	double vdwInternal;


} PlacingRule;
int			PlacingRuleCreate(PlacingRule* pThis,char* fileName);
void		PlacingRuleDestroy(PlacingRule* pThis);
BOOL        PlacingRuleGetDeployedFlag(PlacingRule* pThis);
char* PlacingRuleGetResiName(PlacingRule* pThis);
int         PlacingRuleGetPosInChain(PlacingRule* pThis);
char* PlacingRuleGetChainName(PlacingRule* pThis);
char* PlacingRuleGetRotamerType(PlacingRule* pThis);
double      PlacingRuleGetTruncatedBackboneRange(PlacingRule* pThis);
//The following three functions are all huge
//PlacingRuleReadFile() will call PlacingActionRead_XXX to read every PlacingAction 
int			PlacingRuleReadFile(PlacingRule* pThis,FileReader* pFileReader);
//PlacingRuleDeploy() will call PlacingActionDeploy_XXX to deploy every PlacingAction. 
int         PlacingRuleDeploy(PlacingRule* pThis,
							  Residue* pSmallMol,
							  CataConsSitePairArray* pCataConsArray,
							  int relatedProteinSiteCount,
							  DesignSite** relatedProteinSites,
							  AtomArray* pTruncatedBackbone);
//The following functions begin with "PlacingActionDeploy_XXX" are all sub-routines of "PlacingRuleDepoly()" 
int PlacingActionDeploy_CheckCataCons(PlacingAction* pThis,
									  Residue* pSmallMol,
									  CataConsSitePairArray* pCataConsArray,
									  char* startingSiteChainName,
									  int startingSitePosInChain,
									  int relatedProteinSiteCount,
									  DesignSite* relatedProteinSites);
int PlacingActionDeploy_CheckMultiCons(PlacingAction* pThis,
									   Residue* pSmallMol,
									   CataConsSitePairArray* pCataConsArray,
									   char* startingSiteChainName,
									   char* startingSiteResidueName,
									   int startingSitePosInChain,
									   int relatedProteinSiteCount,
									   DesignSite* relatedProteinSites);
int PlacingActionDeploy_CheckRMSD(PlacingAction* pThis,IntArray* pSmallmolAtomGetXyzByThisStep);
int PlacingActionDeploy_CheckVDWBackbone(PlacingAction* pThis,IntArray* pSmallmolAtomGetXyzByThisStep);
int PlacingActionDeploy_CheckVDWInternal(PlacingAction* pThis,Residue* pSmallMol,
										 IntArray* pSmallmolAtomGetXyzByThisStep);
int PlacingActionDeploy_Calc(PlacingAction* pThis,
							 IntArray* pAtomPosOnSmallMol,
							 IntArray* pSmallmolAtomGetXyzByThisStep);

int PlacingActionDeploy_Load(PlacingAction* pThis,
							 IntArray* pAtomPosOnSmallMol,
							 IntArray* pSmallmolAtomGetXyzByThisStep);
//PlacingRuleProcess() will call the following sub-routines to execute every PlacingAction 
int			PlacingRuleProcess(	PlacingRule* pThis,
								Rotamer* pProteinRotamerOnStartingSite,
								CataConsSitePairArray* pCataConsArray,
								int relatedProteinSiteCount,
								DesignSite* relatedProteinSites,
								RotamerSet* pSmallMolRotSetForOutput,
								int step);
//The following functions begin with "PlacingRuleProcess_" are all sub-routines of "PlacingRuleProcess()" 
int    PlacingRuleProcess_CALC(PlacingRule* pThis,PlacingAction* pAction);
BOOL   PlacingRuleProcess_CHECK_CATA_CONS(PlacingRule* pThis,
										  Rotamer* pProteinRotamerOnStartingSite,										  
										  int relatedProteinSiteCount,
										  DesignSite* relatedProteinSites,
										  PlacingAction* pAction);
BOOL   PlacingRuleProcess_CHECK_MULTI_CONS_Recurrence(PlacingAction* pAction,
													  int relatedProteinSiteCount,
													  DesignSite** relatedProteinSites,
													  int* rotamersSelectedOnPreviousSteps,
													  Rotamer* pSmallMolRotamer,
													  Rotamer* pStartingProteinRotamer,
													  int currentStepIndex);
BOOL   PlacingRuleProcess_CHECK_MULTI_CONS( PlacingRule* pThis,
											Rotamer* pProteinRotamerOnStartingSite,										  
											int relatedProteinSiteCount,
											DesignSite** relatedProteinSites,
											PlacingAction* pAction);
BOOL   PlacingRuleProcess_CHECK_RMSD(PlacingRule* pThis,RotamerSet* pSmallMolRotSet,PlacingAction* pAction);
BOOL   PlacingRuleProcess_CHECK_VDW_BACKBONE(PlacingRule* pThis,PlacingAction* pAction);
BOOL   PlacingRuleProcess_CHECK_VDW_INTERNAL(PlacingRule* pThis,PlacingAction* pAction);
int    PlacingRuleProcess_EVALUATE(PlacingRule* pThis,PlacingAction* pAction);
int    PlacingRuleProcess_LOAD(PlacingRule* pThis,
							   Rotamer* pProteinRotamerOnStartingSite,
							   PlacingAction* pAction);

//PlacingRulePlaceSmallMol() is the ultimate user interface in order to place small molecules 
int    PlacingRulePlaceSmallMol(PlacingRule* pThis,
								CataConsSitePairArray* pCataConsArray,
								DesignSite* pStartingSite,
								int relatedProteinDesignSiteCount,
								DesignSite** relatedProteinDesignSites,
								RotamerSet* pSmallMolRotSetForOutput);

//Methods used for debugging 
int    PlacingRuleShowAtom(PlacingRule* pThis,int index);
int    PlacingRuleShowAtomDistances(PlacingRule* pThis);
int    PlacingRuleShowParam(PlacingRule* pThis,int index);
int    PlacingRuleShowAction(PlacingRule* pThis,int index);
int    PlacingRuleShow(PlacingRule* pThis);
int    PlacingRuleTester(char* filename);

//Supplementary 
int    ScreenSmallmolRotamersByRMSD(char* oriFileName,char* outputFileName,double rmsdThresold);

int    AnalyzeSmallMolRotamers(char* oriFileName,Residue* pOriginalSmallMol);
int    AnalyzeSmallMolRotamersForSpecifiedAtoms(char* oriFileName,Residue* pNativeSmallMol,char* specificAtomFile);

int CompareByInternalEnergy(const void *a, const void *b);
int CompareByBackboneEnergy(const void *a, const void *b);
int SmallMolRotamersGetBothHighRankOfBackboneVdwAndInternalVdw(char *oldRotamersFile,char* newRotamersFile,double percent);

#endif  //SMALL_MOL_H 
