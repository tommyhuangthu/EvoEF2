///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#include "SmallMol.h"
#include <string.h>

#include <ctype.h>
#include <time.h>

//#define DEBUGGING_SMALLMOL


///////////////////////////////////////////
//The Catalytic Constraints Hierarchy
///////////////////////////////////////////
//Item
int    CataConsItemShow(CataConsItem* pThis){
	int i;
	char type[MAX_LENGTH_ONE_LINE_IN_FILE+1];
	char location[MAX_LENGTH_ONE_LINE_IN_FILE+1];

	switch(pThis->type){
	case Type_CataCons_NewPseudoAtom:
		strcpy(type,"PSEUDO_ATOM");break;
	case Type_CataCons_Distance:
		strcpy(type,"DISTANCE   ");break;
	case Type_CataCons_Angle:
		strcpy(type,"ANGLE      ");break;
	case Type_CataCons_Torsion:
		strcpy(type,"TORSION    ");break;
	}
	printf("%s ",type);
	for(i=0;i<4;i++){
		printf("%4.4s",pThis->atomName[i]);
		switch(pThis->atomLocation[i]){
		case Type_CataConsAtomLocation_Pseudo:
			strcpy(location,"(PSU)");break;
		case Type_CataConsAtomLocation_FirstSite:
			strcpy(location,"(1st)");break;
		case Type_CataConsAtomLocation_SecondSite:
			strcpy(location,"(2nd)");break;
		default:
			strcpy(location,"(UKN)");break;
		}
		printf("%s ",location);
		printf("%2d ",pThis->atomIndex[i]);
	}
	printf("%.3f %.3f\n",pThis->min,pThis->max);
	return Success;
}

//Group
int    CataConsGroupDeploy(CataConsGroup* pThis,Rotamer* pFirstSiteRotamer,Rotamer* pSecondSiteRotamer){
	int itemIndex;
	char errMsg[MAX_LENGTH_ERR_MSG+1];
	int  result = Success;
	StringArray pseudoAtomNames;
	StringArrayCreate(&pseudoAtomNames);

	for(itemIndex=0; itemIndex<pThis->consItemsCount; itemIndex++){
		CataConsItem* pCurItem = &pThis->consItems[itemIndex];
		if(pCurItem->type==Type_CataCons_NewPseudoAtom){
			StringArrayAppend(&pseudoAtomNames,pCurItem->atomName[0]);
		}
	}

	for(itemIndex=0; itemIndex<pThis->consItemsCount; itemIndex++){
		int i;
		int atomCount;
		CataConsItem* pCurItem = &pThis->consItems[itemIndex];
		switch(pCurItem->type){
			case Type_CataCons_NewPseudoAtom:
				atomCount = 3;
				break;
			case Type_CataCons_Distance:
				atomCount = 2;
				break;
			case Type_CataCons_Angle:
				atomCount = 3;
				break;
			case Type_CataCons_Torsion:
				atomCount = 4;
				break;
			default:
				atomCount = 4;
		}

		for(i=0;i<atomCount;i++){

			switch(pCurItem->atomLocation[i]){
				case Type_CataConsAtomLocation_FirstSite:
					result = RotamerFindAtom(pFirstSiteRotamer,pCurItem->atomName[i],&pCurItem->atomIndex[i]);
					if( FAILED(result) ){
						sprintf(errMsg,"In CataConsGroupDeploy(), cannot find atom %s on the "
							"first site",pCurItem->atomName[i]);
						TraceError(errMsg,result);
						return result;
					}
					break;

				case Type_CataConsAtomLocation_SecondSite:
					result = RotamerFindAtom(pSecondSiteRotamer,pCurItem->atomName[i],&pCurItem->atomIndex[i]);
					if( FAILED(result) ){
						sprintf(errMsg,"In CataConsGroupDeploy(), cannot find atom %s on the "
							"second site",pCurItem->atomName[i]);
						TraceError(errMsg,result);
						return result;
					}
					break;

				case Type_CataConsAtomLocation_Pseudo:
					result = StringArrayFind(&pseudoAtomNames,pCurItem->atomName[i],&pCurItem->atomIndex[i]);
					if( FAILED(result) ){
						sprintf(errMsg,"In CataConsGroupDeploy(), cannot find atom %s in the pseudo atom set",
							pCurItem->atomName[i]);
						TraceError(errMsg,result);
						return result;
					}
					break;

				default:
					sprintf(errMsg,"In CataConsGroupDeploy(), unspecified location for atom %s",
						pCurItem->atomName[i]);
					TraceError(errMsg,result);
					return result;
			}
		}
	}

	StringArrayDestroy(&pseudoAtomNames);

	return Success;
}

//Site
int         CataConsSitePairCreate(CataConsSitePair* pThis,FileReader* pFile){
	char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
	char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1];
	char errMsg[MAX_LENGTH_ERR_MSG+1];
	char atomName[MAX_LENGTH_ONE_LINE_IN_FILE+1];
	int  result;

	CataConsItem* pCurItem;
	CataConsGroup* pCurGroup;

	BOOL done;

	//the deployed flad is initialized as FALSE;
	pThis->deployedFlag = FALSE;

	pThis->groups = NULL;
	pThis->groupsCount = 0;

	done = FALSE;
	pCurGroup = NULL;
	while(!FAILED(FileReaderGetNextLine(pFile,line))){
		ExtractFirstStringFromSourceString(keyword,line);

		// when the keyword 'SITEPAIR' is first met, the value of flag is false;
		// and the information following the 'SITEPAIR' should be read, and the flag is set as true.
		// when the keyword 'SITEPAIR' is met again, the value of flag is true;
		// the position of filereader should be reset and exit this function;
		if(strcmp(keyword,"SITEPAIR")==0  && done==FALSE){  
				sscanf(line,"%s %s %d %s %s %d",
				pThis->firstSiteChainName,pThis->firstSiteResidueName,&pThis->firstSitePosInChain,
				pThis->secondSiteChainName,pThis->secondSiteResidueName,&pThis->secondSitePosInChain);
			done = TRUE;
		}
		else if(strcmp(keyword,"SITEPAIR")==0 && done==TRUE){
			FileReaderSetCurrentPos(pFile,FileReaderGetCurrentPos(pFile)-1);
			return Success;
		}
		else if(strcmp(keyword,"GROUP")==0){
			//Create New Group
			(pThis->groupsCount)++;
			pThis->groups = (CataConsGroup*)realloc(pThis->groups,sizeof(CataConsGroup)*pThis->groupsCount);
			pCurGroup = &pThis->groups[pThis->groupsCount-1];
			sscanf(line,"%s %s",pCurGroup->firstSiteRotamerType,pCurGroup->secondSiteRotamerType);
			pCurGroup->consItemsCount = 0;
			pCurGroup->consItems = NULL;
			pCurGroup->pseudoAtomsCount = 0;
			pCurGroup->pseudoAtoms = NULL;
		}
		else if(strcmp(keyword,"ITEM")==0){
			//Create New Item
			int atomIndex;
			char constraintType[MAX_LENGTH_ONE_LINE_IN_FILE+1];

			// allocate memories for new catalytic constraint items;
			(pCurGroup->consItemsCount)++;
			pCurGroup->consItems = (CataConsItem*)realloc(pCurGroup->consItems,
				sizeof(CataConsItem)*pCurGroup->consItemsCount);
			pCurItem = &pCurGroup->consItems[pCurGroup->consItemsCount-1];

			// read the catalytic constraint type;
			if(FAILED(ExtractFirstStringFromSourceString(constraintType,line))){
				result = FormatError;
				sprintf(errMsg,"In CataConsSitePairCreate(), unrecognized format:\n%s",line);
				TraceError(errMsg,result);
				return result;
			}

			if(strcmp(constraintType,"PSEUDO_ATOM")==0){
				pCurItem->type = Type_CataCons_NewPseudoAtom;
				(pCurGroup->pseudoAtomsCount)++;
				pCurGroup->pseudoAtoms = (XYZ*)realloc(pCurGroup->pseudoAtoms,
					sizeof(XYZ)*pCurGroup->pseudoAtomsCount);
			}
			else if(strcmp(constraintType,"DISTANCE")==0){
				pCurItem->type = Type_CataCons_Distance;
			}
			else if(strcmp(constraintType,"ANGLE")==0){
				pCurItem->type = Type_CataCons_Angle;
			}
			else if(strcmp(constraintType,"TORSION")==0){
				pCurItem->type = Type_CataCons_Torsion;
			}
			else{
				result = FormatError;
				sprintf(errMsg,"In CataConsSitePairCreate(), unrecognized keyword %s",constraintType);
				TraceError(errMsg,result);
				return result;
			}

			// read the four atoms for each atom;
			for(atomIndex=0; atomIndex<4; atomIndex++){
				if(FAILED(ExtractFirstStringFromSourceString(atomName,line))){
					result = FormatError;
					sprintf(errMsg,"In CataConsSitePairCreate(), unrecognized format:\n%s",line);
					TraceError(errMsg,result);
					return result;
				}
				if(atomName[0] == SECOND_CATACON_SITE_SYMBOL){
					pCurItem->atomLocation[atomIndex] = Type_CataConsAtomLocation_SecondSite;
					// remove the symbol '*' before the atom name;
					strcpy(pCurItem->atomName[atomIndex],atomName+1);
				}
				else if(atomName[0] == PSEUDO_ATOM_SYMBOL){
					pCurItem->atomLocation[atomIndex] = Type_CataConsAtomLocation_Pseudo;
					// remove the symbol '+' before the atom name;
					strcpy(pCurItem->atomName[atomIndex],atomName+1);
				}
				else{
					pCurItem->atomLocation[atomIndex] = Type_CataConsAtomLocation_FirstSite;
					strcpy(pCurItem->atomName[atomIndex],atomName);
				}
			}

			// read the range for constraint;
			sscanf(line,"%lf %lf",&pCurItem->min,&pCurItem->max);

			if( pCurItem->type == Type_CataCons_Angle ||
				pCurItem->type == Type_CataCons_Torsion){
					pCurItem->min = DegToRad(pCurItem->min);
					pCurItem->max = DegToRad(pCurItem->max);
			}

			if(pCurItem->type!=Type_CataCons_NewPseudoAtom){
				// set a disturbance for the upper and lower bound;
				pCurItem->min -= fabs(DISTURBANCE_IN_RANGE_CHECK);
				pCurItem->max += fabs(DISTURBANCE_IN_RANGE_CHECK);
			}

			// initialize the index for the related four atoms, and set the index as -1;
			for(atomIndex=0; atomIndex<4; atomIndex++){
				pCurItem->atomIndex[atomIndex] = -1;
			}
		}

		else{
      ;
		}

	}
	if(done==TRUE){
		return Success;
	}
	else{
		return ValueError;
	}

}
void        CataConsSitePairDestroy(CataConsSitePair* pThis){
	int i;
	for(i=0;i<pThis->groupsCount;i++){
		free(pThis->groups[i].consItems);
		free(pThis->groups[i].pseudoAtoms);
	}
	free(pThis->groups);
}

BOOL        CataConsSitePairGetDeployedFlag(CataConsSitePair* pThis){
	return pThis->deployedFlag;
}

int         CataConsSitePairDeploy(CataConsSitePair* pThis,RotamerSet* pFirstSiteRotSet,RotamerSet* pSecondSiteRotSet)
{
	int i;
	int    result;
	char errMsg[MAX_LENGTH_ERR_MSG+1];

	if(CataConsSitePairGetDeployedFlag(pThis)){
		result = AssertionError;
		TraceError("In CataConsSitePairDeploy(), This Catalytic Constraint Site Has Already "
			"Been Deployed.",result);
		return result;
	}


	for(i=0;i<pThis->groupsCount;i++){
		Rotamer* pFirstSiteRotamer = NULL;
		Rotamer* pSecondSiteRotamer = NULL;
		char* firstSiteRotamerType = pThis->groups[i].firstSiteRotamerType;
		char* secondSiteRotamerType = pThis->groups[i].secondSiteRotamerType;

		pFirstSiteRotamer = RotamerSetGetRepresentative(pFirstSiteRotSet,firstSiteRotamerType);
		pSecondSiteRotamer = RotamerSetGetRepresentative(pSecondSiteRotSet,secondSiteRotamerType);

		if(pFirstSiteRotamer==NULL){
			result = DataNotExistError;
			sprintf(errMsg,"In CataConsSitePairDeploy(), there are no %s type Rotamer in the rotamerSet\n"
				"of the first catalytic constraint site, catalytic constraint may not apply",firstSiteRotamerType);
			TraceError(errMsg,result);
			return result;

		}
		if(pSecondSiteRotamer==NULL){
			result = DataNotExistError;
			sprintf(errMsg,"In CataConsSitePairDeploy(), there are no %s type Rotamer in the rotamerSet\n"
				"of the second catalytic constraint site, catalytic constraint may not apply",secondSiteRotamerType);
			TraceError(errMsg,result);
			return result;
		}


		result = CataConsGroupDeploy(&pThis->groups[i],pFirstSiteRotamer,pSecondSiteRotamer);
		if(FAILED(result)){
			TraceError("In CataConsSitePairDeploy()",result);
			return result;
		}
	}

	pThis->deployedFlag = TRUE;

	return Success;
}
int         CataConsSitePairShow(CataConsSitePair* pThis){
	int groupIndex;
	int itemIndex;
	CataConsGroup* pCurGroup = NULL;

	printf("SITE    %s   %s   %d    %s   %s   %d\n",
		pThis->firstSiteChainName,pThis->firstSiteResidueName,pThis->firstSitePosInChain,
		pThis->secondSiteChainName,pThis->secondSiteResidueName,pThis->secondSitePosInChain);

	for(groupIndex=0; groupIndex<pThis->groupsCount; groupIndex++){
		pCurGroup = &pThis->groups[groupIndex];
		printf("GROUP    %s    %s\n",pCurGroup->firstSiteRotamerType,pCurGroup->secondSiteRotamerType);
		for(itemIndex=0; itemIndex<pCurGroup->consItemsCount; itemIndex++){
			CataConsItemShow( &pCurGroup->consItems[itemIndex]);
		}
	}
	return Success; 
}


//Array
int           CataConsSitePairArrayCreate(CataConsSitePairArray* pThis,char* cataConsFile){
	int    result;
	FileReader file;
	result = FileReaderCreate(&file,cataConsFile);
	if(FAILED(result)){
		TraceError("In CataConsSitePairArrayCreate()",result);
		return result;
	}
	pThis->count = 0;
	pThis->sites = NULL;

	while(TRUE){
		int newCount = pThis->count + 1;
		CataConsSitePair* pCurSite;

		pThis->sites = (CataConsSitePair*)realloc(pThis->sites,sizeof(CataConsSitePair)*newCount);
		pCurSite = &pThis->sites[newCount-1];

		result = CataConsSitePairCreate(pCurSite,&file);
		if(FAILED(result)){
			break;
		}
		else{
			pThis->count = newCount;
		}
	}
	FileReaderDestroy(&file);
	return Success;
}
void          CataConsSitePairArrayDestroy(CataConsSitePairArray* pThis){
	int i;
	for(i=0;i<pThis->count;i++){
		CataConsSitePairDestroy(&pThis->sites[i]);
	}
	free(pThis->sites);
}
int           CataConsSitePairArrayGetCount(CataConsSitePairArray* pThis){
	return pThis->count;
}
CataConsSitePair* CataConsSitePairArrayGet(CataConsSitePairArray* pThis,int index){
	if(index<0 || index>pThis->count){
		return NULL;
	}
	return &pThis->sites[index];
}
CataConsSitePair* CataConsSitePairArrayFind(CataConsSitePairArray* pThis,
									char* firstSiteChainName,
									int firstSitePosInChain,
									char* firstSiteResidueName,
									char* secondSiteChainName,
									int secondSitePosInChain,
									char* secondSiteResidueName)
{
	int i;
	for(i=0;i<pThis->count;i++){
		if(firstSitePosInChain != pThis->sites[i].firstSitePosInChain){
			continue;
		}
		if(secondSitePosInChain != pThis->sites[i].secondSitePosInChain){
			continue;
		}
		if(strcmp(firstSiteChainName,pThis->sites[i].firstSiteChainName)!=0){
			continue;
		}
		if(strcmp(secondSiteChainName,pThis->sites[i].secondSiteChainName)!=0){
			continue;
		}
		if(strcmp(firstSiteResidueName,pThis->sites[i].firstSiteResidueName)!=0){
			continue;
		}
		if(strcmp(secondSiteResidueName,pThis->sites[i].secondSiteResidueName)!=0){
			continue;
		}
		return &pThis->sites[i];
	}
	return NULL;
}
int           CataConsSitePairArrayShow(CataConsSitePairArray* pThis){
	int i;
	for(i=0;i<CataConsSitePairArrayGetCount(pThis);i++){
		CataConsSitePairShow(CataConsSitePairArrayGet(pThis,i));
		printf("\n");
	}
	return Success;
}
int           CataConsSitePairArrayTester(char* filename){
	int i;
	CataConsSitePairArray cataCons;
	CataConsSitePairArrayCreate(&cataCons,filename);
	CataConsSitePairArrayShow(&cataCons);

	for(i=0;i<100;i++){
		CataConsSitePairArrayDestroy(&cataCons);
		CataConsSitePairArrayCreate(&cataCons,filename);
	}

	CataConsSitePairArrayShow(&cataCons);
	CataConsSitePairArrayDestroy(&cataCons);
	return Success;
}


//Catalytic constraint checking operations
BOOL CataConsItemCheck(CataConsItem* pThis,XYZArray* pOnFirstSite,
					   XYZArray* pOnSecondSite,XYZ* pseudoAtoms)
{
	int i;
	double value;
	XYZ* pAtoms[4];
	Atom* pAtomOnProtein;
	Atom* pAtomOnRotamer;
	XYZ* pPseudoAtom;
	XYZ vectorI,vectorJ;

	pAtomOnProtein = NULL;
	pAtomOnRotamer = NULL;
	pPseudoAtom = NULL;

	for(i=0;i<4;i++){
		int atomIndex = pThis->atomIndex[i];
		switch(pThis->atomLocation[i]){
			case Type_CataConsAtomLocation_FirstSite :
				pAtoms[i] = XYZArrayGet(pOnFirstSite,atomIndex);
				break;
			case Type_CataConsAtomLocation_SecondSite:
				pAtoms[i] = XYZArrayGet(pOnSecondSite,atomIndex);
				break;
			case Type_CataConsAtomLocation_Pseudo:
				pPseudoAtom = &pseudoAtoms[atomIndex];
				pAtoms[i] = pPseudoAtom;
				break;
			case Type_CataConsAtomLocation_Undefined:
				break;
		}
	}

	switch(pThis->type){
		case Type_CataCons_NewPseudoAtom:
			pPseudoAtom->X = pThis->min*(pAtoms[1]->X) + pThis->max*(pAtoms[2]->X);
			pPseudoAtom->Y = pThis->min*(pAtoms[1]->Y) + pThis->max*(pAtoms[2]->Y);
			pPseudoAtom->Z = pThis->min*(pAtoms[1]->Z) + pThis->max*(pAtoms[2]->Z);
			return TRUE;

		case Type_CataCons_Distance:
			value = XYZDistance(pAtoms[1],pAtoms[0]);
			//printf("dist = %f\n", value);
			if(value > pThis->max    ||   value < pThis->min ){
				return FALSE;
			}
			return TRUE;

		case Type_CataCons_Angle:
			vectorI = XYZDifference(pAtoms[0],pAtoms[1]);
			vectorJ = XYZDifference(pAtoms[2],pAtoms[1]);
			value = XYZAngle(&vectorI,&vectorJ);
			//printf("angle = %f\n", RadToDeg(value));
			if(value>pThis->max || value<pThis->min ){
				return FALSE;
			}
			return TRUE;

		case Type_CataCons_Torsion:
			value = GetTorsionAngle(pAtoms[0],pAtoms[1],pAtoms[2],pAtoms[3]);
			//printf("torsion = %f\n", RadToDeg(value));
			return RadInRange(value,pThis->min,pThis->max);
	}

	return FALSE;
}

BOOL CataConsGroupCheck(CataConsGroup* pThis,Rotamer* pOnFirstSite,Rotamer* pOnSecondSite){
	int i;
	if( strcmp(pOnFirstSite->type, pThis->firstSiteRotamerType )!=0 ||
		strcmp(pOnSecondSite->type,pThis->secondSiteRotamerType)!=0 )
	{
		return FALSE;
	}

	for(i=0;i<pThis->consItemsCount;i++){
		if( ! CataConsItemCheck(&pThis->consItems[i], 
			&pOnFirstSite->xyzs,&pOnSecondSite->xyzs,pThis->pseudoAtoms)){
			return FALSE;
		}
	}
	return TRUE;
}

BOOL CataConsSitePairCheck(CataConsSitePair* pThis,Rotamer* pOnFirstSite,Rotamer* pOnSecondSite){
	int i;

	if(pThis->deployedFlag == FALSE){
		TraceError("In CataConsSitePairCheck(), the CataConsSitePair has not been deployed yet",AssertionError);
		return FALSE;
	}

	for(i=0;i<pThis->groupsCount;i++){
		if( CataConsGroupCheck(&pThis->groups[i],pOnFirstSite,pOnSecondSite) ){
			return TRUE;
		}
	}
	return FALSE;
}


//////////////////////////////////////////////////////////////////////////
//Methods of Placing Rule                                             
//////////////////////////////////////////////////////////////////////////
//Various Actions In Placing SmallMol
int    PlacingActionCreate(PlacingAction* pThis,Type_PlacingAction type){
	int i;
	pThis->actionType = type;
	switch(type){
		case Type_PlacingAction_Load:
			IntArrayCreate(&pThis->load_atoms,0);
			break;
		case Type_PlacingAction_Evaluate:
			IntArrayCreate(&pThis->evaluate_atoms,0);
			break;
		case Type_PlacingAction_Calc:
			IntArrayCreate(&pThis->calc_atoms,0);
			IntArrayCreate(&pThis->calc_params,0);
			break;
		case Type_PlacingAction_CheckVDW_Backbone:
			IntArrayCreate(&pThis->checkVDW_backbone_smallmolAtomHasXyz,0);
			break;
		case Type_PlacingAction_CheckVDW_Internal:
			pThis->checkVDW_internal_smallMolAtomCount = 0;
			pThis->checkVDW_internal_smallMolAtom13bondedMatrix = NULL;
			IntArrayCreate(&pThis->checkVDW_internal_smallmolAtomHasXyz,0);
			break;
		case Type_PlacingAction_CheckRMSD:
			IntArrayCreate(&pThis->checkRMSD_smallmolAtomHasXyz,0);
			break;
		case Type_PlacingAction_CheckMultiCons:
			pThis->checkMultiCons_cataConsCount = 0;
			StringArrayCreate(&pThis->checkMultiCons_firstSiteChainNames);
			StringArrayCreate(&pThis->checkMultiCons_firstSiteResidueNames);
			IntArrayCreate(&pThis->checkMultiCons_firstSitePosInChains,0);
			StringArrayCreate(&pThis->checkMultiCons_secondSiteChainNames);
			StringArrayCreate(&pThis->checkMultiCons_secondSiteResidueNames);
			IntArrayCreate(&pThis->checkMultiCons_secondSitePosInChains,0);

			for(i=0;i<MAX_COUNT_CHECK_MULTI_CONS;i++){
				pThis->checkMultiCons_rotamerCountOnEachSite[i][0] = 
					pThis->checkMultiCons_rotamerCountOnEachSite[i][1] = 0;
				pThis->checkMultiCons_predeterminedConsRelations[i] = NULL;
			}

			break;
		default:
			break;
	}
	return Success;
}
int    PlacingActionDestroy(PlacingAction* pThis){
	int i;
	switch(pThis->actionType){
		case Type_PlacingAction_Load:
			IntArrayDestroy(&pThis->load_atoms);
			break;
		case Type_PlacingAction_Evaluate:
			IntArrayDestroy(&pThis->evaluate_atoms);
			break;
		case Type_PlacingAction_Calc:
			IntArrayDestroy(&pThis->calc_atoms);
			IntArrayDestroy(&pThis->calc_params);
			break;
		case Type_PlacingAction_CheckVDW_Backbone:
			IntArrayDestroy(&pThis->checkVDW_backbone_smallmolAtomHasXyz);
			break;
		case Type_PlacingAction_CheckVDW_Internal:
			for(i=0;i<pThis->checkVDW_internal_smallMolAtomCount;i++){
				free(pThis->checkVDW_internal_smallMolAtom13bondedMatrix[i]);
			}
			free(pThis->checkVDW_internal_smallMolAtom13bondedMatrix);
			pThis->checkVDW_internal_smallMolAtom13bondedMatrix = NULL;
			pThis->checkVDW_internal_smallMolAtomCount = 0;
			IntArrayDestroy(&pThis->checkVDW_internal_smallmolAtomHasXyz);
			break;
		case Type_PlacingAction_CheckRMSD:
			IntArrayDestroy(&pThis->checkRMSD_smallmolAtomHasXyz);
			break;
		case Type_PlacingAction_CheckMultiCons:
			pThis->checkMultiCons_cataConsCount = 0;
			StringArrayDestroy(&pThis->checkMultiCons_firstSiteChainNames);
			StringArrayDestroy(&pThis->checkMultiCons_firstSiteResidueNames);
			IntArrayDestroy(&pThis->checkMultiCons_firstSitePosInChains);
			StringArrayDestroy(&pThis->checkMultiCons_secondSiteChainNames);
			StringArrayDestroy(&pThis->checkMultiCons_secondSiteResidueNames);
			IntArrayDestroy(&pThis->checkMultiCons_secondSitePosInChains);
			for(i=0;i<MAX_COUNT_CHECK_MULTI_CONS;i++){
				int j;
				if(pThis->checkMultiCons_predeterminedConsRelations[i]!=NULL){
					for(j=0;j<pThis->checkMultiCons_rotamerCountOnEachSite[i][0];j++){
						free(pThis->checkMultiCons_predeterminedConsRelations[i][j]);
					}
					free(pThis->checkMultiCons_predeterminedConsRelations[i]);
					pThis->checkMultiCons_predeterminedConsRelations[i] = NULL;
				}
			}
			break;
		default:
			break;
	}
	return Success;
}


BOOL   PlacingActionValidParamName(char* paramName){
	if( isalpha(paramName[0]) || paramName[0]=='_'){
		return TRUE;
	}
	else{
		char errMsg[MAX_LENGTH_ERR_MSG+1];
		sprintf(errMsg,"Parameter name %s is invalid. It can only begin with a letter or an underscore",paramName);
		TraceError(errMsg,AssertionError);
		return FALSE;
	}
}

int    PlacingActionRead_CALC(PlacingAction* pThis,StringArray* pAtomNames,XYZArray* pAtomXYZs,
							  StringArray* pParamNames,DoubleArray* pParams,StringArray* pContent)
{
	int i;
	int index;
	char* atomName;
	char* paramName;
	char errMsg[MAX_LENGTH_ERR_MSG+1];
	int    result;

	if(StringArrayGetCount(pContent)!=9){
		return FormatError;
	}
	for(i=0;i<3;i++){
		atomName = StringArrayGet(pContent,i);
		//To skip the '*' in case the user uses '*' to denote an improper dihedral
		if(i==2 && atomName[0]=='*'){
			atomName++;
		}
		result = StringArrayFind(pAtomNames,atomName,&index);
		if(FAILED(result)){
			result = AssertionError;
			sprintf(errMsg,"In PlacingActionRead_CALC(), Atom %s does not exist or its coordinated has "
				"not been calculated yet",atomName);
			TraceError(errMsg,result);
			return result;
		}
		IntArrayAppend(&pThis->calc_atoms,index);
	}
	atomName = StringArrayGet(pContent,3);
	result = StringArrayFind(pAtomNames,atomName,&index);
	if(FAILED(result)){
		StringArrayAppend(pAtomNames,atomName);
		index = StringArrayGetCount(pAtomNames)-1;
		XYZArrayResize(pAtomXYZs,XYZArrayGetLength(pAtomXYZs)+1);
	}
	IntArrayAppend(&pThis->calc_atoms,index);

	for(i=0;i<5;i++){
		paramName = StringArrayGet(pContent,4+i);
		if(paramName[0]=='+' || paramName[0]=='-' || isdigit(paramName[0])){
			StringArrayAppend(pParamNames,"");
			if(1<=i && i<=3){
				DoubleArrayAppend(pParams,DegToRad(atof(paramName)));
			}
			else{
				DoubleArrayAppend(pParams,atof(paramName));
			}
			index = StringArrayGetCount(pParamNames) - 1;
			IntArrayAppend(&pThis->calc_params,index);
		}
		else{
			result = StringArrayFind(pParamNames,paramName,&index);
			if(FAILED(result)){
				sprintf(errMsg,"In PlacingActionRead_CALC(), parameter '%s'has not "
					"been defined",paramName);
				result = FormatError;
				TraceError(errMsg,result);
				return result;
			}
			IntArrayAppend(&pThis->calc_params,index);
		}
	}
	return Success;
}

int    PlacingActionRead_CHECK_CATA_CONS(PlacingAction* pThis,StringArray* pContent){
	if(StringArrayGetCount(pContent)!=6)
		return FormatError;

	strcpy(pThis->checkCataCons_firstSiteChainName,StringArrayGet(pContent,0));
	strcpy(pThis->checkCataCons_firstSiteResidueName,StringArrayGet(pContent,1));
	pThis->checkCataCons_firstSitePosInChain = atoi(StringArrayGet(pContent,2));

	strcpy(pThis->checkCataCons_secondSiteChainName,StringArrayGet(pContent,3));
	strcpy(pThis->checkCataCons_secondSiteResidueName,StringArrayGet(pContent,4));
	pThis->checkCataCons_secondSitePosInChain = atoi(StringArrayGet(pContent,5));

	return Success;
}
int    PlacingActionRead_CHECK_MULTI_CONS(PlacingAction* pThis,StringArray* pContent){
	int i;
	pThis->checkMultiCons_cataConsCount = StringArrayGetCount(pContent)/6;
	for(i=0;i<pThis->checkMultiCons_cataConsCount;i++){
		StringArrayAppend(&pThis->checkMultiCons_firstSiteChainNames,  StringArrayGet(pContent,6*i));
		StringArrayAppend(&pThis->checkMultiCons_firstSiteResidueNames,   StringArrayGet(pContent,6*i+1));
		IntArrayAppend(   &pThis->checkMultiCons_firstSitePosInChains, atoi(StringArrayGet(pContent,6*i+2)));
		StringArrayAppend(&pThis->checkMultiCons_secondSiteChainNames, StringArrayGet(pContent,6*i+3));
		StringArrayAppend(&pThis->checkMultiCons_secondSiteResidueNames,  StringArrayGet(pContent,6*i+4));
		IntArrayAppend(   &pThis->checkMultiCons_secondSitePosInChains,atoi(StringArrayGet(pContent,6*i+5)));
	}
	return Success;
}
int    PlacingActionRead_CHECK_RMSD(PlacingAction* pThis,StringArray* pContent){
	char* withHydrogen;
	char errMsg[MAX_LENGTH_ERR_MSG+1];
	int    result;

	if(StringArrayGetCount(pContent)!=2)
		return FormatError;

	withHydrogen = StringArrayGet(pContent,0);
	if(strcmp(withHydrogen,"WITH_HYDROGEN")==0){
		pThis->checkRMSD_withHydrogen = TRUE;
	}
	else if(strcmp(withHydrogen,"WITHOUT_HYDROGEN")==0){
		pThis->checkRMSD_withHydrogen = FALSE;
	}
	else{
		result = FormatError;
		sprintf(errMsg,"In PlacingActionRead_CHECK_RMSD(), Unrecognized keyword %s",withHydrogen);
		TraceError(errMsg,result);
		return result;
	}

	pThis->checkRMSD_minDifference = atof(StringArrayGet(pContent,1));
	return Success;
}

int    PlacingActionRead_CHECK_VDW_BACKBONE(PlacingAction* pThis,StringArray* pContent){
	char* withHydrogen;
	char errMsg[MAX_LENGTH_ERR_MSG+1];
	int    result;


	if(StringArrayGetCount(pContent)!=3)
		return FormatError;

	withHydrogen = StringArrayGet(pContent,0);
	if(strcmp(withHydrogen,"WITH_HYDROGEN")==0){
		pThis->checkVDW_backbone_withHydrogen = TRUE;
	}
	else if(strcmp(withHydrogen,"WITHOUT_HYDROGEN")==0){
		pThis->checkVDW_backbone_withHydrogen = FALSE;
	}
	else{
		result = FormatError;
		sprintf(errMsg,"In PlacingActionRead_CHECK_VDW(), Unrecognized keyword %s",withHydrogen);
		TraceError(errMsg,result);
		return result;
	}

	pThis->checkVDW_backbone_activeRange = atof(StringArrayGet(pContent,1));
	pThis->checkVDW_backbone_maxAllowed = atof(StringArrayGet(pContent,2));

	return Success;
}
int    PlacingActionRead_CHECK_VDW_INTERNAL(PlacingAction* pThis,StringArray* pContent){
	char* withHydrogen;
	char errMsg[MAX_LENGTH_ERR_MSG+1];
	int    result;

	if(StringArrayGetCount(pContent)!=2)
		return FormatError;

	withHydrogen = StringArrayGet(pContent,0);
	if(strcmp(withHydrogen,"WITH_HYDROGEN")==0){
		pThis->checkVDW_internal_withHydrogen = TRUE;
	}
	else if(strcmp(withHydrogen,"WITHOUT_HYDROGEN")==0){
		pThis->checkVDW_internal_withHydrogen = FALSE;
	}
	else{
		result = FormatError;
		sprintf(errMsg,"In PlacingActionRead_CHECK_VDW(), Unrecognized keyword %s",withHydrogen);
		TraceError(errMsg,result);
		return result;
	}

	pThis->checkVDW_internal_maxAllowed = atof(StringArrayGet(pContent,1));

	return Success;
}

int    PlacingActionRead_EVALUATE(PlacingAction* pThis,StringArray* pAtomNames,
								  StringArray* pParamNames,DoubleArray* pParams,
								  StringArray* pContent)
{
	int i;
	char* paramName;
	char* atomName;
	char* typeName;
	int wordCount;
	int atomCount;
	int    result;
	char errMsg[MAX_LENGTH_ERR_MSG+1];
	int index;

	wordCount = StringArrayGetCount(pContent);

	if(wordCount<4)
		return FormatError;

	paramName = StringArrayGet(pContent,0);
	if(!PlacingActionValidParamName(paramName)){
		return FormatError;
	}

	if(FAILED(StringArrayFind(pParamNames,paramName,&index))){
		StringArrayAppend(pParamNames,paramName);
		index = StringArrayGetCount(pParamNames)-1;
		DoubleArrayAppend(pParams,0.0);
		pThis->evaluate_param = index;
	}
	else{
		result = FormatError;
		sprintf(errMsg,"In PlacingActionRead_EVALUATE(), parameter %s has already been defined",paramName);
		TraceError(errMsg,result);
		return result;
	}

	typeName = StringArrayGet(pContent,1);

	atomCount = 0;
	if(strcmp(typeName,"HALF_DISTANCE")==0){
		pThis->evaluate_type = Type_PlacingActionEvaluate_HalfDistance;
		atomCount = 2;
	}
	else if(strcmp(typeName,"DISTANCE")==0){
		pThis->evaluate_type = Type_PlacingActionEvaluate_Distance;
		atomCount = 2;
	}
	else if(strcmp(typeName,"ANGLE")==0){
		pThis->evaluate_type = Type_PlacingActionEvaluate_Angle;
		atomCount = 3;
	}
	else if(strcmp(typeName,"TORSION")==0){
		pThis->evaluate_type = Type_PlacingActionEvaluate_Torsion;
		atomCount = 4;
	}
	else{
		result = FormatError;
		sprintf(errMsg,"In PlacingActionRead_EVALUATE(),keyword '%s' is not defined",typeName);
		TraceError(errMsg,result);
		return result;
	}

	if(atomCount+2 != wordCount){
		return FormatError;
	}

	for(i=0;i<atomCount;i++){
		atomName = StringArrayGet(pContent,i+2);
		if(FAILED(StringArrayFind(pAtomNames,atomName,&index))){
			sprintf(errMsg,"In PlacingActionRead_EVALUATE(), Atom %s has not been defined, or its coordinate has not been calculated yet",atomName);
			result = FormatError;
			TraceError(errMsg,result);
			return result;
		}
		IntArrayAppend(&pThis->evaluate_atoms,index);
	}

	return Success;
}
int    PlacingActionRead_LOAD(PlacingAction* pThis,StringArray* pAtomNames,
							  XYZArray* pAtomXYZs,StringArray* pContent)
{
	int i;
	int index;
	char* newAtomName;
	for(i=0; i<StringArrayGetCount(pContent); i++){
		newAtomName = StringArrayGet(pContent,i);
		if(strlen(newAtomName)>MAX_LENGTH_ATOM_NAME){
			int    result = NameError;
			char errMsg[MAX_LENGTH_ERR_MSG+1];
			sprintf(errMsg,"In PlacingActionRead_LOAD(), Atom name %s is too long",newAtomName);
			TraceError(errMsg,result);
			return result;
		}
		if( FAILED(StringArrayFind(pAtomNames,newAtomName,&index)) ){
			StringArrayAppend(pAtomNames,newAtomName);
			index = StringArrayGetCount(pAtomNames)-1;
			IntArrayAppend(&pThis->load_atoms,index);
			XYZArrayResize(pAtomXYZs,XYZArrayGetLength(pAtomXYZs)+1);
		}
	}
	return Success;
}

int    PlacingActionRead_VARIATE(PlacingAction* pThis,StringArray* pParamNames,
								 DoubleArray* pParams,StringArray* pContent)
{
	char* name;
	char* degreeOrAnsgrom;
	int index;
	char errMsg[MAX_LENGTH_ERR_MSG+1];
	int    result;

	if(StringArrayGetCount(pContent)!=5){
		return FormatError;
	}
	name = StringArrayGet(pContent,0);
	if(    !PlacingActionValidParamName(name) ){
		return FormatError;
	}

	degreeOrAnsgrom = StringArrayGet(pContent,1);
	if( FAILED(StringArrayFind(pParamNames,name,&index)) ){

		StringArrayAppend(pParamNames,name);
		DoubleArrayAppend(pParams,0.0);

		index = StringArrayGetCount(pParamNames)-1;
		pThis->variate_param = index;

		pThis->variate_from      = atof(StringArrayGet(pContent,2));
		pThis->variate_to        = atof(StringArrayGet(pContent,3));
		pThis->variate_increment = atof(StringArrayGet(pContent,4));

		//Add a disturbance, to include the boundary into the interval
		pThis->variate_to += DISTURBANCE_IN_RANGE_CHECK;

		if(strcmp(degreeOrAnsgrom,"DEGREE")==0){
			pThis->variate_from = DegToRad(pThis->variate_from);
			pThis->variate_to = DegToRad(pThis->variate_to);
			pThis->variate_increment = DegToRad(pThis->variate_increment);

			//For variation of degrees, a maximum range of 2*PI is imposed.
			if( pThis->variate_to > pThis->variate_from + 2*PI - DISTURBANCE_IN_RANGE_CHECK)
				pThis->variate_to = pThis->variate_from + 2*PI - DISTURBANCE_IN_RANGE_CHECK;
		}
		else if(strcmp(degreeOrAnsgrom,"ANGSTROM")==0){
			; //do nothing
		}
		else{
			result = FormatError;
			sprintf(errMsg,"In PlacingActionRead_VARIATE(), invalid keyword %s",degreeOrAnsgrom);
			TraceError(errMsg,result);
			return result;
		}
	}
	else{
		result = AssertionError;
		sprintf("In PlacingActionRead_VARIATE(), parameter %s has already been defined",name);
		TraceError(errMsg,result);
		return result;
	}
	return Success;
}


int    PlacingRuleReadFile(PlacingRule* pThis,FileReader* pFileReader){
	BOOL done;
	char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
	char errMsg[MAX_LENGTH_ERR_MSG+1];
	int    result;

	done = FALSE;
	while( ! FAILED(FileReaderGetNextLine(pFileReader,line)) ){
		char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1];
		char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
		PlacingAction* pCurAction;
		StringArray wordsInLine;
		StringArrayCreate(&wordsInLine);

		ExtractFirstStringFromSourceString(keyword,line);
		if(strcmp(keyword,"PLACING")==0   && done==FALSE ){
			ExtractFirstStringFromSourceString(pThis->chainName,line);
			ExtractFirstStringFromSourceString(pThis->residueName,line);
			ExtractFirstStringFromSourceString(buffer,line);
			pThis->posInChain = atoi(buffer);
			ExtractFirstStringFromSourceString(pThis->rotamerType,line);
			done = TRUE;
			StringArrayDestroy(&wordsInLine);
			continue;  //continue 'while'
		}
		else if(strcmp(keyword,"PLACING")==0   && done==TRUE ){
			FileReaderSetCurrentPos(pFileReader,FileReaderGetCurrentPos(pFileReader)-1);
			StringArrayDestroy(&wordsInLine);
			return Success;
		}


		(pThis->actionCount)++;
		pThis->actions = (PlacingAction*)realloc(pThis->actions,
			sizeof(PlacingAction)*pThis->actionCount);
		pCurAction = & pThis->actions[pThis->actionCount-1];

		StringArraySplitString(&wordsInLine,line,' ');        

		if(strcmp(keyword,"LOAD")==0){
			PlacingActionCreate(pCurAction,Type_PlacingAction_Load);
			result = PlacingActionRead_LOAD(pCurAction,&pThis->atomNames,&pThis->atomXYZs,&wordsInLine);
		}
		else if(strcmp(keyword,"VARIATE")==0){
			PlacingActionCreate(pCurAction,Type_PlacingAction_Variate);
			result = PlacingActionRead_VARIATE(pCurAction,&pThis->paramNames,&pThis->params,&wordsInLine);
		}
		else if(strcmp(keyword,"CALC")==0){
			PlacingActionCreate(pCurAction,Type_PlacingAction_Calc);
			result = PlacingActionRead_CALC(pCurAction,&pThis->atomNames,&pThis->atomXYZs,
				&pThis->paramNames,&pThis->params,&wordsInLine);
		}
		else if(strcmp(keyword,"EVALUATE")==0){
			PlacingActionCreate(pCurAction,Type_PlacingAction_Evaluate);
			result = PlacingActionRead_EVALUATE(pCurAction,&pThis->atomNames,&pThis->paramNames,
				&pThis->params,&wordsInLine);
		}
		else if(strcmp(keyword,"CHECK_CATA_CONS")==0){
			PlacingActionCreate(pCurAction,Type_PlacingAction_CheckCataCons);
			result = PlacingActionRead_CHECK_CATA_CONS(pCurAction,&wordsInLine);
		}
		else if(strcmp(keyword,"CHECK_MULTI_CONS_BEGIN")==0){
			StringArray allContent;
			StringArrayCreate(&allContent);
			PlacingActionCreate(pCurAction,Type_PlacingAction_CheckMultiCons);
			while(FileReaderGetNextLine(pFileReader,line)==Success){
				StringArraySplitString(&wordsInLine,line,' ');
				if(StringArrayGetCount(&wordsInLine)!=6){
					break;
				}
				else{
					int i;
					for(i=0;i<6;i++){
						StringArrayAppend(&allContent,StringArrayGet(&wordsInLine,i));
					}
				}
			}
			if(strcmp(line,"CHECK_MULTI_CONS_END")==0){
				int count = StringArrayGetCount(&allContent)/6;
				if(count>MAX_COUNT_CHECK_MULTI_CONS){
					sprintf(errMsg,"In PlacingRuleReadFile(), for a CHECK_MULTI_CONS operation, no more than %d\n"
						"constraints can be checked simultaneously, but the user has provided %d constraints\n"
						"You can modify the macro \"MAX_COUNT_CHECK_MULTI_CONS\" in the source code file\n"
						"\"SmallMol.h\" to allow for a bigger number of constraints",MAX_COUNT_CHECK_MULTI_CONS,count);
					result = IndexError;
					TraceError(errMsg,result);
					return result;
				}
				result = PlacingActionRead_CHECK_MULTI_CONS(pCurAction,&allContent);
			}
			else{
				result = FormatError;
			}
			StringArrayDestroy(&allContent);
		}
		else if(strcmp(keyword,"CHECK_VDW_BACKBONE")==0){
			PlacingActionCreate(pCurAction,Type_PlacingAction_CheckVDW_Backbone);
			result = PlacingActionRead_CHECK_VDW_BACKBONE(pCurAction,&wordsInLine);
		}
		else if(strcmp(keyword,"CHECK_VDW_INTERNAL")==0){
			PlacingActionCreate(pCurAction,Type_PlacingAction_CheckVDW_Internal);
			result = PlacingActionRead_CHECK_VDW_INTERNAL(pCurAction,&wordsInLine);
		}
		else if(strcmp(keyword,"CHECK_RMSD")==0){
			PlacingActionCreate(pCurAction,Type_PlacingAction_CheckRMSD);
			result = PlacingActionRead_CHECK_RMSD(pCurAction,&wordsInLine);
		}
		else{
			sprintf(errMsg,"In PlacingRuleReadFile(), unrecognized keyword %s",keyword);
			result = FormatError;
			TraceError(errMsg,result);
			return result;
		}

		if( FAILED(result) ){
			sprintf(errMsg,"In PlacingReadFile(), when reading : \n"
				"%s  %s",keyword,line);
			TraceError(errMsg,result);
			StringArrayDestroy(&wordsInLine);
			return result;
		}
		StringArrayDestroy(&wordsInLine);
	}

	return Success;
}


//////////////////////////////////////////////////////////////////////////
//Methods of PlacingRule                                              
//////////////////////////////////////////////////////////////////////////
//Creation and destruction
int    PlacingRuleCreate(PlacingRule* pThis,char* fileName){
	int    result;
	FileReader file;
	result = FileReaderCreate(&file,fileName);
	if(FAILED(result)){
		TraceError("In PlacingRuleCreate",result);
		return result;        
	}


	StringArrayCreate(&pThis->atomNames);
	XYZArrayCreate(&pThis->atomXYZs,0);

	IntArrayCreate(&pThis->atomPosOnSmallMol,0);
	XYZArrayCreate(&pThis->smallMolAtomXYZs,0);

	StringArrayCreate(&pThis->paramNames);
	DoubleArrayCreate(&pThis->params,0);

	pThis->actions = NULL;
	pThis->actionCount = 0;

	pThis->deployedFlag = FALSE;

	result = PlacingRuleReadFile(pThis,&file);
	if(FAILED(result)){
		TraceError("In PlacingRuleCreate",result);
		return result;
	}

	FileReaderDestroy(&file);

	return Success;
}

void   PlacingRuleDestroy(PlacingRule* pThis){
	int i;
	for(i=0;i<pThis->actionCount;i++){
		PlacingActionDestroy(&pThis->actions[i]);
	}

	StringArrayDestroy(&pThis->atomNames);
	XYZArrayDestroy(&pThis->atomXYZs);

	IntArrayDestroy(&pThis->atomPosOnSmallMol);
	XYZArrayDestroy(&pThis->smallMolAtomXYZs);

	StringArrayDestroy(&pThis->paramNames);
	DoubleArrayDestroy(&pThis->params);

	free(pThis->actions);
	pThis->actions = NULL;
	pThis->actionCount = 0;

}




//Methods about deployment of PlacingRule
BOOL        PlacingRuleGetDeployedFlag(PlacingRule* pThis){
	return pThis->deployedFlag;
}

char* PlacingRuleGetResiName(PlacingRule* pThis){
	return pThis->residueName;
}
int         PlacingRuleGetPosInChain(PlacingRule* pThis){
	return pThis->posInChain;
}
char* PlacingRuleGetChainName(PlacingRule* pThis){
	return pThis->chainName;
}
char* PlacingRuleGetRotamerType(PlacingRule* pThis){
	return pThis->rotamerType;
}
double      PlacingRuleGetTruncatedBackboneRange(PlacingRule* pThis){
	int i;
	for(i=0;i<pThis->actionCount;i++){
		if(pThis->actions[i].actionType == Type_PlacingAction_CheckVDW_Backbone){
			return pThis->actions[i].checkVDW_backbone_activeRange;
		}
	}
	return 0.0;
}


int PlacingActionDeploy_CheckCataCons(PlacingAction* pThis,
									  Residue* pSmallMol,
									  CataConsSitePairArray* pCataConsArray,
									  char* startingSiteChainName,
									  int startingSitePosInChain,
									  int relatedProteinSiteCount,
									  DesignSite** relatedProteinSites)
{
	int i,j;
	pThis->checkCataCons_pCataCon = CataConsSitePairArrayFind(
		pCataConsArray,
		pThis->checkCataCons_firstSiteChainName,
		pThis->checkCataCons_firstSitePosInChain,
		pThis->checkCataCons_firstSiteResidueName,
		pThis->checkCataCons_secondSiteChainName,
		pThis->checkCataCons_secondSitePosInChain,
		pThis->checkCataCons_secondSiteResidueName);

	if(pThis->checkCataCons_pCataCon == NULL){
		char errMsg[MAX_LENGTH_ERR_MSG+1];
		sprintf(errMsg,"In PlacingActionDeploy_CheckCataCons(), "
			"Cannot find the catalytic constraint between site %s %d %s and site %s %d %s",
			pThis->checkCataCons_firstSiteChainName,
			pThis->checkCataCons_firstSitePosInChain,
			pThis->checkCataCons_firstSiteResidueName,
			pThis->checkCataCons_secondSiteChainName,
			pThis->checkCataCons_secondSitePosInChain,
			pThis->checkCataCons_secondSiteResidueName);
		TraceError(errMsg,DataNotExistError);
		return DataNotExistError;
	}


	//Determine the identity of the first and second site of the catalytic constraint
	for(i=0;i<2;i++){
		int posInChain;
		char* chainName;
		if(i==0){
			posInChain = pThis->checkCataCons_pCataCon->firstSitePosInChain;
			chainName = pThis->checkCataCons_pCataCon->firstSiteChainName;
		}
		else{
			posInChain = pThis->checkCataCons_pCataCon->secondSitePosInChain;
			chainName = pThis->checkCataCons_pCataCon->secondSiteChainName;
		}

		pThis->checkCataCons_pSite[i] = -1;

		if( posInChain==ResidueGetPosInChain(pSmallMol) &&
			strcmp(chainName,ResidueGetChainName(pSmallMol))==0){
				pThis->checkCataCons_pSite[i] = relatedProteinSiteCount;  //the small molecule
		}
		else if( 
			posInChain==startingSitePosInChain &&
			strcmp(chainName,startingSiteChainName)==0){
				pThis->checkCataCons_pSite[i] = relatedProteinSiteCount+1;  //the starting site
		}
		else{
			for(j=0;j<relatedProteinSiteCount;j++){
        //this is probably a smallmol design site, not finished, skip it
        if(relatedProteinSites[j] == NULL) continue;
				if( posInChain == DesignSiteGetPosInChain(relatedProteinSites[j]) &&
					strcmp(chainName,DesignSiteGetChainName(relatedProteinSites[j]))==0){
						pThis->checkCataCons_pSite[i] = j;
						break;
				}
			}
		}

		if(pThis->checkCataCons_pSite[i]==-1){
			char errMsg[MAX_LENGTH_ERR_MSG+1];
			sprintf(errMsg,"In PlacingActionDeploy_CheckCataCons(), "
				"Cannot find design site %s %d needed by checking of catalytic constraint",
				chainName,posInChain);
			TraceError(errMsg,DataNotExistError);
			return DataNotExistError;
		}
	}

	return Success;
}
int PlacingActionDeploy_CheckMultiCons(PlacingAction* pThis,
									  Residue* pSmallMol,
									  CataConsSitePairArray* pCataConsArray,
									  char* startingSiteChainName,
									  char* startingSiteResidueName,
									  int startingSitePosInChain,
									  int relatedProteinSiteCount,
									  DesignSite** relatedProteinSites)
{
	int i,j;
	//Find Identity of every CataConsPair
	for(i=0;i<pThis->checkMultiCons_cataConsCount;i++){
		pThis->checkMultiCons_pCataCons[i] = CataConsSitePairArrayFind(
			pCataConsArray,
			StringArrayGet(&pThis->checkMultiCons_firstSiteChainNames,i),
			IntArrayGet(&pThis->checkMultiCons_firstSitePosInChains,i),
			StringArrayGet(&pThis->checkMultiCons_firstSiteResidueNames,i),
			StringArrayGet(&pThis->checkMultiCons_secondSiteChainNames,i),
			IntArrayGet(&pThis->checkMultiCons_secondSitePosInChains,i),
			StringArrayGet(&pThis->checkMultiCons_secondSiteResidueNames,i));

		if(pThis->checkMultiCons_pCataCons[i] == NULL){
			char errMsg[MAX_LENGTH_ERR_MSG+1];
			sprintf(errMsg,"In PlacingActionDeploy_CheckMultiCons(), "
				"Cannot find the catalytic constraint between site %s %d %s and site %s %d %s",
				StringArrayGet(&pThis->checkMultiCons_firstSiteChainNames,i),
				IntArrayGet(&pThis->checkMultiCons_firstSitePosInChains,i),
				StringArrayGet(&pThis->checkMultiCons_firstSiteResidueNames,i),
				StringArrayGet(&pThis->checkMultiCons_secondSiteChainNames,i),
				IntArrayGet(&pThis->checkMultiCons_secondSitePosInChains,i),
				StringArrayGet(&pThis->checkMultiCons_secondSiteResidueNames,i));
			TraceError(errMsg,DataNotExistError);
			return DataNotExistError;
		}
	}
	//Determine the identity of the first and second site of every catalytic constraint
	for(i=0;i<pThis->checkMultiCons_cataConsCount;i++){
		for(j=0;j<2;j++){
			int posInChain;
			char* chainName;
			if(j==0){
				posInChain = pThis->checkMultiCons_pCataCons[i]->firstSitePosInChain;
				chainName = pThis->checkMultiCons_pCataCons[i]->firstSiteChainName;
			}
			else{
				posInChain = pThis->checkMultiCons_pCataCons[i]->secondSitePosInChain;
				chainName = pThis->checkMultiCons_pCataCons[i]->secondSiteChainName;
			}

			pThis->checkMultiCons_siteIndexes[i][j] = -1;

			if( posInChain==ResidueGetPosInChain(pSmallMol) &&
				strcmp(chainName,ResidueGetChainName(pSmallMol))==0){
					pThis->checkMultiCons_siteIndexes[i][j] = relatedProteinSiteCount;
			}
			else if( 
				posInChain==startingSitePosInChain &&
				strcmp(chainName,startingSiteChainName)==0){
					pThis->checkMultiCons_siteIndexes[i][j] = relatedProteinSiteCount+1;
			}
			else{
				int s;
				for(s=0;s<relatedProteinSiteCount;s++){
          //this is probably a smallmol design site, not finished, skip it
          if(relatedProteinSites[s] == NULL) continue;
					if( posInChain == DesignSiteGetPosInChain(relatedProteinSites[s]) &&
						strcmp(chainName,DesignSiteGetChainName(relatedProteinSites[s]))==0){
							pThis->checkMultiCons_siteIndexes[i][j] = s;
							break;
					}
				}
			}

			if(pThis->checkMultiCons_siteIndexes[i][j]==-1){
				char errMsg[MAX_LENGTH_ERR_MSG+1];
				sprintf(errMsg,"In In PlacingActionDeploy_CheckMultiCons(), "
					"Cannot find design site %s %d needed by checking of catalytic constraint",
					chainName,posInChain);
				TraceError(errMsg,DataNotExistError);
				return DataNotExistError;
			}
		}
	}

	//Determine whether the catalytic constraints between only protein sites are satisfied
	for(i=0;i<pThis->checkMultiCons_cataConsCount;i++){
		int indexOfSite1 = pThis->checkMultiCons_siteIndexes[i][0];
		int indexOfSite2 = pThis->checkMultiCons_siteIndexes[i][1];
		if( indexOfSite1<relatedProteinSiteCount && indexOfSite2<relatedProteinSiteCount)
		{
			int rotOnSite1;
			int rotOnSite2;
			int rotCountOnSite1;
			int rotCountOnSite2;
			pThis->checkMultiCons_rotamerCountOnEachSite[i][0] = rotCountOnSite1 = 
				RotamerSetGetCount(DesignSiteGetRotamers(relatedProteinSites[indexOfSite1]));
			pThis->checkMultiCons_rotamerCountOnEachSite[i][1] = rotCountOnSite2 =
				RotamerSetGetCount(DesignSiteGetRotamers(relatedProteinSites[indexOfSite2]));
			pThis->checkMultiCons_predeterminedConsRelations[i] = (BOOL**)malloc(sizeof(BOOL*)*rotCountOnSite1);
			for(rotOnSite1=0;rotOnSite1<rotCountOnSite1;rotOnSite1++){
				pThis->checkMultiCons_predeterminedConsRelations[i][rotOnSite1] = 
					(BOOL*)malloc(sizeof(BOOL)*rotCountOnSite2);
				for(rotOnSite2=0; rotOnSite2<rotCountOnSite2; rotOnSite2++){
					BOOL checkResult;
					checkResult = CataConsSitePairCheck(pThis->checkMultiCons_pCataCons[i],
						RotamerSetGet(DesignSiteGetRotamers(relatedProteinSites[indexOfSite1]),rotOnSite1),
						RotamerSetGet(DesignSiteGetRotamers(relatedProteinSites[indexOfSite2]),rotOnSite2)
						);
					pThis->checkMultiCons_predeterminedConsRelations[i][rotOnSite1][rotOnSite2] = checkResult;
				}
			}
		}
	}

	//Construct the data structure needed by the recurrence algorithm in PlacingRuleProcess_CHECK_MULTI_CONS()
	pThis->checkMultiCons_stepCount = 0;
	for(i=0;i<pThis->checkMultiCons_cataConsCount;i++){
		for(j=0;j<2;j++){
			int k;
			int siteToFindInPreviousSteps = pThis->checkMultiCons_siteIndexes[i][j];
			BOOL addNewStep = TRUE;
			for(k=0;k<pThis->checkMultiCons_stepCount;k++){
				if(pThis->checkMultiCons_steps[k].designSite == siteToFindInPreviousSteps){
					addNewStep = FALSE;
					break;
				}
			}
			if(addNewStep==TRUE){
				pThis->checkMultiCons_steps[pThis->checkMultiCons_stepCount].designSite = siteToFindInPreviousSteps;
				pThis->checkMultiCons_steps[pThis->checkMultiCons_stepCount].countOfConsToCheckAtThisStep = 0;
				pThis->checkMultiCons_stepCount++;
			}
		}
	}
	for(i=0;i<pThis->checkMultiCons_cataConsCount;i++){
		BOOL site1Found = FALSE;
		BOOL site2Found = FALSE;
		int stepIndex = 0;
		CheckMultiConsStep* pCurStep = NULL;
		for(stepIndex=0; stepIndex<pThis->checkMultiCons_stepCount; stepIndex++){
			pCurStep = &pThis->checkMultiCons_steps[stepIndex];
			if(pCurStep->designSite == pThis->checkMultiCons_siteIndexes[i][0]){
				site1Found = TRUE;
			}
			if(pCurStep->designSite == pThis->checkMultiCons_siteIndexes[i][1]){
				site2Found = TRUE;
			}
			if(site1Found==TRUE && site2Found==TRUE){
				break;
			}
		}
		pCurStep->consToCheckAtThisStep[pCurStep->countOfConsToCheckAtThisStep] = i;
		pCurStep->countOfConsToCheckAtThisStep++;
	}
	
	return Success;
}
int PlacingActionDeploy_CheckRMSD(PlacingAction* pThis,IntArray* pSmallmolAtomGetXyzByThisStep){
	IntArrayCopy(&pThis->checkRMSD_smallmolAtomHasXyz,pSmallmolAtomGetXyzByThisStep);
	return Success;
}

int PlacingActionDeploy_CheckVDWBackbone(PlacingAction* pThis,IntArray* pSmallmolAtomGetXyzByThisStep){
	IntArrayCopy(&pThis->checkVDW_backbone_smallmolAtomHasXyz,pSmallmolAtomGetXyzByThisStep);
	return Success;
}
int PlacingActionDeploy_CheckVDWInternal(PlacingAction* pThis,Residue* pSmallMol,IntArray* pSmallmolAtomGetXyzByThisStep){
	int i,j,k;
	pThis->checkVDW_internal_smallMolAtomCount = ResidueGetAtomCount(pSmallMol);
	pThis->checkVDW_internal_smallMolAtom13bondedMatrix = 
		(BOOL**)malloc(sizeof(BOOL*)*ResidueGetAtomCount(pSmallMol));
	for(i=0;i<ResidueGetAtomCount(pSmallMol);i++){
		pThis->checkVDW_internal_smallMolAtom13bondedMatrix[i] = 
			(BOOL*)malloc(sizeof(BOOL)*ResidueGetAtomCount(pSmallMol));
		memset(pThis->checkVDW_internal_smallMolAtom13bondedMatrix[i],FALSE,
			sizeof(BOOL)*ResidueGetAtomCount(pSmallMol));
	}
	//Find 1-2 bond relation
	for(i=0;i<ResidueGetAtomCount(pSmallMol);i++){
		for(j=0;j<ResidueGetAtomCount(pSmallMol);j++){
			Atom* pAtomI = ResidueGetAtom(pSmallMol,i);
			Atom* pAtomJ = ResidueGetAtom(pSmallMol,j);
			if(BondSetFind(&pSmallMol->bonds,AtomGetName(pAtomI),AtomGetName(pAtomJ)) == Type_Bond_None ){
				pThis->checkVDW_internal_smallMolAtom13bondedMatrix[i][j] = FALSE;
			}
			else{
				pThis->checkVDW_internal_smallMolAtom13bondedMatrix[i][j] = TRUE;
			}
		}
	}
	//Find 1-3 bond relation
	for(i=0;i<ResidueGetAtomCount(pSmallMol);i++){
		for(j=0;j<ResidueGetAtomCount(pSmallMol);j++){
			for(k=0;k<ResidueGetAtomCount(pSmallMol);k++){
				Atom* pAtomI = ResidueGetAtom(pSmallMol,i);
				Atom* pAtomJ = ResidueGetAtom(pSmallMol,j);
				Atom* pAtomK = ResidueGetAtom(pSmallMol,k);
				if(
					(BondSetFind(&pSmallMol->bonds,AtomGetName(pAtomI),AtomGetName(pAtomK)) != Type_Bond_None) &&
					(BondSetFind(&pSmallMol->bonds,AtomGetName(pAtomK),AtomGetName(pAtomJ)) != Type_Bond_None)
					){
						pThis->checkVDW_internal_smallMolAtom13bondedMatrix[i][j] = TRUE;
				}
			}
		}
	}

	IntArrayCopy(&pThis->checkVDW_internal_smallmolAtomHasXyz,pSmallmolAtomGetXyzByThisStep);
	return Success;
}
int PlacingActionDeploy_Calc(PlacingAction* pThis,
							 IntArray* pAtomPosOnSmallMol,
							 IntArray* pSmallmolAtomGetXyzByThisStep){
	int atomIndexInPlacingRule = IntArrayGet(&pThis->calc_atoms,3);
	int atomIndexOnSmallmol = IntArrayGet(pAtomPosOnSmallMol,atomIndexInPlacingRule);
	if(atomIndexOnSmallmol!=-1){
		IntArraySet(pSmallmolAtomGetXyzByThisStep,atomIndexOnSmallmol,1);
	}
	return Success;
}

int PlacingActionDeploy_Load(PlacingAction* pThis,
							 IntArray* pAtomPosOnSmallMol,
							 IntArray* pSmallmolAtomGetXyzByThisStep)
{
	int i;
	for(i=0;i<IntArrayGetLength(&pThis->load_atoms);i++){
		int atomIndexInPlacingRule = IntArrayGet(&pThis->load_atoms,i);
		int atomIndexOnSmallmol = IntArrayGet(pAtomPosOnSmallMol,atomIndexInPlacingRule);
		if(atomIndexOnSmallmol!=-1){
			IntArraySet(pSmallmolAtomGetXyzByThisStep,atomIndexOnSmallmol,1);
		}
	}		
	return Success;
}

int PlacingRuleDeploy(PlacingRule* pThis,Residue* pSmallMol,CataConsSitePairArray* pCataConsArray,int relatedProteinSiteCount,DesignSite** relatedProteinSites,AtomArray* pTruncatedBackbone){
	int i;
	int j;
	int    result;
	char errMsg[MAX_LENGTH_ERR_MSG+1];
	IntArray smallmolAtomGetXyzByThisStep;


	//Make Sure That This Placing Rule Has Not Been Deployed Yet
	if(PlacingRuleGetDeployedFlag(pThis)){
		result = AssertionError;
		sprintf(errMsg,"In PlacingRuleDeploy(),Placing Rule %s %s %d has already been deployed",
			pThis->chainName,pThis->residueName,pThis->posInChain);
		TraceError(errMsg,result);
		return result;
	}else{
		pThis->deployedFlag = TRUE;
		pThis->pSmallMol = pSmallMol;
		pThis->pTruncatedBackbone = pTruncatedBackbone;
	}



	//Maintain the set of small molecule atom xyzs stored in this placing rule
	XYZArrayResize(&pThis->smallMolAtomXYZs,ResidueGetAtomCount(pSmallMol));

	//Find the exact locations of each Atom, special operations are needed if it's on the small molecule
	for(i=0;i<StringArrayGetCount(&pThis->atomNames);i++){
		int index;
		char atomName[MAX_LENGTH_ONE_LINE_IN_FILE];
		strcpy(atomName,StringArrayGet(&pThis->atomNames,i));
		if(atomName[0] == SMALLMOL_ATOM_SYMBOL){
			int result = ResidueFindAtom(pSmallMol,atomName+1,&index); 
			//Use atomName+1 to skip the '#' character used to represent small molecular atoms
			if(FAILED(result)){
				sprintf(errMsg,"In PlacingRuleDeploy(), Cannot find atom %s on the small molecule",atomName);
				TraceError(errMsg,result);
				return result;
			}
			IntArrayAppend(&pThis->atomPosOnSmallMol,index);
		}
		else{
			IntArrayAppend(&pThis->atomPosOnSmallMol,-1);
		}
	}


	IntArrayCreate(&smallmolAtomGetXyzByThisStep,XYZArrayGetLength(&pThis->smallMolAtomXYZs));
	for(j=0; j<IntArrayGetLength(&smallmolAtomGetXyzByThisStep); j++){
		IntArraySet(&smallmolAtomGetXyzByThisStep,j,0);
	}
	for(i=0;i<pThis->actionCount;i++){
		PlacingAction* pCurAction = &pThis->actions[i];
    //printf("deploy placing action %d\n", i);
		switch(pCurAction->actionType){
			case Type_PlacingAction_Calc:
				result = PlacingActionDeploy_Calc(pCurAction,&pThis->atomPosOnSmallMol,&smallmolAtomGetXyzByThisStep);
				break;
			case Type_PlacingAction_CheckCataCons:
				result = PlacingActionDeploy_CheckCataCons(pCurAction,pSmallMol,pCataConsArray,
					pThis->chainName,pThis->posInChain,relatedProteinSiteCount,relatedProteinSites);
				break;
			case Type_PlacingAction_CheckMultiCons:
				result = PlacingActionDeploy_CheckMultiCons(pCurAction,pSmallMol,pCataConsArray,pThis->chainName,
					pThis->residueName,pThis->posInChain,relatedProteinSiteCount,relatedProteinSites);
				break;
			case Type_PlacingAction_CheckRMSD:
				result = PlacingActionDeploy_CheckRMSD(pCurAction,&smallmolAtomGetXyzByThisStep);
				break;
			case Type_PlacingAction_CheckVDW_Internal:
				result = PlacingActionDeploy_CheckVDWInternal(pCurAction,pSmallMol,&smallmolAtomGetXyzByThisStep);
				break;
			case Type_PlacingAction_CheckVDW_Backbone:
				result = PlacingActionDeploy_CheckVDWBackbone(pCurAction,&smallmolAtomGetXyzByThisStep);
				break;
			case Type_PlacingAction_Load:
				result = PlacingActionDeploy_Load(pCurAction,&pThis->atomPosOnSmallMol,&smallmolAtomGetXyzByThisStep);
				break;
			default:
				result = Success;
				break;
		}
		if(FAILED(result)){
			sprintf(errMsg,"%s","In PlacingRuleDeploy(), Error occurred when deploy the following placing action");
			PlacingRuleShowAction(pThis,i);
			TraceError(errMsg,result);
			return result;
		}
	}

	IntArrayDestroy(&smallmolAtomGetXyzByThisStep);

	return Success;
}




//Methods about placing of small molecule
int    PlacingRuleProcess_CALC(PlacingRule* pThis,PlacingAction* pAction){
	XYZ atomXYZs[4];
	double icParam[5];
	int i;
	int atomIndex;
	int paramIndex;
	int atomPosOnSmallMol;
	for(i=0;i<3;i++){
		atomIndex = IntArrayGet(&pAction->calc_atoms,i);
		atomXYZs[i] = *XYZArrayGet(&pThis->atomXYZs,atomIndex);
	}
	for(i=0;i<5;i++){
		paramIndex = IntArrayGet(&pAction->calc_params,i);
		icParam[i] = DoubleArrayGet(&pThis->params,paramIndex);
	}

	GetFourthAtom(&atomXYZs[0],&atomXYZs[1],&atomXYZs[2],icParam,&atomXYZs[3]);
	atomIndex = IntArrayGet(&pAction->calc_atoms,3);
	XYZArraySet(&pThis->atomXYZs,atomIndex,&atomXYZs[3]);

	//If the atom is on small mol, additional operations must be taken to maintain consistency
	atomPosOnSmallMol = IntArrayGet(&pThis->atomPosOnSmallMol,atomIndex);
	if(atomPosOnSmallMol != -1){
		XYZArraySet(&pThis->smallMolAtomXYZs,atomPosOnSmallMol,&atomXYZs[3]);
	}
	return Success;
}
BOOL   PlacingRuleProcess_CHECK_CATA_CONS(PlacingRule* pThis,
										  Rotamer* pProteinRotamerOnStartingSite,										  
										  int relatedProteinSiteCount,
										  DesignSite** relatedProteinSites,
										  PlacingAction* pAction)

{

	int i,j;
	BOOL checkResult = FALSE;
	RotamerSet rotSetOnlySmallMol;
	RotamerSet rotSetOnlyStartingProteinRotamer;
	RotamerSet* pRotSetOnFirstAndSecondSite[2] = {NULL,NULL};
	RotamerSetCreate(&rotSetOnlySmallMol);
	RotamerSetCreate(&rotSetOnlyStartingProteinRotamer);

	RotamerSetAdd(&rotSetOnlyStartingProteinRotamer,pProteinRotamerOnStartingSite);
	if(pThis->pSmallMol != NULL){
		Rotamer tempRotamer;
		RotamerCreate(&tempRotamer);
		RotamerSetType(&tempRotamer,ResidueGetName(pThis->pSmallMol));
		RotamerAddAtoms(&tempRotamer,ResidueGetAllAtoms(pThis->pSmallMol));
		XYZArrayCopy(&tempRotamer.xyzs,&pThis->smallMolAtomXYZs);
		RotamerSetAdd(&rotSetOnlySmallMol,&tempRotamer);
		RotamerDestroy(&tempRotamer);
	}


	//Determine the identity of the first and second site of the catalytic constraint
	for(i=0;i<2;i++){
		int siteIndex = pAction->checkCataCons_pSite[i];
		if(siteIndex == relatedProteinSiteCount){
			pRotSetOnFirstAndSecondSite[i] = &rotSetOnlySmallMol;
		}
		else if(siteIndex == relatedProteinSiteCount+1){
			pRotSetOnFirstAndSecondSite[i] = &rotSetOnlyStartingProteinRotamer;
		}
		else{
			pRotSetOnFirstAndSecondSite[i] = DesignSiteGetRotamers(relatedProteinSites[siteIndex]);
		}
	}


	//check the constraints
	checkResult = FALSE;
	for(i=0;i<RotamerSetGetCount(pRotSetOnFirstAndSecondSite[0]);i++){
		if(checkResult==TRUE){
			break;
		}
		for(j=0;j<RotamerSetGetCount(pRotSetOnFirstAndSecondSite[1]);j++){
			checkResult = CataConsSitePairCheck(pAction->checkCataCons_pCataCon,
				RotamerSetGet(pRotSetOnFirstAndSecondSite[0],i),
				RotamerSetGet(pRotSetOnFirstAndSecondSite[1],j));
#ifdef DEBUGGING_SMALLMOL
			if(checkResult == TRUE){
				printf("(%s, %d) and (%s, %d) cons is true\n",
					RotamerGetType(RotamerSetGet(pRotSetOnFirstAndSecondSite[0],i)),
					i,
					RotamerGetType(RotamerSetGet(pRotSetOnFirstAndSecondSite[1],j)),
					j);
				break;
			}
#endif
		}
	}

	RotamerSetDestroy(&rotSetOnlySmallMol);
	RotamerSetDestroy(&rotSetOnlyStartingProteinRotamer);
	return checkResult;
}
BOOL   PlacingRuleProcess_CHECK_MULTI_CONS_Recurrence(PlacingAction* pAction,
													  int relatedProteinSiteCount,
													  DesignSite** relatedProteinSites,
													  int* rotamersSelectedOnPreviousSteps,
													  Rotamer* pSmallMolRotamer,
													  Rotamer* pStartingProteinRotamer,
													  int currentStepIndex)
{
	int i,j;
	int rotCount;
	CheckMultiConsStep* pCurStep = &pAction->checkMultiCons_steps[currentStepIndex];
	BOOL satisfied = FALSE;

	if(currentStepIndex == pAction->checkMultiCons_stepCount){	
		return TRUE;
	}

	if(pCurStep->designSite<relatedProteinSiteCount){
		rotCount = RotamerSetGetCount(DesignSiteGetRotamers(relatedProteinSites[pCurStep->designSite]));
	}
	else{
		rotCount = 1;
	}

	for(i=0;i<rotCount;i++){
		BOOL checkResult = TRUE;
		
		rotamersSelectedOnPreviousSteps[pCurStep->designSite] = i;

		for(j=0;j<pCurStep->countOfConsToCheckAtThisStep;j++){
			int consIndex = pCurStep->consToCheckAtThisStep[j];
			CataConsSitePair* pCurCons = pAction->checkMultiCons_pCataCons[consIndex];
			int site1Index = pAction->checkMultiCons_siteIndexes[consIndex][0];
			int site2Index = pAction->checkMultiCons_siteIndexes[consIndex][1];
			int rotOnSite1 = rotamersSelectedOnPreviousSteps[site1Index];
			int rotOnSite2 = rotamersSelectedOnPreviousSteps[site2Index];
			Rotamer* pRotOnSite1;
			Rotamer* pRotOnSite2;

			if(site1Index<relatedProteinSiteCount && site2Index<relatedProteinSiteCount){
				checkResult = pAction->checkMultiCons_predeterminedConsRelations[consIndex][rotOnSite1][rotOnSite2];
#ifdef DEBUGGING_SMALLMOL
				if(checkResult){
					printf("(%d %d, %d %d) satisfy constraint %d\n",site1Index,rotOnSite1,site2Index,rotOnSite2, consIndex);
				}
#endif
			}
			else{
				if(site1Index==relatedProteinSiteCount){
					pRotOnSite1 = pSmallMolRotamer;
				}
				else if(site1Index==relatedProteinSiteCount+1){
					pRotOnSite1 = pStartingProteinRotamer;
				}
				else{
					pRotOnSite1 = RotamerSetGet(DesignSiteGetRotamers(relatedProteinSites[site1Index]),rotOnSite1);
				}

				if(site2Index==relatedProteinSiteCount){
					pRotOnSite2 = pSmallMolRotamer;
				}
				else if(site2Index==relatedProteinSiteCount+1){
					pRotOnSite2 = pStartingProteinRotamer;
				}
				else{
					pRotOnSite2 = RotamerSetGet(DesignSiteGetRotamers(relatedProteinSites[site2Index]),rotOnSite2);
				}

				checkResult = CataConsSitePairCheck(pCurCons,pRotOnSite1,pRotOnSite2);
#ifdef DEBUGGING_SMALLMOL
        if(checkResult){
          if(site1Index==relatedProteinSiteCount+1){
            printf("(%d start, %d %d) satisfy constraint %d\n",site1Index,site2Index,rotOnSite2, consIndex);
          }
          else if(site2Index==relatedProteinSiteCount+1){
            printf("(%d %d, %d start) satisfy constraint %d\n",site1Index,rotOnSite1,site2Index, consIndex);
          }
          if(site1Index==relatedProteinSiteCount){
            printf("(%d ligand, %d %d) satisfy constraint %d\n",site1Index,site2Index,rotOnSite2, consIndex);
          }
          else if(site2Index==relatedProteinSiteCount){
            printf("(%d %d, %d ligand) satisfy constraint %d\n",site1Index,rotOnSite1,site2Index, consIndex);
          }
        }

#endif
			}

			if(checkResult==FALSE){
				break;
			}
		}

		if(checkResult==FALSE){
			continue;
		}


		checkResult = PlacingRuleProcess_CHECK_MULTI_CONS_Recurrence(pAction,relatedProteinSiteCount,
			relatedProteinSites,rotamersSelectedOnPreviousSteps,pSmallMolRotamer,pStartingProteinRotamer,
			currentStepIndex+1);
		if(checkResult==TRUE){
			satisfied = TRUE;
#ifdef DEBGGING_PlacingRuleProcess_CHECK_MULTI_CONS_Recurrence
			continue;
#else
			break;
#endif
		}
	}
	return satisfied;
#undef DEBGGING_PlacingRuleProcess_CHECK_MULTI_CONS_Recurrence
}
BOOL   PlacingRuleProcess_CHECK_MULTI_CONS(PlacingRule* pThis,
										   Rotamer* pProteinRotamerOnStartingSite,										  
										   int relatedProteinSiteCount,
										   DesignSite** relatedProteinSites,
										   PlacingAction* pAction)
{
	int i;
	int checkResult;
	int* rotamersSelectedOnPreviousSteps = (int*)malloc(sizeof(int)*(relatedProteinSiteCount+2));
	Rotamer smallMolRotamer;
	RotamerCreate(&smallMolRotamer);
	RotamerSetType(&smallMolRotamer,ResidueGetName(pThis->pSmallMol));
	RotamerAddAtoms(&smallMolRotamer,ResidueGetAllAtoms(pThis->pSmallMol));
	XYZArrayCopy(&smallMolRotamer.xyzs,&pThis->smallMolAtomXYZs);
	for(i=0;i<RotamerGetAtomCount(&smallMolRotamer);i++){
		RotamerGetAtom(&smallMolRotamer,i)->xyz = *XYZArrayGet(&pThis->smallMolAtomXYZs,i);
	}

	for(i=0;i<relatedProteinSiteCount+2;i++){
		rotamersSelectedOnPreviousSteps[i] = -1;
	}

	checkResult = PlacingRuleProcess_CHECK_MULTI_CONS_Recurrence(pAction,relatedProteinSiteCount,
		relatedProteinSites,rotamersSelectedOnPreviousSteps,&smallMolRotamer,pProteinRotamerOnStartingSite,0);

	RotamerDestroy(&smallMolRotamer);
	free(rotamersSelectedOnPreviousSteps);

	return checkResult;

}

BOOL   PlacingRuleProcess_CHECK_RMSD(PlacingRule* pThis,RotamerSet* pSmallMolRotSet,PlacingAction* pAction){
	int rotIndex;
	int atomIndex;
	int rotCount = RotamerSetGetCount(pSmallMolRotSet);
	int atomCount = ResidueGetAtomCount(pThis->pSmallMol);
	BOOL* omitThisAtom;

	if(rotCount==0){
		return TRUE;
	}

	omitThisAtom = (BOOL*)calloc(atomCount,sizeof(BOOL));
	for(atomIndex=0; atomIndex<atomCount; atomIndex++){
		if( AtomIsHydrogen(ResidueGetAtom(pThis->pSmallMol,atomIndex)) && 
			(!pAction->checkRMSD_withHydrogen) ){
				omitThisAtom[atomIndex] = TRUE;
		}
		else if(IntArrayGet(&pAction->checkRMSD_smallmolAtomHasXyz,atomIndex)==FALSE){
			omitThisAtom[atomIndex] = TRUE;
		}
		else{
			omitThisAtom[atomIndex] = FALSE;
		}
	}

	for(rotIndex=0; rotIndex<rotCount; rotIndex++){
		double totalRMSD = 0.0;
		int totalAtomCount = 0;
		Rotamer* pRot = RotamerSetGet(pSmallMolRotSet,rotIndex);
		for(atomIndex=0; atomIndex<atomCount; atomIndex++){
			XYZ* atomInRotSet;
			XYZ* atomInNewRot;
			if(omitThisAtom[atomIndex]){
				continue;
			}
			atomInRotSet = XYZArrayGet(&pRot->xyzs,atomIndex);
			atomInNewRot = XYZArrayGet(&pThis->smallMolAtomXYZs,atomIndex);
			totalRMSD += pow(XYZDistance(atomInRotSet,atomInNewRot),2);
			totalAtomCount++;
		}
		totalRMSD = sqrt(totalRMSD/totalAtomCount);
		//printf("%f ",totalRMSD);
		if(totalRMSD < pAction->checkRMSD_minDifference){
			free(omitThisAtom);
			return FALSE;
		}
	}

	free(omitThisAtom);
	return TRUE;
}
BOOL   PlacingRuleProcess_CHECK_VDW_BACKBONE(PlacingRule* pThis,PlacingAction* pAction){
	Atom* pAtomOnSmallMol;
	Atom* pAtomOnBackbone;
	int i;
	int j;

	pAction->checkVDW_backbone_totalPotential = 0.0;

	for(i=0;i<AtomArrayGetCount(pThis->pTruncatedBackbone);i++){
		pAtomOnBackbone = AtomArrayGet(pThis->pTruncatedBackbone,i);
		if(pAction->checkVDW_backbone_withHydrogen==FALSE && AtomIsHydrogen(pAtomOnBackbone)){
			continue;
		}
		for(j=0;j<ResidueGetAtomCount(pThis->pSmallMol);j++){
			double dist;
			double rminSum;
			double ratio;
			pAtomOnSmallMol = ResidueGetAtom(pThis->pSmallMol,j);
			if(pAction->checkVDW_backbone_withHydrogen == FALSE && AtomIsHydrogen(pAtomOnSmallMol)){
				continue;
			}
			if(IntArrayGet(&pAction->checkVDW_backbone_smallmolAtomHasXyz,j)==0){ 
				continue;
			}
			dist = XYZDistance(&(pAtomOnBackbone->xyz),XYZArrayGet(&pThis->smallMolAtomXYZs,j));
			rminSum = pAtomOnBackbone->vdw_radius + pAtomOnSmallMol->vdw_radius;

			ratio = dist/rminSum;
			if(ratio<1/1.1225){
				pAction->checkVDW_backbone_totalPotential += 10.0 - 11.225*ratio ;
			}
			pThis->vdwBackbone = pAction->checkVDW_backbone_totalPotential;
			
		}

		if(pAction->checkVDW_backbone_totalPotential > pAction->checkVDW_backbone_maxAllowed){
			return FALSE;
		}

	}

	return TRUE;
}
BOOL   PlacingRuleProcess_CHECK_VDW_INTERNAL(PlacingRule* pThis,PlacingAction* pAction){
	Atom* pAtomI;
	Atom* pAtomJ;
	int i;
	int j;

	pAction->checkVDW_internal_totalPotential = 0.0;

	for(i=0;i<pAction->checkVDW_internal_smallMolAtomCount;i++){
		pAtomI = ResidueGetAtom(pThis->pSmallMol,i);
		if(pAction->checkVDW_internal_withHydrogen==FALSE && AtomIsHydrogen(pAtomI)){
			continue;
		}
		if(IntArrayGet(&pAction->checkVDW_internal_smallmolAtomHasXyz,i)==0){
			//This atom's Xyz has not been calculated yet
			continue;
		}

		for(j=i+1;j<pAction->checkVDW_internal_smallMolAtomCount;j++){
			double dist;
			double rminSum;
			double ratio;
			pAtomJ = ResidueGetAtom(pThis->pSmallMol,j);
			if(pAction->checkVDW_internal_withHydrogen == FALSE && AtomIsHydrogen(pAtomJ)){
				continue;
			}
			if(IntArrayGet(&pAction->checkVDW_internal_smallmolAtomHasXyz,j)==0){
				//This atom's Xyz has not been calculated yet
				continue;
			}

			//If atomI and atomJ is bonded or 1-3 bonded, continue
			if(pAction->checkVDW_internal_smallMolAtom13bondedMatrix[i][j] == TRUE){
				continue;
			}
			//A rough estimation of VDW potential
			dist = XYZDistance(XYZArrayGet(&pThis->smallMolAtomXYZs,i),XYZArrayGet(&pThis->smallMolAtomXYZs,j));
			rminSum = pAtomI->vdw_radius + pAtomJ->vdw_radius;
			ratio = dist/rminSum;

			if(ratio < 1/1.1225){
				pAction->checkVDW_internal_totalPotential += 10.0 - 11.225*ratio;
			}
			//Debug
			pThis->vdwInternal = pAction->checkVDW_internal_totalPotential;

			if(pAction->checkVDW_internal_totalPotential > pAction->checkVDW_internal_maxAllowed){
				return FALSE;
			}
		}
	}
	return TRUE;
}
int    PlacingRuleProcess_EVALUATE(PlacingRule* pThis,PlacingAction* pAction){
	XYZ atomXYZs[4];
	int i;
	double newValue = 0.0;
	XYZ vectorI;
	XYZ vectorJ;

	for(i=0; i<IntArrayGetLength(&pAction->evaluate_atoms); i++){
		atomXYZs[i] = *XYZArrayGet(&pThis->atomXYZs,IntArrayGet(&pAction->evaluate_atoms,i));
	}

	switch(pAction->evaluate_type){
		case Type_PlacingActionEvaluate_HalfDistance:
			newValue = XYZDistance(&atomXYZs[0],&atomXYZs[1])*0.5;
			break;
		case Type_PlacingActionEvaluate_Distance:
			newValue = XYZDistance(&atomXYZs[0],&atomXYZs[1]);
			break;
		case Type_PlacingActionEvaluate_Angle:
			vectorI = XYZDifference(&atomXYZs[1],&atomXYZs[0]);
			vectorJ = XYZDifference(&atomXYZs[1],&atomXYZs[2]);
			newValue = XYZAngle(&vectorI,&vectorJ);
			break;
		case Type_PlacingActionEvaluate_Torsion:
			newValue = GetTorsionAngle(
				&atomXYZs[0],&atomXYZs[1],&atomXYZs[2],&atomXYZs[3]);
			break;
		default:
			//Should not get here
			break;
	}
	DoubleArraySet(&pThis->params,pAction->evaluate_param,newValue);
	return Success;    
}
int    PlacingRuleProcess_LOAD(PlacingRule* pThis,
							   Rotamer* pProteinRotamerOnStartingSite,
							   PlacingAction* pAction){
	int i;
	char errMsg[MAX_LENGTH_ERR_MSG+1];
	int    result;
	for(i=0;i<IntArrayGetLength(&pAction->load_atoms);i++){
		char* atomName;
		Atom* pAtomOnProteinOrSmallMol;
		int index = IntArrayGet(&pAction->load_atoms,i);
		atomName = StringArrayGet(&pThis->atomNames,index);

		if(atomName[0]==SMALLMOL_ATOM_SYMBOL){
			//Find Atom On SmallMol
			pAtomOnProteinOrSmallMol = ResidueGetAtomByName(pThis->pSmallMol,atomName+1);
			//use 'atomName+1' instead of 'atomName' to skip the '*' symbol
			if( pAtomOnProteinOrSmallMol != NULL ){
				int atomPosOnSmallMol = IntArrayGet(&pThis->atomPosOnSmallMol,index);
				XYZArraySet(&pThis->atomXYZs,index,&pAtomOnProteinOrSmallMol->xyz);
				//The atom is on small mol, additional operations must be taken to maintain consistency
				XYZArraySet(&pThis->smallMolAtomXYZs,atomPosOnSmallMol,&pAtomOnProteinOrSmallMol->xyz);  
				continue; //continue 'for'
			}
		}
		else{
			//Find Atom On Protein Rotamer
			Atom* pTempAtom = RotamerGetAtomByName(pProteinRotamerOnStartingSite,atomName);
			if(pTempAtom != NULL){
				XYZArraySet(&pThis->atomXYZs,index,&pTempAtom->xyz);
				continue; //continue 'for'
			}
		}

		//Cannot find the atom, report an error
		result = DataNotExistError;
		sprintf(errMsg,"In PlacingRuleProcess_LOAD(),cannot find Atom %s anywhere, when processing the"
			"'LOAD' action below :",atomName);
		TraceError(errMsg,result);
		PlacingRuleShowAction(pThis,(int)(pAction - pThis->actions));
		return result;
	}

	return Success;
}

int    PlacingRuleProcess(PlacingRule* pThis,
						  Rotamer* pProteinRotamerOnStartingSite,
						  CataConsSitePairArray* pCataConsArray,
						  int relatedProteinSiteCount,
						  DesignSite** relatedProteinSites,
						  RotamerSet* pSmallMolRotSetForOutput,
						  int step){

	int    result = Success;
	
	for( ; step<pThis->actionCount; step++){

		PlacingAction* pCurAction = &pThis->actions[step];
		BOOL cataConsSatisfied = FALSE;
		BOOL rmsdSatisfied = FALSE;
		BOOL vdwSatisfied = FALSE;

		if(pCurAction->actionType == Type_PlacingAction_Variate){
			break;
		}

		switch(pCurAction->actionType){
			case Type_PlacingAction_Calc:
				result = PlacingRuleProcess_CALC(pThis,pCurAction);
				break;
			case Type_PlacingAction_CheckCataCons:
				cataConsSatisfied = PlacingRuleProcess_CHECK_CATA_CONS(
					pThis,pProteinRotamerOnStartingSite,
					relatedProteinSiteCount,relatedProteinSites,pCurAction);
				result = Success;
				break;
			case Type_PlacingAction_CheckMultiCons:
				cataConsSatisfied = PlacingRuleProcess_CHECK_MULTI_CONS(
					pThis,pProteinRotamerOnStartingSite,
					relatedProteinSiteCount,relatedProteinSites,pCurAction);
				result = Success;
				break;
			case Type_PlacingAction_CheckRMSD:
				rmsdSatisfied = PlacingRuleProcess_CHECK_RMSD(pThis,pSmallMolRotSetForOutput,pCurAction);
				result = Success;
				break;
			case Type_PlacingAction_CheckVDW_Backbone:
				vdwSatisfied = PlacingRuleProcess_CHECK_VDW_BACKBONE(pThis,pCurAction);
				result = Success;
				break;
			case Type_PlacingAction_CheckVDW_Internal:
				vdwSatisfied = PlacingRuleProcess_CHECK_VDW_INTERNAL(pThis,pCurAction);
				result = Success;
				break;
			case Type_PlacingAction_Evaluate:
				result = PlacingRuleProcess_EVALUATE(pThis,pCurAction);
				break;
			case Type_PlacingAction_Load:
				result = PlacingRuleProcess_LOAD(pThis,pProteinRotamerOnStartingSite,pCurAction);
				break;
			default:
				printf("small molecule process type %d in placing rule can not be handled\n",pCurAction->actionType);
				result = FormatError;
				break;
		}

		if(FAILED(result)){
			return result;
		}
		if(pCurAction->actionType==Type_PlacingAction_CheckCataCons && !cataConsSatisfied){
			//printf("\nPlacing Step %d CataCons doesn't satisfy!\n", step+1);
			return Success;
		}
		if(pCurAction->actionType==Type_PlacingAction_CheckMultiCons && !cataConsSatisfied){
			//printf("\nPlacing Step %d MultiCons doesn't satisfy!\n", step+1);
			return Success;
		}
		if(pCurAction->actionType==Type_PlacingAction_CheckRMSD && !rmsdSatisfied){
			return Success;
		}
		if(pCurAction->actionType==Type_PlacingAction_CheckVDW_Backbone && !vdwSatisfied){
			return Success;
		}
		if(pCurAction->actionType==Type_PlacingAction_CheckVDW_Internal && !vdwSatisfied){
			return Success;
		}

	}

	if(step == pThis->actionCount){
		//Create a new rotamer, take care of atoms, bonds, Xyzs, and RBorns
		Rotamer newRot;
		RotamerCreate(&newRot);
		RotamerAddAtoms(&newRot,ResidueGetAllAtoms(pThis->pSmallMol));
		strcpy(newRot.type,ResidueGetName(pThis->pSmallMol));
		RotamerCopyAtomXYZ(&newRot,&pThis->smallMolAtomXYZs);
		BondSetCopy(&newRot.bonds, &pThis->pSmallMol->bonds );

		newRot.vdwBackbone = pThis->vdwBackbone;
		newRot.vdwInternal = pThis->vdwInternal;

		RotamerSetAdd(pSmallMolRotSetForOutput,&newRot);
		RotamerDestroy(&newRot);
		return Success;
	}

	if(pThis->actions[step].actionType == Type_PlacingAction_Variate){ 
		double curValue;
		for(curValue =  pThis->actions[step].variate_from; 
			curValue <  pThis->actions[step].variate_to;
			curValue += pThis->actions[step].variate_increment
			){
				DoubleArraySet(&pThis->params,pThis->actions[step].variate_param,curValue);
				//Recursively call 'PlacingRuleProcess()'
				result = PlacingRuleProcess(pThis,
					pProteinRotamerOnStartingSite,
					pCataConsArray,
					relatedProteinSiteCount,
					relatedProteinSites,
					pSmallMolRotSetForOutput,step+1);
				if(FAILED(result)){
					return result;
				}
		}
	}

	return Success;
}


int    PlacingRulePlaceSmallMol(PlacingRule* pThis,
								CataConsSitePairArray* pCataConsArray,
								DesignSite* pStartingSite,
								int relatedProteinDesignSiteCount,
								DesignSite** relatedProteinDesignSites,
								RotamerSet* pSmallMolRotSetForOutput)
{
	int    result;
	char errMsg[MAX_LENGTH_ERR_MSG+1];
	Rotamer* pCurProteinRot;
	int i;
	int progressBarWidth;
	time_t begin;
	time_t now;

	//Make Sure the Placing Rule Has Already Been Deployed.
	if(PlacingRuleGetDeployedFlag(pThis)==FALSE){
		result = AssertionError;
		sprintf(errMsg,"In PlacingRulePlaceSmallMol(), Placing Rule on %s %s %d has not been deployed",
			pThis->chainName,pThis->residueName,pThis->posInChain);
		TraceError(errMsg,result);
		return result;
	}


	//For every rotamer in the starting site's rotamer set, if this placing rule
	//applies to the rotamer type, carry out the calculation.
	begin = time(NULL);
	printf("start generating small mol rotamers\n");
	progressBarWidth = 50;
	ShowProgress(progressBarWidth,0.0);
	printf("\n");
	//fflush(stdout);
	for(i=0; i<RotamerSetGetCount(DesignSiteGetRotamers(pStartingSite)); i++){
		int elapsed;
		double percentage;

		pCurProteinRot = RotamerSetGet(DesignSiteGetRotamers(pStartingSite),i);
		RotamerRestore(pCurProteinRot,DesignSiteGetRotamers(pStartingSite));

		//The rotamer on protein design site must be the required type
		if(strcmp(pThis->rotamerType,RotamerGetType(pCurProteinRot))==0){
			//The Following Calculation is Time Consuming :
			result = PlacingRuleProcess(
				pThis,
				pCurProteinRot,
				pCataConsArray,
				relatedProteinDesignSiteCount,
				relatedProteinDesignSites,
				pSmallMolRotSetForOutput,
				0);

			if(FAILED(result)){
				sprintf(errMsg,"In PlacingRulePlaceSmallMol(), when placing the SmallMol from "
					"site %s %s %d, Rotamer type %s", pThis->chainName,
					pThis->residueName,pThis->posInChain,pThis->rotamerType);
				TraceError(errMsg,result);
				return result;
			}
		}

		RotamerExtract(pCurProteinRot);

		//Timing and show progress
		now = time(NULL);
		elapsed = (int)(now-begin);
		percentage = (double)(i+1) / (double)(RotamerSetGetCount(DesignSiteGetRotamers(pStartingSite))) * 100.0;
		//printf("\r");
		ShowProgress(progressBarWidth,percentage);
		printf("Running Time : %dh %dm %ds; %d Rotamers Found\n", 
			elapsed/3600, (elapsed%3600)/60,(elapsed%3600)%60,
			RotamerSetGetCount(pSmallMolRotSetForOutput));
		//fflush(stdout);
	}
	printf("\n");
  //fflush(stdout);
	return Success;
}


//Methods used for debugging
int    PlacingRuleShowAtom(PlacingRule* pThis,int index){
	int value;

	if(index<0 || index>=StringArrayGetCount(&pThis->atomNames))
		return IndexError;

	printf("ATOM  %2d  %5.5s  ",index,StringArrayGet(&pThis->atomNames,index));

	XYZShow(XYZArrayGet(&pThis->atomXYZs,index));
	if(pThis->deployedFlag){
		value = IntArrayGet(&pThis->atomPosOnSmallMol,index);
		printf("  %2d  ",value);
		if(value!=-1){
			printf("%5.5s  ",AtomGetName(ResidueGetAtom(pThis->pSmallMol,value)));
			XYZShow(XYZArrayGet(&pThis->smallMolAtomXYZs,value));
		}
	}
	printf("\n");
	return Success;
}
int    PlacingRuleShowAtomDistances(PlacingRule* pThis){
	int count = XYZArrayGetLength(&pThis->atomXYZs);
	int i;
	int j;
	for(i=-1;i<count;i++){
		if(i==-1){
			printf("%5.5s|","");
		}
		else{
			printf("%5.5s|",StringArrayGet(&pThis->atomNames,i));
		}
		for(j=0;j<count;j++){
			if(i==-1){
				printf("%5.5s|",StringArrayGet(&pThis->atomNames,j));
				continue;
			}

			if(i!=j){
				printf("%5.2f|",XYZDistance(
					XYZArrayGet(&pThis->atomXYZs,i),
					XYZArrayGet(&pThis->atomXYZs,j))
					);
			}
			else{
				printf("  \\  |");
			}
		}
		printf("\n");
	}
	return Success;
}
int    PlacingRuleShowParam(PlacingRule* pThis,int index){
	if(index<0 || index>=StringArrayGetCount(&pThis->paramNames))
		return IndexError;
	printf("PARAM  %2d  %8.8s  %f\n",index,
		StringArrayGet(&pThis->paramNames,index),
		DoubleArrayGet(&pThis->params,index));
	return Success;
}

int    PlacingRuleShowAction(PlacingRule* pThis,int index){
	int j;
	PlacingAction* pCurAction;

	if(index<0 || index>=pThis->actionCount)
		return IndexError;

	pCurAction = &pThis->actions[index];
	switch(pCurAction->actionType){
		case Type_PlacingAction_Load:
			for(j=0;j<IntArrayGetLength(&pCurAction->load_atoms);j++){
				index = IntArrayGet(&pCurAction->load_atoms,j);
				printf("LOAD  %s ",StringArrayGet(&pThis->atomNames,index));
				XYZShow(XYZArrayGet(&pThis->atomXYZs,index));
				printf("\n");
			}
			break;
		case Type_PlacingAction_Variate:
			printf("VARIATE  %s  %f  %f  %f\n",
				StringArrayGet(&pThis->paramNames,pCurAction->variate_param),
				pCurAction->variate_from,
				pCurAction->variate_to,pCurAction->variate_increment);
			break;
		case Type_PlacingAction_Calc:
			printf("CALC   ");
			for(j=0;j<IntArrayGetLength(&pCurAction->calc_atoms);j++){
				index = IntArrayGet(&pCurAction->calc_atoms,j);
				printf("%4.4s ",StringArrayGet(&pThis->atomNames,index));
			}
			for(j=0;j<IntArrayGetLength(&pCurAction->calc_params);j++){
				double paramValue;
				index = IntArrayGet(&pCurAction->calc_params,j);
				paramValue = DoubleArrayGet(&pThis->params,index);
				if(strcmp(StringArrayGet(&pThis->paramNames,index),"")!=0){
					printf("%8.8s",  StringArrayGet(&pThis->paramNames,index) );
				}
				else{
					if(1<=j && j<=3){
						printf("%8.2f",RadToDeg(paramValue));
					}
					else{
						printf("%8.2f",paramValue);
					}
				}
				printf("  ");
			}
			printf("\n");
			break;
		case Type_PlacingAction_Evaluate:
			printf("EVALUATE  ");
			printf("%s  %d  ",
				StringArrayGet(&pThis->paramNames,pCurAction->evaluate_param),
				pCurAction->evaluate_type);
			for(j=0;j<IntArrayGetLength(&pCurAction->evaluate_atoms);j++){
				index = IntArrayGet(&pCurAction->evaluate_atoms,j);
				printf("%4.4s ",StringArrayGet(&pThis->atomNames,index));
			}
			printf("\n");
			break;
		case Type_PlacingAction_CheckCataCons:
			printf("CHECK_CATA_CONS  ");
			printf("%s  %s  %d\t%s	%s	%d\n",
				pCurAction->checkCataCons_firstSiteChainName,
				pCurAction->checkCataCons_firstSiteResidueName,
				pCurAction->checkCataCons_firstSitePosInChain,
				pCurAction->checkCataCons_secondSiteChainName,
				pCurAction->checkCataCons_secondSiteResidueName,
				pCurAction->checkCataCons_secondSitePosInChain);
			break;
		case Type_PlacingAction_CheckMultiCons:
			printf("CHECK_MULTI_CONS_BEGIN\n");
			for(j=0;j<pCurAction->checkMultiCons_cataConsCount;j++){
				printf("%s  %s  %d\t%s	%s	%d\n",
					StringArrayGet(&pCurAction->checkMultiCons_firstSiteChainNames,j),
					StringArrayGet(&pCurAction->checkMultiCons_firstSiteResidueNames,j),
					IntArrayGet(&pCurAction->checkMultiCons_firstSitePosInChains,j),
					StringArrayGet(&pCurAction->checkMultiCons_secondSiteChainNames,j),
					StringArrayGet(&pCurAction->checkMultiCons_secondSiteResidueNames,j),
					IntArrayGet(&pCurAction->checkMultiCons_secondSitePosInChains,j));
			}
			printf("CHECK_MULTI_CONS_END\n");
			break;
		case Type_PlacingAction_CheckVDW_Backbone:
			printf("CHECK_VDW_BACKBONE   %s  %f\n",
				pCurAction->checkVDW_backbone_withHydrogen? "WITH_HYDROGEN" : "WITHOUT_HYDROGEN",
				pCurAction->checkVDW_backbone_maxAllowed);
			for(j=0;j<IntArrayGetLength(&pCurAction->checkVDW_backbone_smallmolAtomHasXyz);j++){
				Atom* pAtom = ResidueGetAtom(pThis->pSmallMol,j);
				BOOL include = IntArrayGet(&pCurAction->checkVDW_backbone_smallmolAtomHasXyz,j);
				if(include){
					printf("%s ",AtomGetName(pAtom));
				}
			}printf("\n");
			break;
		case Type_PlacingAction_CheckVDW_Internal:
			printf("CHECK_VDW_INTERNAL   %s  %f\n",
				pCurAction->checkVDW_internal_withHydrogen? "WITH_HYDROGEN" : "WITHOUT_HYDROGEN",
				pCurAction->checkVDW_internal_maxAllowed);
			for(j=0;j<IntArrayGetLength(&pCurAction->checkVDW_internal_smallmolAtomHasXyz);j++){
				Atom* pAtom = ResidueGetAtom(pThis->pSmallMol,j);
				BOOL include = IntArrayGet(&pCurAction->checkVDW_internal_smallmolAtomHasXyz,j);
				if(include){
					printf("%s ",AtomGetName(pAtom));
				}
			}printf("\n");
			break;
		case Type_PlacingAction_CheckRMSD:
			printf("CHECK_RMSD   %s  %f\n",
				pCurAction->checkRMSD_withHydrogen? "WITH_HYDROGEN" : "WITHOUT_HYDROGEN",
				pCurAction->checkRMSD_minDifference);
			break;
		default:
			printf("small molecule placing action cannot be handled\n");
	}
	return Success;
}
int    PlacingRuleShow(PlacingRule* pThis){
	int i;
	for(i=0;i<StringArrayGetCount(&pThis->atomNames);i++){
		PlacingRuleShowAtom(pThis,i);
	}
	for(i=0;i<StringArrayGetCount(&pThis->paramNames);i++){
		PlacingRuleShowParam(pThis,i);
	}
	for(i=0;i<pThis->actionCount;i++){
		PlacingRuleShowAction(pThis,i); 
	}
	return Success;
}
int    PlacingRuleTester(char* filename){
	int    result;
	PlacingRule rule;

	result = PlacingRuleCreate(&rule,filename);
	if(FAILED(result)){
		return result;
	}

	PlacingRuleDestroy(&rule);
	return Success;
}



int    ScreenSmallmolRotamersByRMSD(char* oriFileName,char* outputFileName,double rmsdThresold){
	FILE* fin = fopen(oriFileName,"r");
	FILE* fout = fopen(outputFileName,"w");
	char errMsg[MAX_LENGTH_ERR_MSG+1];

	if(fin==NULL || fout==NULL){
		sprintf(errMsg,"In ScreenSmallmolRotamersByRMSD, cannot open %s or %s",oriFileName,outputFileName);
		TraceError(errMsg,IOError);
		return IOError;
	}

	int atomCounter = 0;
  char line[MAX_LENGTH_ONE_LINE_IN_FILE];
	while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,fin)){
		char keyword[10];
		ExtractFirstStringFromSourceString(keyword,line);
		if(strcmp(keyword,"ATOM")==0){
			atomCounter++;
		}
		else if(strcmp(keyword,"ENDMDL")==0){
			break;
		}
	}
  Rotamer currentRotamer;
  RotamerCreate(&currentRotamer);
  XYZArrayResize(&currentRotamer.xyzs,atomCounter);

	fseek(fin,0,SEEK_SET);

	int totalRotamerCount = 0;	
	BOOL lastAccepted = FALSE;
  RotamerSet acceptedRotamers;
  RotamerSetCreate(&acceptedRotamers);
  StringArray buffer;
  StringArrayCreate(&buffer);
	while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,fin)){
		char keyword[10]="";
		ExtractTargetStringFromSourceString(keyword,line,0,4);
		if(strcmp(keyword,"ENDM")==0){
			lastAccepted = TRUE;
			for(int i=0;i<RotamerSetGetCount(&acceptedRotamers);i++){
				double rmsd = XYZArrayRMSD(
					&RotamerSetGet(&acceptedRotamers,i)->xyzs,
					&currentRotamer.xyzs);
				if(rmsd < rmsdThresold){
					lastAccepted = FALSE;
					break;
				}
			}
			if(lastAccepted){
				fprintf(fout,"MODEL     %d\n",RotamerSetGetCount(&acceptedRotamers));
				for(int i=0;i<StringArrayGetCount(&buffer);i++){
					fprintf(fout,"%s",StringArrayGet(&buffer,i));
				}
				fprintf(fout,"ENDMDL\n");
				RotamerSetAdd(&acceptedRotamers,&currentRotamer);
			}

			totalRotamerCount++;
			printf("\r %d / %d Rotamers Processed.        ",
				RotamerSetGetCount(&acceptedRotamers),
				totalRotamerCount);
			//fflush(stdout);
		}
		else if(strcmp(keyword,"MODE")==0){
			atomCounter = 0;
			StringArrayDestroy(&buffer);
			StringArrayCreate(&buffer);
		}
		else if(strcmp(keyword,"ATOM")==0){
			XYZ* pXYZ = XYZArrayGet(&currentRotamer.xyzs,atomCounter);
			ExtractTargetStringFromSourceString(keyword,line,31,7);
			pXYZ->X = atof(keyword);
			ExtractTargetStringFromSourceString(keyword,line,39,7);
			pXYZ->Y = atof(keyword);
			ExtractTargetStringFromSourceString(keyword,line,47,7);
			pXYZ->Z = atof(keyword);
			StringArrayAppend(&buffer,line);
			atomCounter++;
		}
		else if(strcmp(keyword,"ENER")==0 && lastAccepted==TRUE){
			fprintf(fout,"%s",line);
		}
	}


	StringArrayDestroy(&buffer);
	RotamerDestroy(&currentRotamer);
	RotamerSetDestroy(&acceptedRotamers);
	fclose(fin);
	fclose(fout);

	return Success;

}


int    AnalyzeSmallMolRotamers(char* oriFileName,Residue* pNativeSmallMol){
	int atomCounter;
	int rotamerCounter;
	char errMsg[MAX_LENGTH_ERR_MSG+1];
	char line[MAX_LENGTH_ONE_LINE_IN_FILE];
	double minInternalVDW = 1e8;
	double minBackboneVDW = 1e8;
	double internalVDWOfMinRMSD = 1e4;
	double backboneVDWOfMinRMSD = 1e4;
	double minRMSD = 1e8;
	int minRMSDindex = -1;

	XYZArray currentRotamerXYZ;
	IntArray atomPosOnNativeSmallMol;
	FILE* fin = fopen(oriFileName,"r");
	if(fin==NULL){
		sprintf(errMsg,"in file %s function %s line %d, can not open file %s",__FILE__,__FUNCTION__,__LINE__,oriFileName);
		TraceError(errMsg,IOError);
		return IOError;
	}
	IntArrayCreate(&atomPosOnNativeSmallMol,0);
	XYZArrayCreate(&currentRotamerXYZ,0);


	atomCounter = 0;
	rotamerCounter = 0;
	while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,fin)){
		char keyword[10];
		char atomName[10];
		ExtractTargetStringFromSourceString(keyword,line,0,6);
		if(strcmp(keyword,"ATOM")==0){
			int posOnNativeSmallMol;
			atomCounter++;
			ExtractTargetStringFromSourceString(atomName,line,12,4);
			ResidueFindAtom(pNativeSmallMol,atomName,&posOnNativeSmallMol);
			if(posOnNativeSmallMol==-1){
				sprintf(errMsg,"In AnalyzeSmallMolRotamers(), cannot find atom %s on smallmol",atomName);
				TraceError(errMsg,DataNotExistError);
				return DataNotExistError;
			}
			IntArrayAppend(&atomPosOnNativeSmallMol,posOnNativeSmallMol);
		}
		else if(strcmp(keyword,"ENDMDL")==0){
			break;
		}
	}

	XYZArrayResize(&currentRotamerXYZ,atomCounter);

	fseek(fin,0,SEEK_SET);

	while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,fin)){
		char keyword[10];
		strcpy(keyword,"");
		ExtractTargetStringFromSourceString(keyword,line,0,4);
		if(strcmp(keyword,"ENDM")==0){
			double rmsd = 0.0;
			int i;
			for(i=0;i<atomCounter;i++){
				Atom* pNativeAtom = 
					ResidueGetAtom(pNativeSmallMol,IntArrayGet(&atomPosOnNativeSmallMol,i));
				double dist = XYZDistance(XYZArrayGet(&currentRotamerXYZ,i),
					&pNativeAtom->xyz);
				rmsd += dist*dist;
			}
			rmsd = sqrt(rmsd/atomCounter);

			if(rmsd<minRMSD){
				minRMSD = rmsd;
				minRMSDindex = rotamerCounter-1;
			}
		}
		else if(strcmp(keyword,"MODE")==0){
			atomCounter = 0;
			rotamerCounter++;
		}
		else if(strcmp(keyword,"ATOM")==0){
			XYZ* pXYZ = XYZArrayGet(&currentRotamerXYZ,atomCounter);
			ExtractTargetStringFromSourceString(keyword,line,31,7);
			pXYZ->X = atof(keyword);
			ExtractTargetStringFromSourceString(keyword,line,39,7);
			pXYZ->Y = atof(keyword);
			ExtractTargetStringFromSourceString(keyword,line,47,7);
			pXYZ->Z = atof(keyword);
			atomCounter++;
		}
		else if(strcmp(keyword,"ENER")==0){
			char internalVDW[20];
			char backboneVDW[20];
			ExtractFirstStringFromSourceString(keyword,line);
      ExtractFirstStringFromSourceString(internalVDW,line);
			ExtractFirstStringFromSourceString(internalVDW,line);
			ExtractFirstStringFromSourceString(backboneVDW,line);
      ExtractFirstStringFromSourceString(backboneVDW,line);
			if(atof(internalVDW)<minInternalVDW){
				minInternalVDW = atof(internalVDW);
			}
			if(atof(backboneVDW)<minBackboneVDW){
				minBackboneVDW = atof(backboneVDW);
			}
			if(rotamerCounter==minRMSDindex+1){
				internalVDWOfMinRMSD = atof(internalVDW);
				backboneVDWOfMinRMSD = atof(backboneVDW);
			}

		}

		else{
			continue;
		}
	}

	printf("Total Rotamer Count    : %d\n",rotamerCounter);
	printf("Min RMSD               : %f\n",minRMSD);
	printf("Min RMSD Rotamer Index : %d\n",minRMSDindex);
	printf("InternalVDW of minRMSD : %f\n",internalVDWOfMinRMSD);
	printf("BackboneVDW of minRMSD : %f\n",backboneVDWOfMinRMSD);
	printf("Min InternalVDW        : %f\n",minInternalVDW);
	printf("Min BackboneVDW        : %f\n",minBackboneVDW);

	fclose(fin);
	XYZArrayDestroy(&currentRotamerXYZ);
	IntArrayDestroy(&atomPosOnNativeSmallMol);


	return Success;
}


int    AnalyzeSmallMolRotamersForSpecifiedAtoms(char* oriFileName,Residue* pNativeSmallMol,char* specificAtomFile){
	int i, atomCounter;
	int rotamerCounter;
	char err_msg[MAX_LENGTH_ERR_MSG+1];
	char line[MAX_LENGTH_ONE_LINE_IN_FILE];
	double minInternalVDW = 1e8;
	double minBackboneVDW = 1e8;
	double internalVDWOfMinRMSD = 1e4;
	double backboneVDWOfMinRMSD = 1e4;
	double minRMSD = 1e8;
	int minRMSDindex = -1;

	FileReader file;
	StringArray specificAtomNames;

	XYZArray currentRotamerXYZ;
	IntArray atomPosOnNativeSmallMol;
	FILE* fin = fopen(oriFileName,"r");
	if(fin==NULL){
		sprintf(err_msg,"In AnalyzeSmallMolRotamers, cannot open %s",oriFileName);
		TraceError(err_msg,IOError);
		return IOError;
	}
	IntArrayCreate(&atomPosOnNativeSmallMol,0);
	XYZArrayCreate(&currentRotamerXYZ,0);


	StringArrayCreate(&specificAtomNames);
	FileReaderCreate(&file, specificAtomFile);
	while(!FAILED(FileReaderGetNextLine(&file, line))){
		StringArray wordsInLine;
		StringArrayCreate(&wordsInLine);
		StringArraySplitString(&wordsInLine, line, ' ');
		for(i=0; i<StringArrayGetCount(&wordsInLine); i++){
			int posOnNativeSmallMol=-1;
			ResidueFindAtom(pNativeSmallMol, StringArrayGet(&wordsInLine, i), &posOnNativeSmallMol);
			if(posOnNativeSmallMol==-1){
				sprintf(err_msg,"In AnalyzeSmallMolRotamers(), cannot find atom %s on smallmol",StringArrayGet(&wordsInLine, i));
				TraceError(err_msg,DataNotExistError);
				return DataNotExistError;
			}
			StringArrayAppend(&specificAtomNames, StringArrayGet(&wordsInLine, i));
		}
		StringArrayDestroy(&wordsInLine);
	}


	atomCounter = 0;
	rotamerCounter = 0;
	while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,fin)){
		char keyword[10];
		char atomName[10];
		ExtractTargetStringFromSourceString(keyword,line,0,6);
		if(strcmp(keyword,"ATOM")==0){
			int posOnNativeSmallMol;
			
			ExtractTargetStringFromSourceString(atomName,line,12,4);
			for(i=0; i<StringArrayGetCount(&specificAtomNames); i++){
				if(strcmp(atomName, StringArrayGet(&specificAtomNames, i))==0){
					atomCounter++;
					ResidueFindAtom(pNativeSmallMol,atomName,&posOnNativeSmallMol);
					if(posOnNativeSmallMol==-1){
						sprintf(err_msg,"In AnalyzeSmallMolRotamers(), cannot find atom %s on smallmol",atomName);
						TraceError(err_msg,DataNotExistError);
						return DataNotExistError;
					}
					IntArrayAppend(&atomPosOnNativeSmallMol,posOnNativeSmallMol);
				}
			}
		}
		else if(strcmp(keyword,"ENDMDL")==0){
			break;
		}
	}

	XYZArrayResize(&currentRotamerXYZ,atomCounter);

	fseek(fin,0,SEEK_SET);

	while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,fin)){
		char keyword[10];
		strcpy(keyword,"");
		ExtractTargetStringFromSourceString(keyword,line,0,4);
		if(strcmp(keyword,"ENDM")==0){
			double rmsd = 0.0;
			int i;
			for(i=0;i<atomCounter;i++){
				Atom* pNativeAtom = 
					ResidueGetAtom(pNativeSmallMol,IntArrayGet(&atomPosOnNativeSmallMol,i));
				double dist = XYZDistance(XYZArrayGet(&currentRotamerXYZ,i),
					&pNativeAtom->xyz);
				rmsd += dist*dist;
			}
			rmsd = sqrt(rmsd/atomCounter);

			if(rmsd<minRMSD){
				minRMSD = rmsd;
				minRMSDindex = rotamerCounter-1;
			}
		}
		else if(strcmp(keyword,"MODE")==0){
			atomCounter = 0;
			rotamerCounter++;
		}
		else if(strcmp(keyword,"ATOM")==0){
			char atomName[MAX_LENGTH_ATOM_NAME+1];

			ExtractTargetStringFromSourceString(atomName,line,12,4);
			for(i=0; i<StringArrayGetCount(&specificAtomNames); i++){
				if(strcmp(atomName, StringArrayGet(&specificAtomNames, i))==0){
					XYZ* pXYZ = XYZArrayGet(&currentRotamerXYZ,atomCounter);
					ExtractTargetStringFromSourceString(keyword,line,31,7);
					pXYZ->X = atof(keyword);
					ExtractTargetStringFromSourceString(keyword,line,39,7);
					pXYZ->Y = atof(keyword);
					ExtractTargetStringFromSourceString(keyword,line,47,7);
					pXYZ->Z = atof(keyword);
					atomCounter++;
				}
			}
			
		}
		else if(strcmp(keyword,"ENER")==0){
			char internalVDW[20];
			char backboneVDW[20];
			ExtractFirstStringFromSourceString(keyword,line);
			ExtractFirstStringFromSourceString(internalVDW,line);
      ExtractFirstStringFromSourceString(internalVDW,line);
			ExtractFirstStringFromSourceString(backboneVDW,line);
      ExtractFirstStringFromSourceString(backboneVDW,line);
			if(atof(internalVDW)<minInternalVDW){
				minInternalVDW = atof(internalVDW);
			}
			if(atof(backboneVDW)<minBackboneVDW){
				minBackboneVDW = atof(backboneVDW);
			}
			if(rotamerCounter==minRMSDindex+1){
				internalVDWOfMinRMSD = atof(internalVDW);
				backboneVDWOfMinRMSD = atof(backboneVDW);
			}

		}

		else{
			continue;
		}
	}

	printf("Total Rotamer Count    : %d\n",rotamerCounter);
	printf("Min RMSD               : %f\n",minRMSD);
	printf("Min RMSD Rotamer Index : %d\n",minRMSDindex);
	printf("InternalVDW of minRMSD : %f\n",internalVDWOfMinRMSD);
	printf("BackboneVDW of minRMSD : %f\n",backboneVDWOfMinRMSD);
	printf("Min InternalVDW        : %f\n",minInternalVDW);
	printf("Min BackboneVDW        : %f\n",minBackboneVDW);

	fclose(fin);
	XYZArrayDestroy(&currentRotamerXYZ);
	IntArrayDestroy(&atomPosOnNativeSmallMol);
	FileReaderDestroy(&file);
	StringArrayDestroy(&specificAtomNames);

	return Success;
}


typedef struct  
{
	double vdwInternal;
	double vdwBackbone;
	int index;
} RotamerEnergy;


int CompareByInternalEnergy(const void *a, const void *b)
{
	double e1 = (*(RotamerEnergy*)a).vdwInternal;
	double e2 = (*(RotamerEnergy*)b).vdwInternal;

	return e1 < e2 ? -1 : 1;
}

int CompareByBackboneEnergy(const void *a, const void *b)
{
	double e1 = (*(RotamerEnergy*)a).vdwBackbone;
	double e2 = (*(RotamerEnergy*)b).vdwBackbone;

	return e1 < e2 ? -1 : 1;
}

int SmallMolRotamersGetBothHighRankOfBackboneVdwAndInternalVdw(char *oldRotamersFile,char* newRotamersFile,double percent){
  printf("select smallmol rotamers that have good internal and backbone VDW energy\n");
  FILE* fin=fopen(oldRotamersFile,"r");
  int smallmolRotamerNum=0;
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
	while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,fin)){
		char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1];
		ExtractFirstStringFromSourceString(keyword,line);
		if(strcmp(keyword, "MODEL") == 0) smallmolRotamerNum++;
	}
  printf("total smallmol rotamer count: %d\n",smallmolRotamerNum);
	RotamerEnergy* pRotamerEnergy = (RotamerEnergy *)malloc(sizeof(RotamerEnergy)*smallmolRotamerNum);
  BOOL* flagRotamerWithinRank=(BOOL*)malloc(sizeof(BOOL)*smallmolRotamerNum);
  for(int i=0;i<smallmolRotamerNum;i++) flagRotamerWithinRank[i]=FALSE;

  fseek(fin,0,SEEK_SET);
	int smallmolRotamerIndex = 0;
	while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,fin)){
		double vdwInternal = 0.0, vdwBackbone = 0.0;
		char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1];
		ExtractFirstStringFromSourceString(keyword, line);
		if(strcmp(keyword, "ENERGY") != 0) continue;
		ExtractFirstStringFromSourceString(keyword, line);
    ExtractFirstStringFromSourceString(keyword, line);
		vdwInternal = atof(keyword);
		ExtractFirstStringFromSourceString(keyword, line);
    ExtractFirstStringFromSourceString(keyword, line);
		vdwBackbone = atof(keyword);
		pRotamerEnergy[smallmolRotamerIndex].vdwInternal = vdwInternal;
		pRotamerEnergy[smallmolRotamerIndex].vdwBackbone = vdwBackbone;
		pRotamerEnergy[smallmolRotamerIndex].index = smallmolRotamerIndex;
		smallmolRotamerIndex++;
	}
	printf("rank smallmol rotamers by vdwInternal and vdwBackbone energy\n");

  int rankThreshold=(int)(smallmolRotamerNum*percent);

	// rank by vdwInternal;
  IntArray highRankIndexByInternal;
	IntArrayCreate(&highRankIndexByInternal, rankThreshold);
	qsort(pRotamerEnergy, smallmolRotamerNum, sizeof(RotamerEnergy), CompareByInternalEnergy);
	smallmolRotamerIndex = 0;
	while(smallmolRotamerIndex < rankThreshold){
		IntArraySet(&highRankIndexByInternal, smallmolRotamerIndex, pRotamerEnergy[smallmolRotamerIndex].index);
		smallmolRotamerIndex++;
	}

	// rank by vdwBackbone;
  IntArray highRankIndexByBackbone;
	IntArrayCreate(&highRankIndexByBackbone, rankThreshold);
	qsort(pRotamerEnergy, smallmolRotamerNum, sizeof(RotamerEnergy), CompareByBackboneEnergy);
	smallmolRotamerIndex = 0;
	while(smallmolRotamerIndex < rankThreshold){
		IntArraySet(&highRankIndexByBackbone, smallmolRotamerIndex, pRotamerEnergy[smallmolRotamerIndex].index);
		smallmolRotamerIndex++;
	}

	// find out the rotamers that are ranked in top 'highRank' by vdwInternal as well as vdwBackbone;
	printf("rotamer rank in top %.0f%% (%d rotamers) by both internal and backbone energy\n", percent*100,rankThreshold);
	for(int i = 0; i < IntArrayGetLength(&highRankIndexByInternal); i++){
		int indexI = IntArrayGet(&highRankIndexByInternal, i);
		for(int j = 0; j < IntArrayGetLength(&highRankIndexByBackbone); j++){
			int indexJ = IntArrayGet(&highRankIndexByBackbone, j);
			if(indexJ == indexI){
				printf("RotamerIndex: %6d, InternalRank: %6d, BackboneRank: %6d\n", indexJ, i, j);
        flagRotamerWithinRank[indexI]=TRUE;
			}
		}
	}

  //write top rotamers to a new file
  FILE* fout=fopen(newRotamersFile,"w");
  if(fout==NULL){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    sprintf(errMsg,"in file %s function %s line %d, can not open file %s for writing",__FILE__,__FUNCTION__,__LINE__,newRotamersFile);
    TraceError(errMsg,IOError);
    return IOError;
  }

  fseek(fin,0,SEEK_SET);
  smallmolRotamerIndex=0;
  while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,fin)){
    char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    ExtractTargetStringFromSourceString(keyword,line,0,4);
    if(strcmp(keyword, "MODE") == 0) smallmolRotamerIndex++;
    if(flagRotamerWithinRank[smallmolRotamerIndex-1]==TRUE){
      fprintf(fout,"%s",line);
    }
  }
  fclose(fout);

	IntArrayDestroy(&highRankIndexByBackbone);
	IntArrayDestroy(&highRankIndexByInternal);
	free(pRotamerEnergy);
  free(flagRotamerWithinRank);

	return Success;
}

