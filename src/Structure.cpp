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

#include "Structure.h"
#include "EnergyFunction.h"
#include <string.h>

//#define DEBUGGING_STRUCTURE


int StructureCreate(Structure* pThis){
  strcpy(pThis->name, "");
  pThis->chainNum = 0;
  pThis->chains = NULL;
  pThis->designSiteCount = 0;
  pThis->designSites = NULL;
  return Success;
}
int StructureDestroy(Structure* pThis){
  //first free the memory for design sites
  for(int i=0; i<pThis->designSiteCount; i++){
    DesignSiteDestroy(&pThis->designSites[i]);
  }
  pThis->designSites=NULL;
  pThis->designSiteCount=0;
  for(int i=0;i<pThis->chainNum;i++){
    ChainDestroy(&pThis->chains[i]);
  }
  free(pThis->chains);
  strcpy(pThis->name, "");
  pThis->chainNum = 0;
  pThis->chains = NULL;
  return Success;
}

char* StructureGetName(Structure* pThis){
  return pThis->name;
}

int StructureSetName(Structure* pThis, char* newName){
  if(strlen(newName)>MAX_LENGTH_STRUCTURE_NAME){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->name, newName);
  return Success;
}

int StructureGetChainCount(Structure* pThis){
  return pThis->chainNum;
}

Chain* StructureGetChain(Structure* pThis, int index){
  if(index<0 || index>=StructureGetChainCount(pThis)){
    return NULL;
  }
  return &pThis->chains[index];
}

Chain* StructureFindChainByName(Structure* pThis, char* chainName){
  int index = -1;
  int result = StructureFindChainIndex(pThis, chainName, &index);
  if(FAILED(result)){
    return NULL;
  }
  else{
    return StructureGetChain(pThis, index);
  }
}

int StructureFindChainIndex(Structure* pThis, char* chainName, int* index){
  int i;
  for(i=0;i<pThis->chainNum;i++){
    if(strcmp(ChainGetName(&pThis->chains[i]), chainName)==0){
      *index = i;
      return Success;
    }
  }
  return DataNotExistError;
}


int StructureAddChain(Structure* pThis, Chain* newChain){
  int index = -1;
  int result = StructureFindChainIndex(pThis, ChainGetName(newChain), &index);
  if(FAILED(result)){
    (pThis->chainNum)++;
    pThis->chains = (Chain*)realloc(pThis->chains, sizeof(Chain)*pThis->chainNum);
    ChainCreate(&pThis->chains[pThis->chainNum-1]);
    return ChainCopy(&pThis->chains[pThis->chainNum-1], newChain);
  }
  else{
    return ChainCopy(&pThis->chains[index], newChain);
  }
}

int StructureDeleteChain(Structure* pThis, char* chainName){
  int index;
  int result = StructureFindChainIndex(pThis, chainName, &index);
  if(FAILED(result))
    return result;
  for(int i=index;i<pThis->chainNum-1;i++){
    ChainCopy(&pThis->chains[i], &pThis->chains[i+1]);
  }
  ChainDestroy(&pThis->chains[pThis->chainNum-1]);
  (pThis->chainNum)--;
  return Success;
}

int StructureShowInPDBFormat(Structure* pThis, BOOL showHydrogen, FILE* pFile){
  int atomIndex=1;
  for(int i=0;i<StructureGetChainCount(pThis);i++){
    Chain* pChain = StructureGetChain(pThis, i);
    for(int j = 0; j < ChainGetResidueCount(pChain); j++){
      Residue *pResi = ChainGetResidue(pChain,j);
      ResidueShowInPDBFormat(pResi, "ATOM", ResidueGetChainName(pResi), atomIndex, ResidueGetPosInChain(pResi), showHydrogen, pFile);
      atomIndex += ResidueGetAtomCount(pResi);
    }
  }
  return Success;
}

int StructureConfig(Structure* pStructure, char* pdbFile, AtomParamsSet* pAtomParams, ResiTopoSet* pTopos){
  char initChainID[MAX_LENGTH_CHAIN_NAME+1];
  char initResPos[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  strcpy(initChainID, "UNK");
  initChainID[strlen(initChainID)]='\0';
  strcpy(initResPos, "UNKNOWN");
  initResPos[strlen(initResPos)]='\0';
  BOOL firstResidueInChain = TRUE;

  FileReader file;
  if(FAILED(FileReaderCreate(&file, pdbFile))){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    sprintf(errMsg, "in file %s function %s() line %d, when opening:\n%s",__FILE__, __FUNCTION__, __LINE__, pdbFile);
    TraceError(errMsg, IOError);
    return IOError;
  }
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(!FAILED(FileReaderGetNextLine(&file, line))){
    char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1] = "";
    char strAtomName[MAX_LENGTH_ATOM_NAME+1] = "";
    char strResName[MAX_LENGTH_RESIDUE_NAME+1] = "";
    char strChainID[MAX_LENGTH_CHAIN_NAME+1] = "";
    char strResPos[MAX_LENGTH_ONE_LINE_IN_FILE+1] = "";

    ExtractTargetStringFromSourceString(keyword, line, 0, 4);
    if(strcmp(keyword, "ATOM") != 0) continue;
    ExtractTargetStringFromSourceString(strAtomName, line, 12, 4);
    ExtractTargetStringFromSourceString(strResName, line, 17, 4);
    ExtractTargetStringFromSourceString(strChainID, line, 21, 1);
    ExtractTargetStringFromSourceString(strResPos, line, 22, 5);
    if(strcmp(strChainID, "") == 0 || strcmp(strChainID, " ") == 0){
      //int result = Warning;
      //char errMsg[MAX_LENGTH_ERR_MSG+1];
      //sprintf(errMsg, "in file %s function %s() line %d, no chain ID identified from line %s, we set it as 'A' in default",  __FILE__, __FUNCTION__, __LINE__, line);
      //TraceError(errMsg, Warning);
      strcpy(strChainID, "A");
    }
    // new chain
    if(strcmp(initChainID, strChainID) != 0){
      //patch C-ter of the old chain if necessary
      Chain *pChain = StructureGetChain(pStructure,StructureGetChainCount(pStructure)-1);
      if(pChain != NULL){
        if(ChainGetType(pChain) == Type_Chain_Protein){
          Residue *pFirsResidueInChain = ChainGetResidue(pChain,0);
          Residue *pLastResidueInChain = ChainGetResidue(pChain,ChainGetResidueCount(pChain)-1);
          if(ResidueGetAtomByName(pFirsResidueInChain, "HT1") != NULL || ResidueGetAtomByName(pFirsResidueInChain, "HN1") !=NULL){
            // if the first residue is patched with NTER, the last residue will be patched CTER to keep charge balance
            ResiduePatchCTER(pLastResidueInChain, "CTER", pAtomParams, pTopos);
            pLastResidueInChain->resiTerm = Type_ResidueIsCter;
          }
        }
        ChainCalcAllAtomXYZ(pChain, pTopos);
      }

      //deal with new chains
      strcpy(initChainID, strChainID);
      Chain newChain;
      Type_Chain chainType;
      ChainCreate(&newChain);
      //set the chain type
      chainType = ChainTypeIdentifiedFromResidueName(strResName);
      ChainSetType(&newChain, chainType);
      ChainSetName(&newChain, strChainID);
      StructureAddChain(pStructure, &newChain);
      ChainDestroy(&newChain);
      FileReaderSetCurrentPos(&file, FileReaderGetCurrentPos(&file)-1);
      firstResidueInChain = TRUE;
    }
    else{
      // new residue
      if(strcmp(initResPos, strResPos) != 0){
        Residue newResi;
        ResidueCreate(&newResi);
        if(strcmp(strResName, "HIS")==0) strcpy(strResName, "HSD");
        ResidueSetName(&newResi, strResName);
        ResidueSetChainName(&newResi,strChainID);
        ResidueSetPosInChain(&newResi, atoi(strResPos));
        ResidueAddAtomsFromAtomParams(&newResi, pAtomParams);
        ResidueAddBondsFromResiTopos(&newResi, pTopos);
        if(firstResidueInChain && StructureGetChain(pStructure, StructureGetChainCount(pStructure)-1)->type == Type_Chain_Protein){
          ResiduePatchNTERorCTER(&newResi, "NTER", pAtomParams, pTopos);
          newResi.resiTerm = Type_ResidueIsNter;
          firstResidueInChain = FALSE;
        }
        FileReaderSetCurrentPos(&file, FileReaderGetCurrentPos(&file)-1);
        // Residue read XYZ from FileReader
        ResidueReadXYZFromPDB(&newResi, &file);
        ResidueCheckAtomCoordinateValidity(&newResi);
        if(StructureGetChain(pStructure, StructureGetChainCount(pStructure)-1) == NULL){
          char errMsg[MAX_LENGTH_ERR_MSG+1];
          sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
          TraceError(errMsg, ValueError);
          return ValueError;
        }
        else{
          ChainAppendResidue(StructureGetChain(pStructure, StructureGetChainCount(pStructure)-1), &newResi);
        }
        ResidueDestroy(&newResi);
        strcpy(initResPos, strResPos);
      }
    }
  }

  //patch C-ter of the last chain if necessary
  Chain *pChain = StructureGetChain(pStructure,StructureGetChainCount(pStructure)-1);
  if(pChain != NULL){
    if(ChainGetType(pChain) == Type_Chain_Protein){
      Residue *pFirsResidueInChain = ChainGetResidue(pChain,0);
      Residue *pLastResidueInChain = ChainGetResidue(pChain,ChainGetResidueCount(pChain)-1);
      if(ResidueGetAtomByName(pFirsResidueInChain, "HT1") != NULL || ResidueGetAtomByName(pFirsResidueInChain, "HN1") !=NULL){
        // if the first residue is patched with NTER, the last residue will be patched CTER to keep charge balance
        ResiduePatchCTER(pLastResidueInChain, "CTER", pAtomParams, pTopos);
        pLastResidueInChain->resiTerm = Type_ResidueIsCter;
      }
    }
    ChainCalcAllAtomXYZ(pChain, pTopos);
  }

  FileReaderDestroy(&file);
  return Success;
}


int StructureConfigLigand(Structure* pStructure, char* mol2file, AtomParamsSet* pAtomParams, ResiTopoSet* pTopos){
  FILE* fpmol2=fopen(mol2file,"r");

  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  BOOL readingAtom = FALSE;
  BOOL readingBond = FALSE;

  char resiName[MAX_LENGTH_RESIDUE_NAME+1]="#";

  //add a new chain
  Chain newChain;
  Type_Chain chainType;
  ChainCreate(&newChain);
  chainType = Type_Chain_SmallMol;
  ChainSetType(&newChain, chainType);

  Residue newResi;
  ResidueCreate(&newResi);
  while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,fpmol2)!=NULL){
    strcpy(keyword,"");
    sscanf(line,"%s",keyword);
    if(strcmp(keyword,"@<TRIPOS>ATOM")==0){
      readingAtom = TRUE;
      readingBond = FALSE;
      continue;
    }
    else if(strcmp(keyword,"@<TRIPOS>BOND")==0){
      readingAtom = FALSE;
      readingBond = TRUE;
      continue;
    }
    else if(keyword[0] == '@'){
      readingAtom = readingBond = FALSE;
      continue;
    }

    if(readingAtom){
      char strAtomName[MAX_LENGTH_ATOM_NAME+1] = "";
      char strResName[MAX_LENGTH_RESIDUE_NAME+1] = "";
      char strChainID[MAX_LENGTH_CHAIN_NAME+1] = "";
      char strResPos[MAX_LENGTH_ONE_LINE_IN_FILE+1] = "";
      double charge=0.0;
      double x=0,y=0,z=0;
      char type[MAX_LENGTH_ATOM_TYPE+1]="";

      int atomId;
      sscanf(line,"%d %s %lf %lf %lf %s %s %s %lf",&atomId,strAtomName,&x,&y,&z,type,strChainID,strResName,&charge);
      int posInChain=atoi(strChainID);
      if(resiName[0]=='#'){
        strcpy(resiName,strResName);
        if(strlen(resiName)>3) resiName[3]='\0';
        //check if the ligand residue name conflicts with amino acid residue name
        if(LigandResidueNameConflictWithAminoAcid(resiName)){
          strcpy(resiName,"LIG");
        }
      }
      //deploy the residue first if it's not ready
      if(ResidueGetAtomCount(&newResi)==0){
        ResidueSetName(&newResi, resiName);
        ResidueSetPosInChain(&newResi, posInChain);
        ResidueAddAtomsFromAtomParams(&newResi, pAtomParams);
        ResidueAddBondsFromResiTopos(&newResi, pTopos);
        //set chain name
        ChainSetName(&newChain, strChainID);
      }
      Atom* pAtom=ResidueGetAtomByName(&newResi,strAtomName);
      pAtom->xyz.X=x;
      pAtom->xyz.Y=y;
      pAtom->xyz.Z=z;
      pAtom->isXyzValid=TRUE;
    }
  }
  ChainAppendResidue(&newChain,&newResi);
  StructureAddChain(pStructure, &newChain);
  ChainDestroy(&newChain);
  ResidueDestroy(&newResi);
  fclose(fpmol2);

  return Success;
}



int ChainComputeResiduePosition(Structure *pStructure, int chainIndex){
  Chain *pChainI = StructureGetChain(pStructure, chainIndex);
  for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
    Residue *pResiIR = ChainGetResidue(pChainI,ir);
    pResiIR->nCbIn8A = 0;
  }
  for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
    Residue *pResiIR = ChainGetResidue(pChainI,ir);
    Atom *pAtomCAorCB1 = ResidueGetAtomByName(pResiIR, "CB");
    if(pAtomCAorCB1 == NULL) pAtomCAorCB1 = ResidueGetAtomByName(pResiIR, "CA");
    for(int is = ir+1; is < ChainGetResidueCount(pChainI); is++){
      Residue *pResiIS = ChainGetResidue(pChainI,is);
      Atom *pAtomCAorCB2 = ResidueGetAtomByName(pResiIS, "CB");
      if(pAtomCAorCB2 == NULL) pAtomCAorCB2 = ResidueGetAtomByName(pResiIS, "CA");
      if(XYZDistance(&pAtomCAorCB1->xyz, &pAtomCAorCB2->xyz) < 8.0){
        pResiIR->nCbIn8A++;
        pResiIS->nCbIn8A++;
      }
    }
    //printf("Residue: %s %d %s, NumCBwithin8AtoCurResi: %d\n",ResidueGetChainName(pResiIR), ResidueGetPosInChain(pResiIR), ResidueGetName(pResiIR), pResiIR->numCBwithin8AtoCurResidue);
  }

  return Success;
}

int StructureComputeResiduePosition(Structure *pStructure){
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChainI =StructureGetChain(pStructure,i);
    for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
      Residue *pResiIR = ChainGetResidue(pChainI,ir);
      pResiIR->nCbIn8A = 0;
    }
  }
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChainI =StructureGetChain(pStructure,i);
    if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
    for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
      Residue *pResiIR = ChainGetResidue(pChainI,ir);
      Atom *pAtomCAorCB1 = ResidueGetAtomByName(pResiIR, "CB");
      if(pAtomCAorCB1==NULL) pAtomCAorCB1 = ResidueGetAtomByName(pResiIR, "CA");
      for(int is = ir+1; is < ChainGetResidueCount(pChainI); is++){
        Residue *pResiIS = ChainGetResidue(pChainI,is);
        Atom *pAtomCAorCB2 = ResidueGetAtomByName(pResiIS, "CB");
        if(pAtomCAorCB2==NULL) pAtomCAorCB2 = ResidueGetAtomByName(pResiIS, "CA");
        if(XYZDistance(&pAtomCAorCB1->xyz, &pAtomCAorCB2->xyz) < 8.0){
          pResiIR->nCbIn8A++;
          pResiIS->nCbIn8A++;
        }
      }
      for(int k = i+1; k < StructureGetChainCount(pStructure); k++){
        Chain *pChainK = StructureGetChain(pStructure,k);
        if(ChainGetType(pChainK) != Type_Chain_Protein) continue;
        for(int ks = 0; ks < ChainGetResidueCount(pChainK); ks++){
          Residue *pResiKS = ChainGetResidue(pChainK,ks);
          Atom *pAtomCAorCB2 = pAtomCAorCB2 = ResidueGetAtomByName(pResiKS, "CB");
          if(pAtomCAorCB2==NULL) pAtomCAorCB2 = ResidueGetAtomByName(pResiKS, "CA");
          if(XYZDistance(&pAtomCAorCB1->xyz, &pAtomCAorCB2->xyz) < 8.0){
            pResiIR->nCbIn8A++;
            pResiKS->nCbIn8A++;
          }
        }
      }
      //printf("Residue: %s %d %s, NumCBwithin8AtoCurResi: %d\n",ResidueGetChainName(pResiIR), ResidueGetPosInChain(pResiIR), ResidueGetName(pResiIR), pResiIR->numCBwithin8AtoCurResidue);
    }
  }

  return Success;
}

int StructureGetAminoAcidComposition(Structure* pStructure, int *aas){
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain* pChain = StructureGetChain(pStructure, i);
    for(int j = 0; j < ChainGetResidueCount(pChain); j++){
      Residue* pResi = ChainGetResidue(pChain, j);
      if(strcmp(ResidueGetName(pResi), "ALA") == 0) aas[0]++;
      else if(strcmp(ResidueGetName(pResi), "CYS") == 0) aas[1]++;
      else if(strcmp(ResidueGetName(pResi), "ASP") == 0) aas[2]++;
      else if(strcmp(ResidueGetName(pResi), "GLU") == 0) aas[3]++;
      else if(strcmp(ResidueGetName(pResi), "PHE") == 0) aas[4]++;
      else if(strcmp(ResidueGetName(pResi), "GLY") == 0) aas[5]++;
      else if(strcmp(ResidueGetName(pResi), "HSD") == 0) aas[6]++;
      else if(strcmp(ResidueGetName(pResi), "HSE") == 0) aas[6]++;
      else if(strcmp(ResidueGetName(pResi), "HIS") == 0) aas[6]++;
      else if(strcmp(ResidueGetName(pResi), "ILE") == 0) aas[7]++;
      else if(strcmp(ResidueGetName(pResi), "LYS") == 0) aas[8]++;
      else if(strcmp(ResidueGetName(pResi), "LEU") == 0) aas[9]++;
      else if(strcmp(ResidueGetName(pResi), "MET") == 0) aas[10]++;
      else if(strcmp(ResidueGetName(pResi), "ASN") == 0) aas[11]++;
      else if(strcmp(ResidueGetName(pResi), "PRO") == 0) aas[12]++;
      else if(strcmp(ResidueGetName(pResi), "GLN") == 0) aas[13]++;
      else if(strcmp(ResidueGetName(pResi), "ARG") == 0) aas[14]++;
      else if(strcmp(ResidueGetName(pResi), "SER") == 0) aas[15]++;
      else if(strcmp(ResidueGetName(pResi), "THR") == 0) aas[16]++;
      else if(strcmp(ResidueGetName(pResi), "VAL") == 0) aas[17]++;
      else if(strcmp(ResidueGetName(pResi), "TRP") == 0) aas[18]++;
      else if(strcmp(ResidueGetName(pResi), "TYR") == 0) aas[19]++;
    }
  }
  printf("\nAmino acid composition of structures:\n");
  printf("ALA =            %d\n", aas[0]);
  printf("CYS =            %d\n", aas[1]);
  printf("ASP =            %d\n", aas[2]);
  printf("GLU =            %d\n", aas[3]);
  printf("PHE =            %d\n", aas[4]);
  printf("GLY =            %d\n", aas[5]);
  printf("HIS =            %d\n", aas[6]);
  printf("ILE =            %d\n", aas[7]);
  printf("LYS =            %d\n", aas[8]);
  printf("LEU =            %d\n", aas[9]);
  printf("MET =            %d\n", aas[10]);
  printf("ASN =            %d\n", aas[11]);
  printf("PRO =            %d\n", aas[12]);
  printf("GLN =            %d\n", aas[13]);
  printf("ARG =            %d\n", aas[14]);
  printf("SER =            %d\n", aas[15]);
  printf("THR =            %d\n", aas[16]);
  printf("VAL =            %d\n", aas[17]);
  printf("TRP =            %d\n", aas[18]);
  printf("TYR =            %d\n", aas[19]);
  return Success;
}


////////////////////////////////////////////////////////////////////////////////////////////
// functional methods to check bugs in EvoEF
///////////////////////////////////////////////////////////////////////////////////////////

int StructureShowAtomParameter(Structure* pStructure){
  for(int i=0; i<StructureGetChainCount(pStructure); i++){
    ChainShowAtomParameter(StructureGetChain(pStructure,i));
  }
  return Success;
}

int StructureShowBondInformation(Structure* pStructure){
  for(int i=0; i<StructureGetChainCount(pStructure); i++){
    ChainShowBondInformation(StructureGetChain(pStructure,i));
  }
  return Success;
}


int StructureCheckIntraBondType(Structure *pStructure){
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChain=StructureGetChain(pStructure,i);
    for(int j = 0; j < ChainGetResidueCount(pChain); j++){
      Residue *pResidue = ChainGetResidue(pChain,j);
      for(int k = 0; k < ResidueGetAtomCount(pResidue); k++){
        Atom *pAtom = ResidueGetAtom(pResidue,k);
        for(int m = k+1; m < ResidueGetAtomCount(pResidue); m++){
          Atom *pAtom2 = ResidueGetAtom(pResidue,m);
          //int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtom), AtomGetName(pAtom2), pResidue);
          int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtom), AtomGetName(pAtom2), ResidueGetBonds(pResidue));
          //printf("Residue: %s %d %s, Atom1: %s, Atom2: %s, BondType: %d\n",ResidueGetName(pResidue), ResidueGetPosInChain(pResidue), ResidueGetChainName(pResidue),AtomGetName(pAtom),AtomGetName(pAtom2),bondType);
        }
      }
    }
  }
  return Success;
}

int StructureCheckNeighbouringBondType(Structure *pStructure){
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChain=StructureGetChain(pStructure,i);
    for(int j = 0; j < ChainGetResidueCount(pChain)-1; j++){
      Residue *pResi1 = ChainGetResidue(pChain,j);
      Residue *pResi2 = ChainGetResidue(pChain,j+1);
      for(int p = 0; p < pResi1->atoms.atomNum; p++){
        Atom *pAtom1 = ResidueGetAtom(pResi1,p);
        for(int q = 0; q < pResi2->atoms.atomNum; q++){
          Atom *pAtom2 = ResidueGetAtom(pResi2,q);
          int bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(pAtom1->name, pAtom2->name, pResi2->name);
          printf("Residue: %s %d %s, Atom1: %s, Residue: %s %d %s, Atom2: %s, BondType: %d\n",ResidueGetName(pResi1), ResidueGetPosInChain(pResi1), ResidueGetChainName(pResi1), AtomGetName(pAtom1), ResidueGetName(pResi2), ResidueGetPosInChain(pResi2), ResidueGetChainName(pResi2), AtomGetName(pAtom2),bondType);
        }
      }
    }
  }

  return Success;
}

///////////////////////////////////////////////
//functions for dealing with design sites
//////////////////////////////////////////////
int StructureGetDesignSiteCount(Structure* pThis){
  return pThis->designSiteCount;
}

DesignSite* StructureGetDesignSite(Structure* pThis, int index){
  if(index<0 || index>=pThis->designSiteCount) return NULL;
  else return &pThis->designSites[index];
}

DesignSite* StructureFindDesignSite(Structure* pThis, int chainIndex, int resiIndex){
  if(chainIndex<0 || chainIndex>=StructureGetChainCount(pThis)) return NULL;
  Chain* pChain = StructureGetChain(pThis, chainIndex);
  if(resiIndex<0 || resiIndex>=ChainGetResidueCount(pChain)) return NULL;
  for(int i=0; i < StructureGetDesignSiteCount(pThis); i++){
    DesignSite* pDesignSite = StructureGetDesignSite(pThis, i);
    if(pDesignSite->chainIndex == chainIndex && pDesignSite->resiIndex == resiIndex){
      return pDesignSite;
    }
  }
  return NULL;
}

DesignSite * StructureFindDesignSiteByChainName(Structure *pStructure, char *chainName, int posInChain){
  for(int i=0; i<StructureGetDesignSiteCount(pStructure); i++){
    DesignSite *pDesignSite = StructureGetDesignSite(pStructure, i);
    if(strcmp(chainName, DesignSiteGetChainName(pDesignSite)) == 0 &&
      DesignSiteGetPosInChain(pDesignSite) == posInChain){
        return pDesignSite;
    }
  }
  return NULL;
}

int DesignSiteIndexFind(Structure *pStructure, char *chainName, int posInChain){
  for(int i=0; i<StructureGetDesignSiteCount(pStructure); i++){
    DesignSite *pDesignSite = StructureGetDesignSite(pStructure, i);
    if(strcmp(DesignSiteGetChainName(pDesignSite), chainName) == 0 &&
      DesignSiteGetPosInChain(pDesignSite) == posInChain){
        return i;
    }
  }
  return -1;
}

int DesignSiteIndexGet(Structure *pStructure,int chainIndex, int resiIndex){
  for(int i=0; i<StructureGetDesignSiteCount(pStructure); i++){
    DesignSite *pDesignSite = StructureGetDesignSite(pStructure, i);
    if(pDesignSite->chainIndex==chainIndex && pDesignSite->resiIndex == resiIndex){
        return i;
    }
  }
  return -1;
}


int StructureShowDesignSites(Structure* pThis, FILE* pFile){
  if(pFile == NULL) pFile = stdout;
  double conformationSpace = 0.0;
  int totalRotamerCount = 0;
  for(int i=0; i<pThis->designSiteCount; i++){
    DesignSite* pDesignSite = StructureGetDesignSite(pThis, i);
    int rotamerCountOfResidue = RotamerSetGetCount(DesignSiteGetRotamers(pDesignSite));
    fprintf(pFile,"design site %3d : %3s %s %4d, %6d rotamers.  ",
      i,ResidueGetName(pDesignSite->pResidue),ResidueGetChainName(pDesignSite->pResidue),ResidueGetPosInChain(pDesignSite->pResidue),rotamerCountOfResidue);
    totalRotamerCount += rotamerCountOfResidue;
    if(rotamerCountOfResidue>0){
      conformationSpace += log((double)rotamerCountOfResidue)/log(10.0);
    }
    switch(pDesignSite->pResidue->designSiteType){
      case Type_ResidueDesignType_Catalytic:
        fprintf(pFile,"catalytic\n"); break;
      case Type_ResidueDesignType_Fixed:
        fprintf(pFile,"fixed\n"); break;
      case Type_ResidueDesignType_Mutated:
        fprintf(pFile,"mutated\n"); break;
      case Type_ResidueDesignType_SmallMol:
        fprintf(pFile,"smallmol\n"); break;
      case Type_ResidueDesignType_Rotameric:
        fprintf(pFile,"rotameric\n"); break;
      default:
        break;
    }
#ifdef DEBUGGING_STRUCTURE
    switch(ChainGetType(StructureGetChain(pThis, pDesignSite->chainIndex))){
      case Type_Chain_Protein:
        fprintf(pFile,"protein\n"); break;
      case Type_Chain_SmallMol:
        fprintf(pFile,"small molecule\n"); break;
      default:
        break;
    }
#endif
  }
  fprintf(pFile, "total rotamer count: %d, conformation space: 1e^%.1f\n", totalRotamerCount, conformationSpace);
  return Success;
}

int ProteinSiteAddDesignSite(Structure* pThis, int chainIndex, int resiIndex){
  //this function add a new empty design site
  int result = Success;
  DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if(pCurrentDesignSite != NULL){
    return Success;
  }
  (pThis->designSiteCount)++;
  pThis->designSites = (DesignSite*)realloc(pThis->designSites, sizeof(DesignSite)*pThis->designSiteCount);
  DesignSiteCreate(&pThis->designSites[pThis->designSiteCount-1]);
  pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
  pCurrentDesignSite->pResidue = ChainGetResidue(StructureGetChain(pThis, chainIndex), resiIndex);
  pCurrentDesignSite->chainIndex = chainIndex;
  pCurrentDesignSite->resiIndex = resiIndex;
  return result;
}

int StructureRemoveAllDesignSites(Structure* pThis){
  for(int i = 0; i < pThis->designSiteCount; i++){
    DesignSite* pDesignSite = &pThis->designSites[i];
    DesignSiteDestroy(pDesignSite);
    pDesignSite = NULL;
  }
  free(pThis->designSites);
  pThis->designSites = NULL;
  pThis->designSiteCount=0;
  return Success;
}

int ProteinSiteRemoveDesignSite(Structure* pThis, int chainIndex, int resiIndex){
  DesignSite* pDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if(pDesignSite != NULL){
    DesignSiteDestroy(pDesignSite);
    pDesignSite=NULL;
  }
  (pThis->designSiteCount)--;
  return Success;
}



int StructureFindSmallMol(Structure* pThis, Residue** ppSmallMol){
  int result = DataNotExistError;
  for(int i=0;i<pThis->chainNum;i++){
    if(pThis->chains[i].type == Type_Chain_SmallMol){
      *ppSmallMol = ChainGetResidue(&pThis->chains[i], 0);
      if(*ppSmallMol != NULL){
        result = Success;
        break;
      } 
    }
  }
  return result;
}


int StructureCalcProteinResidueSidechainTorsion(Structure* pThis, ResiTopoSet* pResiTopos){
  for(int i=0;i<StructureGetChainCount(pThis);i++){
    Chain* pChain=StructureGetChain(pThis,i);
    if(ChainGetType(pChain)==Type_Chain_Protein){
      for(int j=0;j<ChainGetResidueCount(pChain);j++){
        Residue* pResidue=ChainGetResidue(pChain,j);
        if(strcmp(ResidueGetName(pResidue),"ALA")==0 || strcmp(ResidueGetName(pResidue),"GLY")==0) continue;
        ResidueCalcSidechainTorsion(pResidue,pResiTopos);
      }
    }
  }
  return Success;
}

int StructureCopy(Structure* pThis, Structure* pOther){
  for(int i=0;i<StructureGetChainCount(pOther);i++){
    //Chain tempChain;
    //ChainCreate(&tempChain);
    //ChainCopy(&tempChain,&pOther->chains[i]);
    StructureAddChain(pThis,&pOther->chains[i]);
    //ChainDestroy(&tempChain);
  }
  pThis->chainNum=pOther->chainNum;
  pThis->designSiteCount=pOther->designSiteCount;
  strcpy(pThis->name,pOther->name);
  for(int i=0;i<pOther->designSiteCount;i++){
    DesignSiteCopy(&pThis->designSites[i],&pOther->designSites[i]);
  }
  return Success;
}