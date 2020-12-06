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

#pragma warning(disable:4996)
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "Getopt.h"
#include "ProgramFunction.h"
#include "RotamerBuilder.h"
#include "Structure.h"
#include "EnergyMatrix.h"
#include "Evolution.h"
#include "EnergyOptimization.h"
#include "WeightOpt.h"
//#include "EEF1Ligand.h"//for ligand parameterization


/////////////////////////////////////////////////////////////////////////////
//GLOBAL VARIABLES AND PARAMETERS
/////////////////////////////////////////////////////////////////////////////

//global variables, file paths
char PROGRAM_PATH[MAX_LENGTH_ONE_LINE_IN_FILE+1]=".";
char PROGRAM_NAME[MAX_LENGTH_FILE_NAME+1]="EvoEF2";

//target design protein/ligand 
char PDB[MAX_LENGTH_FILE_NAME+1] = "protein.pdb";
char PDBPATH[MAX_LENGTH_ONE_LINE_IN_FILE+1]=".";
char PDBNAME[MAX_LENGTH_FILE_NAME+1];
char PDBID[MAX_LENGTH_FILE_NAME+1];
char DES_CHAINS[10] = "";

//default parameters for monomer evolution energy
char TGT_PRF[MAX_LENGTH_FILE_NAME+1]="prf.txt";
char TGT_SA[MAX_LENGTH_FILE_NAME+1]="sa.txt";
char TGT_SS[MAX_LENGTH_FILE_NAME+1]="ss.txt";
char TGT_SEQ[MAX_LENGTH_FILE_NAME+1]="seq.txt";
char TGT_PHIPSI[MAX_LENGTH_FILE_NAME+1]="phipsi.txt";

//default parameters for physical energy -> protein
char FILE_AMINOATOMPAR[MAX_LENGTH_FILE_NAME+1] = "library/param_charmm19_lk.prm";
char FILE_AMINOTOP[MAX_LENGTH_FILE_NAME+1]     = "library/top_polh19_prot.inp";
char FILE_ROTLIB[MAX_LENGTH_FILE_NAME+1] = "library/dun2010bb3per.lib";
char FILE_AAPROPENSITY[MAX_LENGTH_FILE_NAME+1] = "library/aapropensity.txt";
char FILE_RAMACHANDRAN[MAX_LENGTH_FILE_NAME+1] = "library/ramachandran.txt";
char FILE_WEIGHT_READ[MAX_LENGTH_FILE_NAME+1]= "wread/weight_EvoEF2.txt";

//other user-specified parameters
char MUTANT_FILE[MAX_LENGTH_FILE_NAME+1] = "individual_list.txt";
double PPI_DIST_CUTOFF = 5.0;
double TORSION_DEVIATION_CUTOFF=20.0;
char INTERFACE_RESIDUE_FILE[MAX_LENGTH_FILE_NAME+1]="interfaceresidue.txt";

//design result files
char FILE_SELF_ENERGY[MAX_LENGTH_FILE_NAME+1]="selfenergy.txt";
char ROT_LIST_ALL[MAX_LENGTH_FILE_NAME+1]="rotlist.txt";
char ROT_LIST_SEC[MAX_LENGTH_FILE_NAME+1]="rotlistsec.txt";
char ROTINDEXFILE[MAX_LENGTH_FILE_NAME+1]="desrots.txt";
char SEQFILE[MAX_LENGTH_FILE_NAME+1]="desseqs.txt";
char BEST_SEQ[MAX_LENGTH_FILE_NAME+1]="bestseq.txt";
char BEST_STRUCT[MAX_LENGTH_FILE_NAME+1]="beststruct.pdb";


//for computing binding affinity of multi-chain protein complexes
char splitchains[MAX_LENGTH_FILE_NAME+1]="";
char split_part1[MAX_LENGTH_FILE_NAME+1]="";
char split_part2[MAX_LENGTH_FILE_NAME+1]="";
BOOL chainSplitFlag=FALSE;


#define DEAL_WITH_FLAGS
/////////////////////////////////////////////////////////////////////////
//DESIGN TYPES AND FLAGS
/////////////////////////////////////////////////////////////////////////
//scoring type
BOOL FLAG_EVOLUTION=FALSE;
BOOL FLAG_PHYSICS=TRUE;

//design type
BOOL FLAG_MONOMER=TRUE;
BOOL FLAG_PPI=FALSE;
BOOL FLAG_WILDTYPE_ONLY=FALSE;

//backbone dependent potential
BOOL FLAG_BBDEP_ROTLIB=TRUE;
//add crystal residue as a rotamer
BOOL FLAG_ADD_CRYSTAL_ROT=FALSE;
//expand -OH rotmaers (SER/THR/TYR)
BOOL FLAG_EXPAND_HYDROXYL_ROT=TRUE;

//important parameters
int TOT_SEQ_LEN=1;

//user specified rotamer library
BOOL FLAG_USER_ROTLIB=FALSE;
char USER_ROTLIB_NAME[MAX_LENGTH_FILE_NAME+1]="";

//options supported in EvoDesign
const char* short_opts = "hv";
struct option long_opts[] = {
  {"help",          no_argument,       NULL,  1},
  {"version",       no_argument,       NULL,  2},
  {"command",       required_argument, NULL,  3},
  {"pdb",           required_argument, NULL,  4},
  {"split_chains",  required_argument, NULL,  5},
  {"mutant_file",   required_argument, NULL,  6},
  {"bbdep",         required_argument, NULL,  7},
  {"use_input_sc",  required_argument, NULL,  8},
  {"ppi_cutoff",    required_argument, NULL,  9},
  {"evolution",     no_argument,       NULL, 11},
  {"physics",       no_argument,       NULL, 12},
  {"monomer",       no_argument,       NULL, 13},
  {"ppint",         no_argument,       NULL, 14},
  {"design_chains", required_argument, NULL, 17},
  {"rotate_hydroxyl",required_argument,NULL, 19},
  {"xdeviation",    required_argument, NULL, 20},
  {"wread",        required_argument,  NULL, 21},
  {"wwrite",       required_argument,  NULL, 22},
  {"wildtype_only", required_argument, NULL, 23},
  {"rotlib",       required_argument,  NULL, 24},
  {NULL,            no_argument,       NULL, 0}
};


BOOL CheckCommandName(char* queryname){
  int MAX_CMD_NUM = 100;
  char *supportedcmd[] = {
    "RepairStructure", 
    "ComputeStability", 
    "ComputeBinding", 
    "BuildMutant",
    "ComputeResiEnergy",
    "ComputeRotamersEnergy",
    "OptimizeHydrogen",
    "AddHydrogen",
    "FindInterfaceResidue",
    "ShowResiComposition",
    "ProteinDesign",
    "OptimizeWeight",
    "GetPhiPsi",
    "TestBBdepRotLib",
    "CheckResiInLab",
    "SideChainRepack",
    "GetResiMinRMSD",
    NULL
  };

  BOOL exist = FALSE;
  for(int i = 0; i < MAX_CMD_NUM; i++){
    if(supportedcmd[i] == NULL) break;
    else{
      if(strcmp(queryname, supportedcmd[i]) == 0){
        exist = TRUE;
        break;
      }
    }
  }
  return exist;
}


int ExtractPathAndName(char* fullpath, char* path, char* name){
  int i=0;
  BOOL slash=FALSE;
  for(i=strlen(fullpath); i>=0; i--){
    if(fullpath[i]=='/' || fullpath[i]=='\\'){
      strcpy(name,fullpath+i+1);
      strncpy(path, fullpath, i);
      path[i]='\0';
      slash=TRUE;
      break;
    }
  }
  if(slash==FALSE){
    strcpy(name,fullpath);
    strcpy(path,".");
  }
  return Success;
}

int GetPDBID(char* pdbfile, char* pdbid){
  for(int i=0; i<(int)strlen(pdbfile);i++){
    if(pdbfile[i]=='.'){
      strncpy(pdbid,pdbfile,i);
      pdbid[i]='\0';
      break;
    }
  }
  return Success;
}


int GenerateProfile(char* programpath,char* pdb){
  char cmd[MAX_LENGTH_ONE_LINE_IN_FILE+1]="";
  sprintf(cmd,"%s/extbin/profile/GenerateProfile.pl %s",programpath,pdb);
  system(cmd);
  return Success;
}

int GenerateSingleSequenceFeatures(char* programpath,char* pdb){
  char cmd[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  printf("secondary structure prediction of target will be done by DSSP\n");
  printf("running DSSP ... \n");
  sprintf(cmd,"%s/extbin/dsspcmbi %s > dssp.txt 2>/dev/null",programpath,pdb);
  system(cmd);
  printf("done\n");
  printf("extract seondary structure information ... \n");
  sprintf(cmd,"%s/extbin/dssp-ver2hor dssp.txt > dssph.txt",programpath);
  system(cmd);
  printf("done\n");
  printf("convert to I-TASSER like secondary structure format ... \n");
  sprintf(cmd,"%s/extbin/convSS.pl dssph.txt > ss.txt",programpath);
  system(cmd);
  printf("done\n");
  printf("generated: secondary structure, surface's solvent-accessibility, and phi-psi files.\n");
  return Success;
}


int main(int argc, char* argv[]){
  clock_t timestart = clock();
  setvbuf(stdout, NULL, _IONBF, 0);
  //deal with filepaths related to program EVOEF
  //you MUST leave this paragraph here to make sure the EvoEF2 program can find the correct path and files
  ExtractPathAndName(argv[0], PROGRAM_PATH, PROGRAM_NAME);
  sprintf(FILE_AMINOATOMPAR,"%s/library/param_charmm19_lk.prm",PROGRAM_PATH);
  sprintf(FILE_AMINOTOP,"%s/library/top_polh19_prot.inp",PROGRAM_PATH);
  sprintf(FILE_ROTLIB,"%s/library/bbdep3per.lib",PROGRAM_PATH);
  sprintf(FILE_AAPROPENSITY,"%s/library/aapropensity.txt",PROGRAM_PATH);
  sprintf(FILE_RAMACHANDRAN,"%s/library/ramachandran.txt",PROGRAM_PATH);
  sprintf(FILE_WEIGHT_READ,"%s/wread/weight_EvoEF2.txt",PROGRAM_PATH);

  //show interface
  EVOEF_interface();
  
  char* cmdname="";
  //char *cmdname = "RepairStructure";
  //char *cmdname = "ComputeStability";
  //char *cmdname = "ComputeBinding";
  //char *cmdname = "BuildMutant";
  //char *cmdname = "OptimizeHydrogen";
  //char *cmdname = "ComputeResiEnergy";
  //char *cmdname = "AddHydrogen";
  //char *cmdname = "FindInterfaceResidue";
  //char *cmdname = "ComputeRotamersEnergy";
  //char *cmdname = "ShowResiComposition";
  //char *cmdname = "ProteinDesign";
  //char *cmdname = "OptimizeWeight";
  //char* cmdname = "GetPhiPsi";
  //char* cmdname = "CheckResiInLab";
  //char* cmdname = "SideChainRepack";
  //char* cmdname = "GetResiMinRMSD";

  while(TRUE){
    int opt=getopt_long(argc, argv, short_opts, long_opts, NULL);
    if(opt == -1){
      break;
    }
    switch(opt){
      //deal with short options
      case 'h':
        EVOEF_help();
        exit(Success);
      case 'v':
        EVOEF_version();
        exit(Success);
        // deal with long options
      case 1:
        EVOEF_help();
        exit(Success);
      case 2:
        EVOEF_version();
        exit(Success);
      case 3:
        cmdname = optarg;
        if(!CheckCommandName(cmdname)){
          printf("command %s is not supported, EVOEF exits\n", cmdname);
          exit(ValueError);
        }
        else printf("command %s works\n", cmdname);
        break;
      case 4:
        strcpy(PDB,optarg);
        break;
      case 5:
        strcpy(splitchains,optarg);
        sscanf(splitchains,"%[^,],%s",split_part1,split_part2);
        chainSplitFlag=TRUE;
        //check if the two parts contain identical chains
        for(int i=0; i<(int)strlen(split_part1); i++){
          char tmp[2]={split_part1[i],'\0'};
          if(strstr(split_part2,tmp)!=NULL){
            printf("the split two parts contain identical chains,please check, EVOEF exits\n");
            exit(FormatError);
          }
        }
        break;
      case 6:
        strcpy(MUTANT_FILE,optarg);
        break;
      case 7:
        if(strcmp(optarg,"enable")==0){
          printf("backbone-dependent rotamer library enabled by user\n");
          FLAG_BBDEP_ROTLIB=TRUE;
        }
        else if(strcmp(optarg,"disable")==0){
          printf("backbone-independent rotamer library enabled by user\n");
          FLAG_BBDEP_ROTLIB=FALSE;
        }
        else{
          printf("backbone-dependent rotamer library enabled by default\n");
          FLAG_BBDEP_ROTLIB=TRUE;
        }
        break;
      case 8:
        if(strcmp(optarg,"enable")==0) FLAG_ADD_CRYSTAL_ROT=TRUE;
        else if(strcmp(optarg,"disable")==0) FLAG_ADD_CRYSTAL_ROT=FALSE;
        else FLAG_ADD_CRYSTAL_ROT=TRUE;
        break;
      case 9:
        PPI_DIST_CUTOFF = atof(optarg);
        break;
      case 11:
        FLAG_EVOLUTION = TRUE;
        break;
      case 12:
        FLAG_PHYSICS = TRUE;
        break;
      case 13:
        FLAG_MONOMER  = TRUE;
        FLAG_PPI      = FALSE;
        break;
      case 14:
        FLAG_PPI      = TRUE;
        FLAG_MONOMER  = FALSE;
        break;
      case 17:
        strcpy(DES_CHAINS,optarg);
        break;
      case 19:
        if(strcmp(optarg,"enable")==0) FLAG_EXPAND_HYDROXYL_ROT=TRUE;
        else if(strcmp(optarg,"disable")==0) FLAG_EXPAND_HYDROXYL_ROT=FALSE;
        else FLAG_EXPAND_HYDROXYL_ROT=TRUE;
        break;
      case 20:
        TORSION_DEVIATION_CUTOFF=atof(optarg);
        break;
      case 21:
        strcpy(FILE_WEIGHT_READ,optarg);
        break;
      case 23:
        FLAG_WILDTYPE_ONLY=TRUE;
      case 24:
        FLAG_USER_ROTLIB=TRUE;
        strcpy(USER_ROTLIB_NAME,optarg);
        break;
      default:
        char errMsg[MAX_LENGTH_ERR_MSG+1];
        sprintf(errMsg, "in file %s function %s() line %d, unknown option, EVOEF exits", __FILE__, __FUNCTION__, __LINE__);
        TraceError(errMsg, ValueError);
        exit(ValueError);
        break;
    }
  }

  // Early termination if no command detected
  if(!strcmp(cmdname,"")){
    printf("no command name detected, please select a command for execution, program exits\n");
    printf("run './%s -h' to show help document\n",PROGRAM_NAME);
    exit(Success);
  }


  //specify the rotamer library
  if(FLAG_BBDEP_ROTLIB==TRUE){
    sprintf(FILE_ROTLIB,"%s/library/dun2010bb3per.lib",PROGRAM_PATH);
  }
  else{
    sprintf(FILE_ROTLIB,"%s/library/honig984.lib",PROGRAM_PATH);
  }
  if(FLAG_USER_ROTLIB==TRUE){
    sprintf(FILE_ROTLIB,"%s/library/%s.lib",PROGRAM_PATH,USER_ROTLIB_NAME);
    if(strstr(USER_ROTLIB_NAME,"honig")!=NULL){
      FLAG_BBDEP_ROTLIB=FALSE;
    }
    else{
      FLAG_BBDEP_ROTLIB=TRUE;
    }
  }

  // deal with pdbid and filenames
  ExtractPathAndName(PDB, PDBPATH, PDBNAME);
  GetPDBID(PDBNAME, PDBID);
  sprintf(FILE_SELF_ENERGY,"%s_selfenergy.txt",PDBID);
  sprintf(ROT_LIST_ALL,"%s_rotlist.txt",PDBID);
  sprintf(ROT_LIST_SEC,"%s_rotlistSEC.txt",PDBID);
  sprintf(ROTINDEXFILE,"%s_desrots.txt",PDBID);
  sprintf(SEQFILE,"%s_desseqs.txt",PDBID);
  sprintf(BEST_SEQ,"%s_bestseq.txt",PDBID);
  sprintf(BEST_STRUCT,"%s_beststruct.pdb",PDBID);

  //read in data and pre-process
  AtomParamsSet atomParam;
  ResiTopoSet resiTopo;
  Structure structure;
  AtomParamsSetCreate(&atomParam);
  ResiTopoSetCreate(&resiTopo);
  AtomParameterRead(&atomParam, FILE_AMINOATOMPAR);
  ResiTopoSetRead(&resiTopo, FILE_AMINOTOP);
  StructureCreate(&structure);
  //read in protein scaffold with PDB format
  StructureConfig(&structure, PDB, &atomParam, &resiTopo);

  ///////////////////////////////////
  //read in the energy weights
  ///////////////////////////////////
  EnergyWeightRead(FILE_WEIGHT_READ);

  //set whole sequence length for energy normalization
  int resiCount=0;
  for(int i=0;i<StructureGetChainCount(&structure);i++){
    if(ChainGetType(StructureGetChain(&structure,i))==Type_Chain_Protein)
      resiCount+=ChainGetResidueCount(StructureGetChain(&structure,i));
  }
  TOT_SEQ_LEN=resiCount>0 ? resiCount:100;

  //set design chains
  if(strcmp(DES_CHAINS,"")==0){
    strcpy(DES_CHAINS,ChainGetName(StructureGetChain(&structure,0)));
  }
  else{
    printf("protein chains %s were selected for design by default\n",DES_CHAINS);
  }



  //calculate phi and psi angles for protein residues
  StructureCalcPhiPsi(&structure);
  //calculate protein sidechain dihedral angles
  StructureCalcProteinResidueSidechainTorsion(&structure,&resiTopo);

  //generate evolution profile
  if(FLAG_EVOLUTION==TRUE){
    GenerateProfile(PROGRAM_PATH,PDB);
    GenerateSingleSequenceFeatures(PROGRAM_PATH,PDB);
  }

  //deal with the specified command
  if(!strcmp(cmdname, "ComputeStability")){
    AAppTable aapptable;
    RamaTable ramatable;
    AApropensityTableReadFromFile(&aapptable,FILE_AAPROPENSITY);
    RamaTableReadFromFile(&ramatable,FILE_RAMACHANDRAN);
    if(FLAG_BBDEP_ROTLIB==TRUE){
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate(&bbrotlib,FILE_ROTLIB);
      double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
      EVOEF_ComputeStabilityWithBBdepRotLib(&structure,&aapptable,&ramatable,&bbrotlib,energyTerms);
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else{
      double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
      EVOEF_ComputeStability(&structure,&aapptable,&ramatable,energyTerms);
    }
  }
  else if(!strcmp(cmdname, "ComputeBinding")){
    if(StructureGetChainCount(&structure)<=2){
      EVOEF_ComputeBinding(&structure);
    }
    else{
      if(chainSplitFlag==FALSE){
        EVOEF_ComputeBinding(&structure);
      }
      else{
        //EVOEF_ComputeBindingWithSplitting(&structure,split_part1,split_part2);
        EVOEF_ComputeBindingWithSplittingNew(&structure,split_part1,split_part2);
      }
    }
  }
  else if(!strcmp(cmdname, "RepairStructure")){
    AAppTable aapptable;
    RamaTable ramatable;
    AApropensityTableReadFromFile(&aapptable,FILE_AAPROPENSITY);
    RamaTableReadFromFile(&ramatable,FILE_RAMACHANDRAN);
    if(FLAG_BBDEP_ROTLIB==TRUE){
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate(&bbrotlib,FILE_ROTLIB);
      EVOEF_RepairStructureWithBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo,PDBID);
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else{
      BBindRotamerLib rotlib;
      BBindRotamerLibCreate(&rotlib,FILE_ROTLIB);
      EVOEF_RepairStructure(&structure, &rotlib, &atomParam, &resiTopo,PDBID);
      BBindRotamerLibDestroy(&rotlib);
    }
  }
  else if(!strcmp(cmdname, "BuildMutant")){
    AAppTable aapptable;
    RamaTable ramatable;
    AApropensityTableReadFromFile(&aapptable,FILE_AAPROPENSITY);
    RamaTableReadFromFile(&ramatable,FILE_RAMACHANDRAN);
    if(FLAG_BBDEP_ROTLIB==TRUE){
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate(&bbrotlib,FILE_ROTLIB);
      EVOEF_BuildMutantWithBBdepRotLib(&structure,MUTANT_FILE,&bbrotlib,&atomParam,&resiTopo,PDBID);
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else{
      BBindRotamerLib rotlib;
      BBindRotamerLibCreate(&rotlib,FILE_ROTLIB);
      EVOEF_BuildMutant(&structure, MUTANT_FILE, &rotlib, &atomParam, &resiTopo,PDBID);
      BBindRotamerLibDestroy(&rotlib);
    }
  }
  else if(!strcmp(cmdname,"ShowResiComposition")){
    int aas[20]={0};
    StructureGetAminoAcidComposition(&structure,aas);
  }
  else if(!strcmp(cmdname,"OptimizeHydrogen")){
    EVOEF_OptimizeHydrogen(&structure,&atomParam, &resiTopo,PDBID);
  }
  else if(!strcmp(cmdname,"AddHydrogen")){
    EVOEF_AddHydrogen(&structure,PDBID);
  }
  else if(!strcmp(cmdname,"FindInterfaceResidue")){
    sprintf(INTERFACE_RESIDUE_FILE,"%s_interfaceresidue.txt",PDBID);
    EVOEF_StructureFindInterfaceResidues(&structure,PPI_DIST_CUTOFF,INTERFACE_RESIDUE_FILE);
  }
  else if(!strcmp(cmdname,"SideChainRepack")){
    if(FLAG_BBDEP_ROTLIB==TRUE){
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate(&bbrotlib,FILE_ROTLIB);
      StructureGenerateWildtypeRotamersByBBdepRotLib(&structure,&bbrotlib,&atomParam,&resiTopo,DES_CHAINS);
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else{
      BBindRotamerLib rotlib;
      BBindRotamerLibCreate(&rotlib, FILE_ROTLIB);
      StructureGenerateWildtypeRotamers(&structure,&rotlib,&atomParam,&resiTopo,DES_CHAINS);
      BBindRotamerLibDestroy(&rotlib);
    }
    AAppTable aapptable;
    RamaTable ramatable;
    AApropensityTableReadFromFile(&aapptable,FILE_AAPROPENSITY);
    RamaTableReadFromFile(&ramatable,FILE_RAMACHANDRAN);
    StructureShowDesignSites(&structure, stdout);
    SelfEnergyGenerateWithBBdepEnergy(&structure,&aapptable,&ramatable,FILE_SELF_ENERGY);
    RotamerList rotList;
    RotamerListCreateFromStructure(&rotList, &structure);
    RotamerListWrite(&rotList,ROT_LIST_ALL);
    SelfEnergyReadAndCheck(&structure,&rotList,FILE_SELF_ENERGY);
    RotamerListWrite(&rotList,ROT_LIST_SEC);
    RotamerListRead(&rotList,ROT_LIST_SEC);
    StructureShowDesignSitesAfterRotamerDelete(&structure,&rotList,stdout);
    SimulatedAnnealingOptimizationForSCP(&structure,&rotList,DES_CHAINS,PROGRAM_PATH,TGT_PRF,TGT_SS,TGT_SA,TGT_PHIPSI,ROTINDEXFILE,SEQFILE);
  }

  else if(!strcmp(cmdname, "ProteinDesign")){
    if(FLAG_BBDEP_ROTLIB==TRUE){
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate(&bbrotlib,FILE_ROTLIB);
      StructureGenerateAllRotamersByBBdepRotLib(&structure,&bbrotlib,&atomParam,&resiTopo,DES_CHAINS);
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else{
      BBindRotamerLib rotlib;
      BBindRotamerLibCreate(&rotlib, FILE_ROTLIB);
      StructureGenerateAllRotamers(&structure,&rotlib,&atomParam,&resiTopo,DES_CHAINS);
      BBindRotamerLibDestroy(&rotlib);
    }

    AAppTable aapptable;
    RamaTable ramatable;
    AApropensityTableReadFromFile(&aapptable,FILE_AAPROPENSITY);
    RamaTableReadFromFile(&ramatable,FILE_RAMACHANDRAN);
    StructureShowDesignSites(&structure, stdout);
    SelfEnergyGenerateWithBBdepEnergy(&structure,&aapptable,&ramatable,FILE_SELF_ENERGY);
    RotamerList rotList;
    RotamerListCreateFromStructure(&rotList, &structure);
    RotamerListWrite(&rotList,ROT_LIST_ALL);
    SelfEnergyReadAndCheck(&structure,&rotList,FILE_SELF_ENERGY);
    RotamerListWrite(&rotList,ROT_LIST_SEC);
    RotamerListRead(&rotList,ROT_LIST_SEC);
    StructureShowDesignSitesAfterRotamerDelete(&structure,&rotList,stdout);
    SimulatedAnnealingOptimization(&structure,&rotList,DES_CHAINS,PROGRAM_PATH,TGT_PRF,TGT_SS,TGT_SA,TGT_PHIPSI,ROTINDEXFILE,SEQFILE);
  }
  else{
    printf("unknown command name: '%s', try './EvoEF2 -h' to get help. Exits\n", cmdname);
    exit(ValueError);
  }

  StructureDestroy(&structure);
  ResiTopoSetDestroy(&resiTopo);
  AtomParamsSetDestroy(&atomParam);

  clock_t timeend = clock();
  SpentTimeShow(timestart,timeend);
  
  return Success;
}



