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
#include "DEE.h"
#include "EnergyOptimization.h"
#include "WeightOpt.h"
#include "EEF1Ligand.h"//for ligand parameterization


/////////////////////////////////////////////////////////////////////////////
//GLOBAL VARIABLES AND PARAMETERS
/////////////////////////////////////////////////////////////////////////////

//global variables, file paths
char PROGRAM_PATH[MAX_LENGTH_ONE_LINE_IN_FILE+1]=".";
char PROGRAM_NAME[MAX_LENGTH_FILE_NAME+1]="EvoEF2";

//target design protein/ligand 
char PDB[MAX_LENGTH_FILE_NAME+1] = "example/1aac.pdb";
char PDBPATH[MAX_LENGTH_ONE_LINE_IN_FILE+1]=".";
char PDBNAME[MAX_LENGTH_FILE_NAME+1];
char PDBID[MAX_LENGTH_FILE_NAME+1];
char MOL2[MAX_LENGTH_FILE_NAME+1]="example/ligand.mol2";
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
char FILE_WEIGHT_WRITE[MAX_LENGTH_FILE_NAME+1]= "wwrite/weightout.txt";

//temp parameters for physical energy -> ligand
char LIG_ATOMPAR[MAX_LENGTH_FILE_NAME+1]="ligatompar.txt";
char LIG_TOPFILE[MAX_LENGTH_FILE_NAME+1]="ligtoph.txt";


//other user-specified parameters
char MUTANT_FILE[MAX_LENGTH_FILE_NAME+1] = "individual_list.txt";
double BIND_RADIUS = 7.5;
double PPI_DIST_CUTOFF = 5.0;
double TORSION_DEVIATION_CUTOFF=20.0;
char INTERFACE_RESIDUE_FILE[MAX_LENGTH_FILE_NAME+1]="interfaceresidue.txt";

//parameters for enzyme design
char FILE_SPECIFIED_SITES[MAX_LENGTH_FILE_NAME+1]="design/enzyme/IT1647/sites_catalytic.txt";
char FILE_SITES_CATALYTIC[MAX_LENGTH_FILE_NAME+1]="design/enzyme/IT1647/sites_catalytic.txt";
char FILE_SITES_MUTATED[MAX_LENGTH_FILE_NAME+1]="design/enzyme/IT1647/sites_mutated.txt";
char FILE_SITES_ROTAMERIC[MAX_LENGTH_FILE_NAME+1]="design/enzyme/IT1647/sites_rotameric.txt";
double dist_cutoff_shell1=5.0;
double dist_cutoff_shell2=7.5;


char FILE_CATALYTIC_CONSTRAINT[MAX_LENGTH_FILE_NAME+1]="design/enzyme/IT1647/catacons.txt";
char FILE_LIG_PLACEMENT[MAX_LENGTH_FILE_NAME+1]="design/enzyme/IT1647/ligplacement.txt";
char FILE_LIG_ENSEMBLE[MAX_LENGTH_FILE_NAME+1]="design/enzyme/IT1647/ligands.pdb";
char FILE_LIG_ENSEMBLE_TOPVDW[MAX_LENGTH_FILE_NAME+1]="design/enzyme/IT1647/ligands_top.pdb";
char FILE_LIG_ENSEMBLE_RMSDSCREEN[MAX_LENGTH_FILE_NAME+1]="design/enzyme/IT1647/ligands_rmsd.pdb";



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
BOOL FLAG_PROT_LIG=FALSE;
BOOL FLAG_ENZYME=FALSE;
BOOL FLAG_WILDTYPE_ONLY=FALSE;

//backbone dependent potential
BOOL FLAG_BBDEP_ROTLIB=TRUE;
//add crystal residue as a rotamer
BOOL FLAG_ADD_CRYSTAL_ROT=FALSE;
//expand -OH rotmaers (SER/THR/TYR)
BOOL FLAG_EXPAND_HYDROXYL_ROT=TRUE;

//important parameters
int TOT_SEQ_LEN=100;

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
  {"split",         required_argument, NULL,  5},
  {"mutant_file",   required_argument, NULL,  6},
  {"bbdep",         required_argument, NULL,  7},
  {"addcrystal",    required_argument, NULL,  8},
  {"ppi_cutoff",    required_argument, NULL,  9},
  {"bind_radius",   required_argument, NULL, 10},
  {"evolution",     no_argument,       NULL, 11},
  {"physics",       no_argument,       NULL, 12},
  {"monomer",       no_argument,       NULL, 13},
  {"ppi",           no_argument,       NULL, 14},
  {"protlig",       no_argument,       NULL, 15},
  {"enzyme",        no_argument,       NULL, 16},
  {"design_chains", required_argument, NULL, 17},
  {"mol2",          required_argument, NULL, 18},
  {"expand_hydroxyl",required_argument,NULL, 19},
  {"xdeviation",    required_argument, NULL, 20},
  {"wread",        required_argument,  NULL, 21},
  {"wwrite",       required_argument,  NULL, 22},
  {"wildtype",     required_argument,  NULL, 23},
  {"rotlib",       required_argument,  NULL, 24},
  {NULL,            no_argument,       NULL, 0}
};


BOOL CheckCommandName(char* queryname);
int ExtractPathAndName(char* fullpath, char* path, char* name);
int GetPDBID(char* pdbfile, char* pdbid);
int GenerateProfile(char* programpath,char* pdb);
int GenerateSingleSequenceFeatures(char* programpath,char* pdb);
int GetAttachedAtomTypeNum(AtomArray* pAtoms,int bonds[][2],int bondNum,int atomIndex,int attachAtomTypeNum[10],int *attachCount,int attachIndex[10]);
int GenerateParameterAndTopologyFromMol2(char* mol2file,char* parfile,char* topfile);


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
  
  //char *cmdname = "RepairStructure";
  //char *cmdname = "ComputeStability";
  //char *cmdname = "ComputeBinding";
  char *cmdname = "BuildMutant";
  //char *cmdname = "OptimizeHydrogen";
  //char *cmdname = "ComputeResiEnergy";
  //char *cmdname = "AddHydrogen";
  //char *cmdname = "FindInterfaceResidue";
  //char *cmdname = "ComputeRotamersEnergy";
  //char *cmdname = "ComputeWtRotamersEnergy";
  //char *cmdname = "ShowResiComposition";
  //char *cmdname = "MakeLigandEnsemble";
  //char *cmdname = "ScreenLigandEnsemble";
  //char *cmdname = "ComputeEnergyMatrix";
  //char *cmdname = "ComputeEvolutionScore";
  //char *cmdname = "ProteinDesign";
  //char *cmdname = "ProtLigDesign";
  //char *cmdname = "EnzymeDesign";
  //char *cmdname = "OptimizeWeight";
  //char* cmdname = "GetPhiPsi";
  //char* cmdname = "TestBBdepRotLib";
  //char* cmdname = "CheckResiInLab";
  //char* cmdname = "SidechainRepack";
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
      case 10:
        BIND_RADIUS = atof(optarg);
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
        FLAG_PROT_LIG = FALSE;
        FLAG_ENZYME   = FALSE;
        break;
      case 14:
        FLAG_PPI      = TRUE;
        FLAG_MONOMER  = FALSE;
        FLAG_PROT_LIG = FALSE;
        FLAG_ENZYME   = FALSE;
        break;
      case 15:
        FLAG_PROT_LIG = TRUE;
        FLAG_PPI      = FALSE;
        FLAG_MONOMER  = FALSE;
        FLAG_ENZYME   = FALSE;
        break;
      case 16:
        FLAG_ENZYME   = TRUE;
        FLAG_PROT_LIG = FALSE;
        FLAG_PPI      = FALSE;
        FLAG_MONOMER  = FALSE;
        break;
      case 17:
        strcpy(DES_CHAINS,optarg);
        break;
      case 18:
        strcpy(MOL2,optarg);
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
      case 22:
        strcpy(FILE_WEIGHT_WRITE,optarg);
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
  sprintf(LIG_ATOMPAR,"%s_ligatompar.txt",PDBID);
  sprintf(LIG_TOPFILE,"%s_ligtoph.txt",PDBID);
  sprintf(FILE_WEIGHT_WRITE,"%s_weight.txt",PDBID);


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
  printf("protein structure file %s.pdb was successfully read by EVOEF\n", PDBID);
  //read in small molecules: generate parameters and read structure with mol2 format
  if(FLAG_PROT_LIG==TRUE || FLAG_ENZYME==TRUE){
    GenerateParameterAndTopologyFromMol2(MOL2,LIG_ATOMPAR,LIG_TOPFILE);
    AtomParameterRead(&atomParam,LIG_ATOMPAR);
    ResiTopoSetRead(&resiTopo,LIG_TOPFILE);
    StructureConfigLigand(&structure,MOL2,&atomParam,&resiTopo);
  }

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
    printf("no chains was designated for design, the 1st chain %s was selected by default\n",DES_CHAINS);
  }
  else{
    printf("protein chains %s were selected for design by default\n",DES_CHAINS);
  }

  printf("current rotlib: %s\n", FILE_ROTLIB);


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
  else if(!strcmp(cmdname, "ComputeResiEnergy")){
    for(int i=0; i<StructureGetChainCount(&structure); ++i){
      Chain* pChain=StructureGetChain(&structure,i);
      for(int j=0; j<ChainGetResidueCount(pChain); ++j){
        Residue* pResidue=ChainGetResidue(pChain,j);
        printf("residue %s%d%s energy details:\n",ChainGetName(pChain), ResidueGetPosInChain(pResidue), ResidueGetName(pResidue));
        EVOEF_StructureComputeResidueInteractionWithFixedSurroundingResidues(&structure, i, j);
      }
    }
  }
  else if(!strcmp(cmdname, "ComputeRotamersEnergy")){
    AAppTable aapptable;
    RamaTable ramatable;
    AApropensityTableReadFromFile(&aapptable,FILE_AAPROPENSITY);
    RamaTableReadFromFile(&ramatable,FILE_RAMACHANDRAN);
    if(FLAG_BBDEP_ROTLIB==TRUE){
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate(&bbrotlib,FILE_ROTLIB);
      if(FLAG_WILDTYPE_ONLY==TRUE){
        EVOEF_ComputeWildtypeRotamersEnergyByBBdepRotLib(&structure,&bbrotlib,&aapptable,&ramatable,&atomParam,&resiTopo,PDBID);
      }
      else{
        EVOEF_ComputeRotamersEnergyByBBdepRotLib(&structure,&bbrotlib,&aapptable,&ramatable,&atomParam,&resiTopo,PDBID);
      }
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else{
      BBindRotamerLib rotlib;
      BBindRotamerLibCreate(&rotlib,FILE_ROTLIB);
      if(FLAG_WILDTYPE_ONLY==TRUE){
        EVOEF_ComputeWildtypeRotamersEnergy(&structure,&rotlib,&aapptable,&ramatable,&atomParam,&resiTopo,PDBID);
      }
      else{
        EVOEF_ComputeRotamersEnergy(&structure,&rotlib,&aapptable,&ramatable,&atomParam,&resiTopo,PDBID);
      }
      BBindRotamerLibDestroy(&rotlib);
    }
  }
  else if(!strcmp(cmdname,"ShowResiComposition")){
    int aas[20]={0};
    StructureGetAminoAcidComposition(&structure,aas);
  }
  else if(!strcmp(cmdname,"OptimizeHydrogen")){
    EVOEF_OptimizeHydrogen(&structure,&atomParam, &resiTopo,PDBNAME);
  }
  else if(!strcmp(cmdname,"AddHydrogen")){
    EVOEF_AddHydrogen(&structure,PDBID);
  }
  else if(!strcmp(cmdname,"FindInterfaceResidue")){
    sprintf(INTERFACE_RESIDUE_FILE,"%s_interfaceresidue.txt",PDBID);
    EVOEF_StructureFindInterfaceResidues(&structure,PPI_DIST_CUTOFF,INTERFACE_RESIDUE_FILE);
  }
  else if(!strcmp(cmdname, "MakeLigandEnsemble")){
    if(FLAG_BBDEP_ROTLIB==TRUE){
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate(&bbrotlib,FILE_ROTLIB);
      StructureGenerateSpecifiedProteinRotamersByBBdepRotLib(&structure,&bbrotlib,&atomParam,&resiTopo,FILE_SPECIFIED_SITES);
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else{
      BBindRotamerLib rotlib;
      BBindRotamerLibCreate(&rotlib, FILE_ROTLIB);
      StructureGenerateSpecifiedProteinRotamers(&structure,&rotlib,&atomParam,&resiTopo,FILE_SPECIFIED_SITES);
      BBindRotamerLibDestroy(&rotlib);
    }
    AAppTable aapptable;
    RamaTable ramatable;
    AApropensityTableReadFromFile(&aapptable,FILE_AAPROPENSITY);
    RamaTableReadFromFile(&ramatable,FILE_RAMACHANDRAN);
    StructureShowDesignSites(&structure, stdout);
    StructureGenerateSmallMolRotamers(&structure,FILE_CATALYTIC_CONSTRAINT,FILE_LIG_PLACEMENT);
    StructureWriteSmallMolRotamers(&structure,FILE_LIG_ENSEMBLE);
    StructureShowDesignSites(&structure, stdout);
  }
  else if(!strcmp(cmdname,"ScreenLigandEnsemble")){
    //StructureReadSmallMolRotamers(&structure,&resiTopo,FILE_LIG_ENSEMBLE);
    //StructureSmallmolOrientationScreen(&structure,&resiTopo,"./cel5a/ligands.pdb","./cel5a/ligands_orientaion.pdb","./cel5a/orientation.txt");
    //StructureReadSmallMolRotamers(&structure,&resiTopo,"./cel5a/ligands.pdb");
    Residue* pSmallmol=NULL;
    StructureFindSmallMol(&structure,&pSmallmol);
    //AnalyzeSmallMolRotamers(FILE_LIG_ENSEMBLE,pSmallmol);
    SmallMolRotamersGetBothHighRankOfBackboneVdwAndInternalVdw(FILE_LIG_ENSEMBLE,FILE_LIG_ENSEMBLE_TOPVDW,0.25);
    ScreenSmallmolRotamersByRMSD(FILE_LIG_ENSEMBLE_TOPVDW,FILE_LIG_ENSEMBLE_RMSDSCREEN,0.05);
    //StructureWriteSmallMolRotamers(&structure,"design/enzyme/IT1647/ligands_screen.pdb");
  }
  else if(!strcmp(cmdname, "MakeSelfEnergyMatrix")){
    if(FLAG_BBDEP_ROTLIB==TRUE){
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate(&bbrotlib,FILE_ROTLIB);
      StructureGenerateSpecifiedProteinRotamersByBBdepRotLib(&structure,&bbrotlib,&atomParam,&resiTopo,FILE_SPECIFIED_SITES);
      StructureGenerateBindingSiteRotamersByBBdepRotLib(&structure,&bbrotlib,&atomParam,&resiTopo,BIND_RADIUS);
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else{
      BBindRotamerLib rotlib;
      BBindRotamerLibCreate(&rotlib, FILE_ROTLIB);
      StructureGenerateSpecifiedProteinRotamers(&structure,&rotlib,&atomParam,&resiTopo,FILE_SPECIFIED_SITES);
      StructureGenerateBindingSiteRotamers(&structure,&rotlib,&atomParam,&resiTopo,BIND_RADIUS);
      BBindRotamerLibDestroy(&rotlib);
    }
    if(FLAG_PROT_LIG==TRUE || FLAG_ENZYME==TRUE){
      StructureReadSmallMolRotamers(&structure,&resiTopo,FILE_LIG_ENSEMBLE);
    }
    AAppTable aapptable;
    RamaTable ramatable;
    AApropensityTableReadFromFile(&aapptable,FILE_AAPROPENSITY);
    RamaTableReadFromFile(&ramatable,FILE_RAMACHANDRAN);
    StructureShowDesignSites(&structure, stdout);
    SelfEnergyGenerateWithBBdepEnergy(&structure,&aapptable,&ramatable,FILE_SELF_ENERGY);
  }
  else if(!strcmp(cmdname,"SidechainRepack")){
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
    //ProteinSiteWriteRotamers(&structure,0,17,"1acc_a17_bbdeprots.pdb");
    //ProteinSiteWriteRotamers(&structure,0,17,"1acc_a17_bbindrots.pdb");
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
    //MonteCarloOptimization(&structure,&rotList,DES_CHAINS,PROGRAM_PATH,TGT_PRF,TGT_SS,TGT_SA,TGT_PHIPSI,ROTINDEXFILE,SEQFILE);
    //ReannealingOptimization(&structure,&rotList,DES_CHAINS,PROGRAM_PATH,TGT_PRF,TGT_SS,TGT_SA,TGT_PHIPSI,ROTINDEXFILE,SEQFILE);
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
    //MonteCarloOptimization(&structure,&rotList,DES_CHAINS,PROGRAM_PATH,TGT_PRF,TGT_SS,TGT_SA,TGT_PHIPSI,ROTINDEXFILE,SEQFILE);
    //ReannealingOptimization(&structure,&rotList,DES_CHAINS,PROGRAM_PATH,TGT_PRF,TGT_SS,TGT_SA,TGT_PHIPSI,ROTINDEXFILE,SEQFILE);
    SimulatedAnnealingOptimization(&structure,&rotList,DES_CHAINS,PROGRAM_PATH,TGT_PRF,TGT_SS,TGT_SA,TGT_PHIPSI,ROTINDEXFILE,SEQFILE);
  }
  else if(!strcmp(cmdname, "ProtLigDesign")){
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
    //MonteCarloOptimization(&structure,&rotList,DES_CHAINS,PROGRAM_PATH,TGT_PRF,TGT_SS,TGT_SA,TGT_PHIPSI,ROTINDEXFILE,SEQFILE);
    //ReannealingOptimization(&structure,&rotList,DES_CHAINS,PROGRAM_PATH,TGT_PRF,TGT_SS,TGT_SA,TGT_PHIPSI,ROTINDEXFILE,SEQFILE);
    SimulatedAnnealingOptimization(&structure,&rotList,DES_CHAINS,PROGRAM_PATH,TGT_PRF,TGT_SS,TGT_SA,TGT_PHIPSI,ROTINDEXFILE,SEQFILE);
    RotamerListDestroy(&rotList);
  }
  else if(!strcmp(cmdname, "EnzymeDesign")){
    if(FLAG_BBDEP_ROTLIB==TRUE){
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate(&bbrotlib,FILE_ROTLIB);
      StructureGenerateCatalyticSiteProteinRotamersByBBdepRotLib(&structure,&bbrotlib,&atomParam,&resiTopo,FILE_SITES_CATALYTIC);
      StructureGenerateFirstShellProteinRotamersByBBdepRotLib(&structure,&bbrotlib,&atomParam,&resiTopo,FILE_SITES_MUTATED,dist_cutoff_shell1);
      StructureGenerateSecondShellProteinRotamersByBBdepRotLib(&structure,&bbrotlib,&atomParam,&resiTopo,FILE_SITES_ROTAMERIC,dist_cutoff_shell2);
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else{
      BBindRotamerLib rotlib;
      BBindRotamerLibCreate(&rotlib, FILE_ROTLIB);
      StructureGenerateSpecifiedProteinRotamers(&structure,&rotlib,&atomParam,&resiTopo,FILE_SPECIFIED_SITES);
      StructureGenerateBindingSiteRotamers(&structure,&rotlib,&atomParam,&resiTopo,BIND_RADIUS);
      BBindRotamerLibDestroy(&rotlib);
    }
    AAppTable aapptable;
    RamaTable ramatable;
    AApropensityTableReadFromFile(&aapptable,FILE_AAPROPENSITY);
    RamaTableReadFromFile(&ramatable,FILE_RAMACHANDRAN);
    //read ligand ensembles
    StructureReadSmallMolRotamers(&structure,&resiTopo,FILE_LIG_ENSEMBLE_RMSDSCREEN);
    StructureShowDesignSites(&structure, stdout);
    //SelfEnergyGenerate(&structure,SELF_ENERGY_FILE);
    SelfEnergyGenerateWithBBdepEnergy(&structure,&aapptable,&ramatable,FILE_SELF_ENERGY);
    RotamerList rotList;
    RotamerListCreateFromStructure(&rotList, &structure);
    RotamerListWrite(&rotList,ROT_LIST_ALL);
    SelfEnergyReadAndCheck(&structure,&rotList,FILE_SELF_ENERGY);
    RotamerListWrite(&rotList,ROT_LIST_SEC);
    RotamerListRead(&rotList,ROT_LIST_SEC);
    StructureShowDesignSitesAfterRotamerDelete(&structure,&rotList,stdout);
    SimulatedAnnealingOptimization(&structure,&rotList,DES_CHAINS,PROGRAM_PATH,TGT_PRF,TGT_SS,TGT_SA,TGT_PHIPSI,ROTINDEXFILE,SEQFILE);
    //ReannealingOptimization(&structure,&rotList,DES_CHAINS,PROGRAM_PATH,TGT_PRF,TGT_SS,TGT_SA,TGT_PHIPSI,ROTINDEXFILE,SEQFILE);
    RotamerListDestroy(&rotList);
  }
  else if(!strcmp(cmdname, "ComputeEvolutionScore")){
    //EvolutionScore(programpath, "test_evo/prf.txt", "test_evo/ss.txt", "test_evo/sa.txt", "test_evo/phi-psi.txt", "test_evo/seq.txt");
    char seq[]="SSAKEELLTWIQEMTAGYPNVNVSNFKTAAKDGMALCALIHKIRPDLIDFSKLKPQDALENVQMALNIAEKHLGIKMLIDPEDIVEPHPDPKSIITYLMALMHFFKQM";
    //EvolutionScore2(programpath, "test_evo/prf.txt", "test_evo/ss.txt", "test_evo/sa.txt", "test_evo/phi-psi.txt", seq);
    double evoscore = EvolutionScore3(PROGRAM_PATH, "test_evo/prf.txt",seq);
    printf("sequence: %s\nevoscore: %f\n",seq,evoscore);
  }
  else if(!strcmp(cmdname, "OptimizeWeight")){
    char *pdblist="./pdblist.txt";
    PNATAA_WeightOptByGradientDescent(pdblist);
    //PNATROT_WeightOptByGradientDescent(pdblist);
  }
  else if(!strcmp(cmdname,"TestBBdepRotLib")){
    BBdepRotamerLib bbrotlib;
    BBdepRotamerLibCreate(&bbrotlib,FILE_ROTLIB);
    BBdepRotamerLibDestroy(&bbrotlib);
  }
  else if(!strcmp(cmdname,"GetPhiPsi")){
    char PHI_PSI_FILE[MAX_LENGTH_ONE_LINE_IN_FILE+1]="";
    sprintf(PHI_PSI_FILE,"%s_phipsi.txt",PDBID);
    EVOEF_StructureShowPhiPsi(&structure,PHI_PSI_FILE);
  }
  else if(!strcmp(cmdname,"CheckResiInLab")){
    if(FLAG_BBDEP_ROTLIB==TRUE){
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate(&bbrotlib,FILE_ROTLIB);
      EVOEF_CheckRotamerInBBdepRotLib(&structure,&bbrotlib,&resiTopo,TORSION_DEVIATION_CUTOFF,PDBID);
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else{
      BBindRotamerLib rotlib;
      BBindRotamerLibCreate(&rotlib,FILE_ROTLIB);
      EVOEF_CheckRotamerInBBindRotLib(&structure,&rotlib,&resiTopo,TORSION_DEVIATION_CUTOFF,PDBID);
      BBindRotamerLibDestroy(&rotlib);
    }
  }
  else if(!strcmp(cmdname,"GetResiMinRMSD")){
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
    EVOEF_GetResiMinRmsdRotFromLab(&structure,PDBID);
  }
  else{
    printf("unknown command name: %s\n, EVOEF exits\n", cmdname);
    exit(ValueError);
  }

  StructureDestroy(&structure);
  ResiTopoSetDestroy(&resiTopo);
  AtomParamsSetDestroy(&atomParam);

  //output the energy weights
  EnergyWeightWrite(FILE_WEIGHT_WRITE);

  clock_t timeend = clock();
  SpentTimeShow(timestart,timeend);
  
  return Success;
}



BOOL CheckCommandName(char* queryname){
  int MAX_CMD_NUM = 100;
  char *supportedcmd[] = {
    "RepairStructure", 
    "ComputeStability", 
    "ComputeBinding", 
    "BuildMutant",
    "ComputeResiEnergy",
    "ComputeRotamersEnergy",
    "ComputeWtRotamersEnergy",
    "OptimizeHydrogen",
    "AddHydrogen",
    "FindInterfaceResidue",
    "ShowResiComposition",
    "MakeLigandEnsemble",
    "ScreenLigandEnsemble",
    "MakeSelfEnergyMatrix",
    "ComputeEvolutionScore",
    "EnzymeDesign",
    "ProteinDesign",
    "EvoDesign",
    "ProtLigDesign",
    "OptimizeWeight",
    "GetPhiPsi",
    "TestBBdepRotLib",
    "CheckResiInLab",
    "SidechainRepack",
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

int GetAttachedAtomTypeNum(AtomArray* pAtoms,int bonds[][2],int bondNum,int atomIndex,int attachAtomTypeNum[10],int *attachCount,int attachIndex[10]){
  //attach[0]:C
  //attach[1]:H
  //attach[2]:O
  //attach[3]:N
  //attach[4]:P
  //attach[5]:S
  //attach[6];F
  //attach[7]:Cl
  //attach[8]:Br
  //attach[9]:I
  //attach[10]:B
  //attach[11]:Fe/Zn/X,metal
  for(int i=0;i<bondNum;i++){
    if(atomIndex==bonds[i][0]-1){
      attachIndex[*attachCount]=bonds[i][1]-1;
      (*attachCount)++;
      Atom* pAtom = AtomArrayGet(pAtoms,bonds[i][1]-1);
      if(pAtom->name[0]=='C') attachAtomTypeNum[0]++;
      else if(pAtom->name[0]=='H')attachAtomTypeNum[1]++;
      else if(pAtom->name[0]=='O')attachAtomTypeNum[2]++;
      else if(pAtom->name[0]=='N')attachAtomTypeNum[3]++;
      else attachAtomTypeNum[4]++;
    }
    else if(atomIndex==bonds[i][1]-1){
      attachIndex[*attachCount]=bonds[i][0]-1;
      (*attachCount)++;
      Atom* pAtom = AtomArrayGet(pAtoms,bonds[i][0]-1);
      if(pAtom->name[0]=='C') attachAtomTypeNum[0]++;
      else if(pAtom->name[0]=='H')attachAtomTypeNum[1]++;
      else if(pAtom->name[0]=='O')attachAtomTypeNum[2]++;
      else if(pAtom->name[0]=='N')attachAtomTypeNum[3]++;
      else attachAtomTypeNum[4]++;
    }
  }

  return Success;
}

int GenerateParameterAndTopologyFromMol2(char* mol2file,char* parfile,char* topfile){
  FILE* fpmol2=fopen(mol2file,"r");
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  BOOL readingAtom = FALSE;
  BOOL readingBond = FALSE;

  char resiName[MAX_LENGTH_RESIDUE_NAME+1]="#";

  int bondNum=0;
  int bonds[500][2];

  AtomArray atoms;
  Atom atom;
  AtomArrayCreate(&atoms);
  AtomCreate(&atom);
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
      int atomId;
      char mol2Resi[MAX_LENGTH_RESIDUE_NAME+1];
      sscanf(line,"%d %s %lf %lf %lf %s %d %s %lf",&atomId,atom.name,&atom.xyz.X,&atom.xyz.Y,&atom.xyz.Z,atom.type,&atom.posInChain,mol2Resi,&atom.charge);
      if(resiName[0]=='#'){
        strcpy(resiName,mol2Resi);
        if(strlen(resiName)>3) resiName[3]='\0';
        if(LigandResidueNameConflictWithAminoAcid(resiName)){
          strcpy(resiName,"LIG");
        }
      }
      AtomArrayAppend(&atoms,&atom);
    }

    if(readingBond){
      int id;
      char type[5];
      sscanf(line,"%d %d %d %s",&id,&bonds[bondNum][0],&bonds[bondNum][1],type);
      bondNum++;
    }
  }
  AtomDestroy(&atom);

  //create atom parameter for ligand
  for(int i=0; i<AtomArrayGetCount(&atoms); i++){
    int attachAtomTypeNum[10];
    int attachCount=0;
    int attachIndex[10];
    for(int j=0;j<10;j++) attachAtomTypeNum[j]=0;
    for(int j=0;j<10;j++) attachIndex[j]=-1;
    Atom* pAtom=AtomArrayGet(&atoms,i);
    pAtom->isBBAtom=FALSE;
    strcpy(pAtom->hbB2,"-");
    strcpy(pAtom->hbDorB,"-");
    strcpy(pAtom->hbHorA,"-");
    pAtom->hybridType=Type_AtomHybridType_None;
    pAtom->polarity=Type_AtomPolarity_NPAliphatic;

    if(pAtom->name[0]=='C' && strstr(pAtom->type,"C.")!=NULL){
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      if(attachAtomTypeNum[1]>=3){
        pAtom->EEF1_atType=EEF1_AtomType_CH3E;
      }
      else if(attachAtomTypeNum[1]==2){
        pAtom->EEF1_atType=EEF1_AtomType_CH2E;
      }
      else if(attachAtomTypeNum[1]==1){
        pAtom->EEF1_atType=EEF1_AtomType_CH1E;
        if(strcmp(pAtom->type,"C.ar")==0){
          pAtom->EEF1_atType=EEF1_AtomType_CR1E;
        }
      }
      else{
        if(strcmp(pAtom->type,"C.ar")==0){
          pAtom->EEF1_atType=EEF1_AtomType_CR;
        }
        else if(strcmp(pAtom->type,"C.cat")==0){
          pAtom->EEF1_atType=EEF1_AtomType_Carg;
        }
        else if(strcmp(pAtom->type,"C.2")==0){
          pAtom->EEF1_atType=EEF1_AtomType_CR;
          if(attachAtomTypeNum[2]>0){
            pAtom->EEF1_atType=EEF1_AtomType_CO;
          }
          for(int k=0;k<attachCount;k++){
            Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
            if(strcmp(pAttachAtom->type,"O.co2")==0 || pAttachAtom->EEF1_atType==EEF1_AtomType_OOC){
              pAtom->EEF1_atType=EEF1_AtomType_COO;
              break;
            }
          }
        }
        else if(strcmp(pAtom->type,"C.3")==0){
          if(attachCount==4){
            pAtom->EEF1_atType=EEF1_AtomType_CH1E;
          }
          else if(attachCount==3){
            pAtom->EEF1_atType=EEF1_AtomType_CH1E;
          }
          else if(attachCount==2){
            pAtom->EEF1_atType=EEF1_AtomType_CH2E;
          }
          else if(attachCount==1){
            pAtom->EEF1_atType=EEF1_AtomType_CH3E;
          }
        }
        else{
          pAtom->EEF1_atType=EEF1_AtomType_CR;
        }
      }
    }

    else if(pAtom->name[0]=='H'){
      pAtom->EEF1_atType=EEF1_AtomType_Hnpl;
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      if(attachAtomTypeNum[2]>0 || attachAtomTypeNum[3]>0){
        pAtom->EEF1_atType=EEF1_AtomType_Hpol;
        pAtom->polarity=Type_AtomPolarity_P;
        strcpy(pAtom->hbHorA,"H");
        strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
      }
    }

    else if(pAtom->name[0]=='O'){
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      pAtom->polarity=Type_AtomPolarity_P;
      if(strcmp(pAtom->type,"O.co2")==0){
        pAtom->EEF1_atType=EEF1_AtomType_OOC;
        pAtom->polarity=Type_AtomPolarity_C;
        strcpy(pAtom->hbHorA,"A");
        pAtom->hybridType=Type_AtomHybridType_SP2;
        for(int k=0;k<attachCount;k++){
          Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
          if(pAttachAtom->name[0]!='H'){
            strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
            break;
          }
        }
      }//COO-
      else if(strcmp(pAtom->type,"O.2")==0){
        pAtom->EEF1_atType=EEF1_AtomType_OC;
        strcpy(pAtom->hbHorA,"A");
        pAtom->hybridType=Type_AtomHybridType_SP2;
        if(attachAtomTypeNum[1]==1){
          pAtom->EEF1_atType=EEF1_AtomType_OH1;
        }
        for(int k=0;k<attachCount;k++){
          Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
          if(pAttachAtom->name[0]!='H'){
            strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
            break;
          }
        }
      }
      else if(strcmp(pAtom->type,"O.3")==0){
        pAtom->EEF1_atType=EEF1_AtomType_Oest;
        if(attachAtomTypeNum[1]==1){
          pAtom->EEF1_atType=EEF1_AtomType_OH1;
          strcpy(pAtom->hbHorA,"A");
          pAtom->hybridType=Type_AtomHybridType_SP3;
          for(int k=0;k<attachCount;k++){
            Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
            if(pAttachAtom->name[0]!='H'){
              strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
              break;
            }
          }
        }
        else if(attachAtomTypeNum[1]==2){
          pAtom->EEF1_atType=EEF1_AtomType_OH2;
          strcpy(pAtom->hbHorA,"A");
          pAtom->hybridType=Type_AtomHybridType_SP3;
          for(int k=0;k<attachCount;k++){
            Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
            if(pAttachAtom->name[0]!='H'){
              strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
              break;
            }
          }
        }
      }
      else{
        pAtom->EEF1_atType=EEF1_AtomType_Oest;
      }
    }

    else if(pAtom->name[0]=='N' && strstr(pAtom->type,"N.")!=NULL){
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      pAtom->polarity=Type_AtomPolarity_P;
      if(strcmp(pAtom->type,"N.4")==0) pAtom->polarity=Type_AtomPolarity_C;
      if(attachAtomTypeNum[1]==3){
        pAtom->EEF1_atType=EEF1_AtomType_NH3;
        pAtom->hybridType=Type_AtomHybridType_SP3;
      }
      else if(attachAtomTypeNum[1]==2){
        pAtom->EEF1_atType=EEF1_AtomType_NH2;
        pAtom->hybridType=Type_AtomHybridType_SP3;
        if(strcmp(pAtom->type,"N.2")==0) pAtom->hybridType=Type_AtomHybridType_SP2;
        for(int k=0;k<attachCount;k++){
          Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
          if(strcmp(pAttachAtom->type,"C.cat")==0 || pAttachAtom->EEF1_atType==EEF1_AtomType_Carg){
            pAtom->EEF1_atType=EEF1_AtomType_Narg;
            break;
          }
        }
      }
      else if(attachAtomTypeNum[1]==1){
        pAtom->EEF1_atType=EEF1_AtomType_NH1;
        pAtom->hybridType=Type_AtomHybridType_SP3;
        if(strcmp(pAtom->type,"N.2")==0) pAtom->hybridType=Type_AtomHybridType_SP2;
        for(int k=0;k<attachCount;k++){
          Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
          if(strcmp(pAttachAtom->type,"C.cat")==0 || pAttachAtom->EEF1_atType==EEF1_AtomType_Carg){
            pAtom->EEF1_atType=EEF1_AtomType_Narg;
            break;
          }
        }
      }
      else{
        if(strcmp(pAtom->type,"N.ar")==0 || strcmp(pAtom->type,"N.2")==0){
          pAtom->EEF1_atType=EEF1_AtomType_NR;
          if(attachCount<3){
            strcpy(pAtom->hbHorA,"A");
            pAtom->hybridType=Type_AtomHybridType_SP2;
            for(int k=0;k<attachCount;k++){
              Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
              if(pAttachAtom->name[0]!='H'){
                strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
                break;
              }
            }
          }
          else{
            pAtom->EEF1_atType=EEF1_AtomType_Npro;
          }
        }
        else if(strcmp(pAtom->type,"N.am")==0){
          pAtom->EEF1_atType=EEF1_AtomType_Npro;
        }
        else if(strcmp(pAtom->type,"N.3")==0 || strcmp(pAtom->type,"N.pl3")==0){
          pAtom->EEF1_atType=EEF1_AtomType_Npro;
        }
        else if(attachCount==3){
          pAtom->EEF1_atType=EEF1_AtomType_Npro;
        }
        else if(strcmp(pAtom->type,"N.1")==0 && attachCount==1){
          pAtom->EEF1_atType=EEF1_AtomType_NR;
        }
        else{
          pAtom->EEF1_atType=EEF1_AtomType_NR;
        }
      }
    }

    else if(pAtom->name[0]=='P'){
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      pAtom->EEF1_atType=EEF1_AtomType_P;
    }

    else if(pAtom->name[0]=='S'){
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      if(attachAtomTypeNum[1]==1){
        pAtom->EEF1_atType=EEF1_AtomType_SH1E;
      }
      else{
        pAtom->EEF1_atType=EEF1_AtomType_S;
      }
    }

    else if(pAtom->name[0]=='B' && strcmp(pAtom->type,"B")==0){
      pAtom->EEF1_atType=EEF1_AtomType_B;
    }
    else if(pAtom->name[0]=='F' && strcmp(pAtom->type,"F")==0){
      pAtom->EEF1_atType=EEF1_AtomType_F;
    }
    else if(pAtom->name[0]=='C' && strcmp(pAtom->type,"Cl")==0){
      pAtom->EEF1_atType=EEF1_AtomType_Cl;
    }
    else if(pAtom->name[0]=='B' && strcmp(pAtom->type,"Br")==0){
      pAtom->EEF1_atType=EEF1_AtomType_Br;
    }
    else if(pAtom->name[0]=='I' && strcmp(pAtom->type,"I")==0){
      pAtom->EEF1_atType=EEF1_AtomType_I;
    }
    else if(pAtom->name[0]=='F' && strcmp(pAtom->type,"Fe")==0){
      pAtom->EEF1_atType=EEF1_AtomType_Fe;
    }
    else if(pAtom->name[0]=='Z' && strcmp(pAtom->type,"F")==0){
      pAtom->EEF1_atType=EEF1_AtomType_Zn;
    }
    else if(pAtom->name[0]=='A' && strcmp(pAtom->type,"Al")==0){
      pAtom->EEF1_atType=EEF1_AtomType_Al;
    }
    else if(pAtom->name[0]=='N' && strcmp(pAtom->type,"Na")==0){
      pAtom->EEF1_atType=EEF1_AtomType_Na;
    }
    else if(pAtom->name[0]=='K' && strcmp(pAtom->type,"K")==0){
      pAtom->EEF1_atType=EEF1_AtomType_K;
    }
    else if(pAtom->name[0]=='C' && strcmp(pAtom->type,"Ca")==0){
      pAtom->EEF1_atType=EEF1_AtomType_Ca;
    }
    else if(pAtom->name[0]=='M' && strcmp(pAtom->type,"Mg")==0){
      pAtom->EEF1_atType=EEF1_AtomType_Mg;
    }
    else{
      pAtom->EEF1_atType=EEF1_AtomType_Other;
    }

    AssignAtomParameterByEEF1Type(pAtom,pAtom->EEF1_atType);
  }

  //output the parameter file
  FILE* fpar=fopen(parfile,"w");
  fprintf(fpar,"!%s\n",resiName);
  fprintf(fpar,"!Resi   Name   Type   Backbone   Polar   epsilon   Rmin     charge  HbH/A  HbD/B  HbB2  Hybrid  DG_free  Volume  Lambda\n");
  for(int i=0; i<AtomArrayGetCount(&atoms); i++){
    Atom* pAtom=AtomArrayGet(&atoms,i);
    char polarity[10];
    char hybrid[5];
    if(pAtom->polarity==Type_AtomPolarity_P) strcpy(polarity,"P");
    else if(pAtom->polarity==Type_AtomPolarity_C) strcpy(polarity,"C");
    else if(pAtom->polarity==Type_AtomPolarity_NPAliphatic) strcpy(polarity,"NP1");
    else strcpy(polarity,"NP2");
    if(pAtom->hybridType==Type_AtomHybridType_SP) strcpy(hybrid,"SP1");
    else if(pAtom->hybridType==Type_AtomHybridType_SP2) strcpy(hybrid,"SP2");
    else if(pAtom->hybridType==Type_AtomHybridType_SP3) strcpy(hybrid,"SP3");
    else strcpy(hybrid,"-");
    fprintf(fpar,"%-6s  %-6s %-6s %-10s %-7s %-8.4f %7.4f  %5.2f    %-5s  %-5s  %-5s %-4s   %7.2f     %4.1f    %3.1f\n",
      resiName,pAtom->name,pAtom->type,"N",polarity,
      pAtom->vdw_epsilon,pAtom->vdw_radius,pAtom->charge,
      pAtom->hbHorA,pAtom->hbDorB,pAtom->hbB2,hybrid,
      pAtom->EEF1_freeDG,pAtom->EEF1_volume,pAtom->EEF1_lamda_);
  }
  fclose(fpar);

  //output the topology file
  FILE* ftop=fopen(topfile,"w");
  fprintf(ftop,"RESI %s\n",resiName);
  for(int i=0;i<AtomArrayGetCount(&atoms);i++){
    fprintf(ftop,"ATOM %s\n",AtomGetName(AtomArrayGet(&atoms,i)));
  }
  fprintf(ftop,"\n");
  for(int i=0;i<bondNum;i++){
    int from=bonds[i][0];
    int to=bonds[i][1];
    fprintf(ftop,"BOND %s %s\n",AtomGetName(AtomArrayGet(&atoms,from-1)),AtomGetName(AtomArrayGet(&atoms,to-1)));
  }
  fprintf(ftop,"\n");
  //calculate charmmIC for hydrogen atoms
  for(int i=0;i<AtomArrayGetCount(&atoms);i++){
    Atom* pAtom=AtomArrayGet(&atoms,i);
    if(AtomIsHydrogen(pAtom)==FALSE) continue;
    //get atomC
    int atomcIndex=-1;
    for(int j=0;j<bondNum;j++){
      if(bonds[j][0]-1==i){atomcIndex=bonds[j][1]-1; break;}
      else if(bonds[j][1]-1==i){atomcIndex=bonds[j][0]-1; break;}
    }
    //get atomB's
    int atomBnum=0;
    int atomBindex[3]={-1,-1,-1};
    for(int j=0;j<bondNum;j++){
      if(bonds[j][0]-1==atomcIndex && bonds[j][1]-1!=i && AtomIsHydrogen(AtomArrayGet(&atoms,bonds[j][1]-1))==FALSE){atomBindex[atomBnum++]=bonds[j][1]-1;}
      else if(bonds[j][1]-1==atomcIndex && bonds[j][0]-1!=i && AtomIsHydrogen(AtomArrayGet(&atoms,bonds[j][0]-1))==FALSE){atomBindex[atomBnum++]=bonds[j][0]-1;}
    }
    //check if improper dihedral angle can be used
    if(atomBnum>=2){
      Atom* pAtomA=AtomArrayGet(&atoms,atomBindex[0]);
      Atom* pAtomB=AtomArrayGet(&atoms,atomBindex[1]);
      Atom* pAtomC=AtomArrayGet(&atoms,atomcIndex);
      Atom* pAtomD=AtomArrayGet(&atoms,i);
      XYZ AB=XYZDifference(&pAtomA->xyz,&pAtomB->xyz);
      XYZ BC=XYZDifference(&pAtomB->xyz,&pAtomC->xyz);
      XYZ CD = XYZDifference(&pAtomC->xyz,&pAtomD->xyz);
      XYZ AC = XYZDifference(&pAtomA->xyz,&pAtomC->xyz);

      double dAB=XYZDistance(&pAtomA->xyz,&pAtomB->xyz);
      double aABC=XYZAngle(&AC,&BC);
      double xABCD=GetTorsionAngle(&pAtomA->xyz,&pAtomB->xyz,&pAtomC->xyz,&pAtomD->xyz);
      double aBCD=PI-XYZAngle(&BC,&CD);
      double dCD=XYZDistance(&pAtomC->xyz,&pAtomD->xyz);
      fprintf(ftop,"IC   %-7.7s %-7.7s *%-7.7s %-7.7s  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
        AtomGetName(pAtomA),AtomGetName(pAtomB),AtomGetName(pAtomC),AtomGetName(pAtomD),dAB,RadToDeg(aABC),RadToDeg(xABCD),RadToDeg(aBCD),dCD);
    }
    else{
      //get atomA's
      int atomaIndex=-1;
      for(int j=0;j<bondNum;j++){
        if(bonds[j][0]-1==atomBindex[0] && bonds[j][1]-1!=atomcIndex){
          atomaIndex=bonds[j][1]-1;
          break;
        }
        else if(bonds[j][1]-1==atomBindex[0] && bonds[j][0]-1!=atomcIndex){
          atomaIndex=bonds[j][0]-1;
          break;
        }
      }
      Atom* pAtomA=AtomArrayGet(&atoms,atomaIndex);
      Atom* pAtomB=AtomArrayGet(&atoms,atomBindex[0]);
      Atom* pAtomC=AtomArrayGet(&atoms,atomcIndex);
      Atom* pAtomD=AtomArrayGet(&atoms,i);
      XYZ AB=XYZDifference(&pAtomA->xyz,&pAtomB->xyz);
      XYZ BC=XYZDifference(&pAtomB->xyz,&pAtomC->xyz);
      XYZ CD = XYZDifference(&pAtomC->xyz,&pAtomD->xyz);
      XYZ AC = XYZDifference(&pAtomA->xyz,&pAtomC->xyz);

      double dAB=XYZDistance(&pAtomA->xyz,&pAtomB->xyz);
      double aABC=PI-XYZAngle(&AB,&BC);
      double xABCD=GetTorsionAngle(&pAtomA->xyz,&pAtomB->xyz,&pAtomC->xyz,&pAtomD->xyz);
      double aBCD=PI-XYZAngle(&BC,&CD);
      double dCD=XYZDistance(&pAtomC->xyz,&pAtomD->xyz);
      fprintf(ftop,"IC   %-7.7s %-7.7s  %-7.7s %-7.7s  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
        AtomGetName(pAtomA),AtomGetName(pAtomB),AtomGetName(pAtomC),AtomGetName(pAtomD),dAB,RadToDeg(aABC),RadToDeg(xABCD),RadToDeg(aBCD),dCD);
    }
  }
  fclose(ftop);

  AtomArrayDestroy(&atoms);

  return Success;
}

