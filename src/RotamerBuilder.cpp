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

#include "RotamerBuilder.h"
#include "EnergyFunction.h"
#include <string.h>

extern double PPI_DIST_CUTOFF;
extern BOOL FLAG_ADD_CRYSTAL_ROT;
extern BOOL FLAG_EXPAND_HYDROXYL_ROT;



#define DEAL_WITH_PROTEIN_ROTAMERS_BBIND
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// functions used to deal with sidechain rotamers, sidechain repacking and protein design
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ProteinSiteBuildAllRotamers(Structure* pThis,int chainIndex,int resiIndex,BBindRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos){
  int result = Success;
  DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if(pCurrentDesignSite != NULL){
    DesignSiteRemoveRotamers(pCurrentDesignSite);
  }
  else{
    (pThis->designSiteCount)++;
    pThis->designSites=(DesignSite*)realloc(pThis->designSites, sizeof(DesignSite)*pThis->designSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->designSiteCount-1]);
    pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
    Chain* pDestChain = StructureGetChain(pThis, chainIndex);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
    pCurrentDesignSite->pResidue = pDestResidue;
    pCurrentDesignSite->chainIndex = chainIndex;
    pCurrentDesignSite->resiIndex = resiIndex;
  }

  // set design types - 20 AA types;
  StringArray designTypes;
  StringArray patchTypes;
  StringArrayCreate(&designTypes);
  StringArrayCreate(&patchTypes);
  StringArrayAppend(&designTypes, "ALA"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ARG"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ASN"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ASP"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "CYS"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "GLN"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "GLU"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "GLY"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "HSD"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "HSE"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ILE"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "LEU"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "LYS"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "MET"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "PHE"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "PRO"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "SER"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "THR"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "TRP"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "TYR"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "VAL"); StringArrayAppend(&patchTypes, "");
  result = RotamerSetOfProteinGenerate(DesignSiteGetRotamers(pCurrentDesignSite),pCurrentDesignSite->pResidue,&designTypes,&patchTypes,rotlib,atomParams,resiTopos);
  //ResidueSetDesignSiteFlag(pCurrentDesignSite->pResidue, Type_ResidueDesignType_Mutated);
  StringArrayDestroy(&patchTypes);
  StringArrayDestroy(&designTypes);

  return result;
}


int ProteinSiteBuildSpecifiedRotamers(Structure* pThis, int chainIndex, int resiIndex, BBindRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos, Type_ResidueDesignType type){
  int result = Success;
  DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if(pCurrentDesignSite != NULL){
    DesignSiteRemoveRotamers(pCurrentDesignSite);
  }
  else{
    (pThis->designSiteCount)++;
    pThis->designSites=(DesignSite*)realloc(pThis->designSites, sizeof(DesignSite)*pThis->designSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->designSiteCount-1]);
    pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
    Chain* pDestChain = StructureGetChain(pThis, chainIndex);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
    pCurrentDesignSite->pResidue = pDestResidue;
    pCurrentDesignSite->chainIndex = chainIndex;
    pCurrentDesignSite->resiIndex = resiIndex;
  }
  // set design types - 20 AA types;
  StringArray designTypes;
  StringArray patchTypes;
  StringArrayCreate(&designTypes);
  StringArrayCreate(&patchTypes);
  if(type == Type_ResidueDesignType_Mutated){
    StringArrayAppend(&designTypes, "ALA"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "ARG"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "ASN"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "ASP"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "CYS"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "GLN"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "GLU"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "GLY"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "HSD"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "HSE"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "ILE"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "LEU"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "LYS"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "MET"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "PHE"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "PRO"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "SER"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "THR"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "TRP"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "TYR"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "VAL"); StringArrayAppend(&patchTypes, "");
  }
  else if(type == Type_ResidueDesignType_Rotameric || type == Type_ResidueDesignType_Catalytic){
    StringArrayAppend(&designTypes, ResidueGetName(pCurrentDesignSite->pResidue)); StringArrayAppend(&patchTypes, "");
  }
  result = RotamerSetOfProteinGenerate(DesignSiteGetRotamers(pCurrentDesignSite),pCurrentDesignSite->pResidue,&designTypes,&patchTypes,rotlib,atomParams,resiTopos);
  StringArrayDestroy(&patchTypes);
  StringArrayDestroy(&designTypes);

  return result;
}




int ProteinSiteBuildMutatedRotamers(Structure* pThis, int chainIndex, int resiIndex, BBindRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos, StringArray *pDesignTypes, StringArray *pPatchTypes){
  int result = Success;
  DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if(pCurrentDesignSite != NULL){
    DesignSiteRemoveRotamers(pCurrentDesignSite);
  }
  else{
    (pThis->designSiteCount)++;
    pThis->designSites=(DesignSite*)realloc(pThis->designSites, sizeof(DesignSite)*pThis->designSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->designSiteCount-1]);
    pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
    Chain* pDestChain = StructureGetChain(pThis, chainIndex);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
    pCurrentDesignSite->pResidue = pDestResidue;
    pCurrentDesignSite->chainIndex = chainIndex;
    pCurrentDesignSite->resiIndex = resiIndex;
  }

  result = RotamerSetOfProteinGenerate(DesignSiteGetRotamers(pCurrentDesignSite),pCurrentDesignSite->pResidue,pDesignTypes,pPatchTypes,rotlib,atomParams,resiTopos);
  ResidueSetDesignSiteFlag(pCurrentDesignSite->pResidue,Type_ResidueDesignType_Mutated);
  return result;
}

int ProteinSiteBuildWildtypeRotamers(Structure* pThis, int chainIndex, int resiIndex, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos){
  int result = Success;
  DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if(pCurrentDesignSite != NULL){
    DesignSiteRemoveRotamers(pCurrentDesignSite);
  }
  else{
    (pThis->designSiteCount)++;
    pThis->designSites=(DesignSite*)realloc(pThis->designSites, sizeof(DesignSite)*pThis->designSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->designSiteCount-1]);
    pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
    Chain* pDestChain = StructureGetChain(pThis, chainIndex);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
    pCurrentDesignSite->pResidue = pDestResidue;
    pCurrentDesignSite->chainIndex = chainIndex;
    pCurrentDesignSite->resiIndex = resiIndex;
  }

  // set native design types;
  StringArray designTypes;
  StringArray patchTypes;
  StringArrayCreate(&designTypes);
  StringArrayCreate(&patchTypes);
  StringArrayAppend(&designTypes, ResidueGetName(pCurrentDesignSite->pResidue));
  StringArrayAppend(&patchTypes, "");
  if(strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "HSD") == 0){
    StringArrayAppend(&designTypes, "HSE");
    StringArrayAppend(&patchTypes, "");
  }
  else if(strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "HSE") == 0){
    StringArrayAppend(&designTypes, "HSD");
    StringArrayAppend(&patchTypes, "");
  }

  result = RotamerSetOfProteinGenerate(DesignSiteGetRotamers(pCurrentDesignSite),pCurrentDesignSite->pResidue,&designTypes,&patchTypes,rotlib,atomParams,resiTopos);
  ResidueSetDesignSiteFlag(pCurrentDesignSite->pResidue, Type_ResidueDesignType_Rotameric);

  StringArrayDestroy(&patchTypes);
  StringArrayDestroy(&designTypes);

  return result;
}

// this function can be used to build crystal rotamers for every amino acid type
int ProteinSiteAddCrystalRotamer(Structure* pThis, int chainIndex, int resiIndex, ResiTopoSet *pResiTopos){
  Chain* pDestChain = StructureGetChain(pThis, chainIndex);
  Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
  if(pDestChain->type == Type_Chain_Protein){
    DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
    if(pCurrentDesignSite == NULL){
      (pThis->designSiteCount)++;
      pThis->designSites=(DesignSite*)realloc(pThis->designSites, sizeof(DesignSite)*pThis->designSiteCount);
      DesignSiteCreate(&pThis->designSites[pThis->designSiteCount-1]);
      pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
      pCurrentDesignSite->pResidue = pDestResidue;
      pCurrentDesignSite->chainIndex = chainIndex;
      pCurrentDesignSite->resiIndex = resiIndex;
    }

    //do not add rotamer for residue ala and gly
    if(strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "ALA") == 0 || strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "GLY") == 0){
      return Success;
    }
    //else add a crystal rotamer for the residue
    RotamerSet* pSetI = DesignSiteGetRotamers(pCurrentDesignSite);
    Rotamer tempRotamer;
    Rotamer* pRotamerRepresentative;
    RotamerCreate(&tempRotamer);
    RotamerSetType(&tempRotamer,ResidueGetName(pCurrentDesignSite->pResidue));
    RotamerSetChainName(&tempRotamer,ResidueGetChainName(pCurrentDesignSite->pResidue));
    RotamerSetPosInChain(&tempRotamer,ResidueGetPosInChain(pCurrentDesignSite->pResidue));
    pRotamerRepresentative = RotamerSetGetRepresentative(pSetI, RotamerGetType(&tempRotamer));
    if(pRotamerRepresentative != NULL){
      AtomArrayCopy(&tempRotamer.atoms, &pRotamerRepresentative->atoms);
      for(int i = 0; i < RotamerGetAtomCount(&tempRotamer); i++){
        Atom* pAtom = RotamerGetAtom(&tempRotamer, i);
        pAtom->xyz = ResidueGetAtomByName(pCurrentDesignSite->pResidue, AtomGetName(pAtom))->xyz;
      }
      BondSetCopy(&tempRotamer.bonds,&pRotamerRepresentative->bonds);
      XYZArrayResize(&tempRotamer.xyzs, RotamerGetAtomCount(pRotamerRepresentative));
      for(int i = 0; i < XYZArrayGetLength(&tempRotamer.xyzs); i++){
        XYZArraySet(&tempRotamer.xyzs, i, &AtomArrayGet(&tempRotamer.atoms,i)->xyz);
      }
      DoubleArrayCopy(&tempRotamer.xtorsions,&pDestResidue->xtorsions);
      RotamerSetAdd(&pCurrentDesignSite->rotamers, &tempRotamer);
    }
    else{
      AtomArrayCopy(&tempRotamer.atoms, &pDestResidue->atoms);
      BondSetCopy(&tempRotamer.bonds,&pDestResidue->bonds);
      XYZArrayResize(&tempRotamer.xyzs, AtomArrayGetCount(&tempRotamer.atoms));
      for(int i = 0; i < XYZArrayGetLength(&tempRotamer.xyzs); i++){
        XYZArraySet(&tempRotamer.xyzs, i, &AtomArrayGet(&tempRotamer.atoms,i)->xyz);
      }
      DoubleArrayCopy(&tempRotamer.xtorsions,&pDestResidue->xtorsions);
      RotamerSetAdd(&pCurrentDesignSite->rotamers, &tempRotamer);
    }
    RotamerDestroy(&tempRotamer);
    //if residue is histidine, add a flipped rotamer
    if(strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "HSD") == 0){
      if((pRotamerRepresentative = RotamerSetGetRepresentative(pSetI, "HSE")) != NULL){
        Residue newResi;
        ResidueCreate(&newResi);
        ResidueSetName(&newResi, "HSE");
        AtomArrayCopy(&newResi.atoms, &pRotamerRepresentative->atoms);
        for(int i = 0; i < ResidueGetAtomCount(&newResi); i++){
          Atom* pAtom = ResidueGetAtom(&newResi, i);
          if(pAtom->isBBAtom == FALSE && AtomIsHydrogen(pAtom) == TRUE){
            pAtom->isXyzValid = FALSE;
            continue;
          }
          pAtom->xyz = ResidueGetAtomByName(pCurrentDesignSite->pResidue, AtomGetName(pAtom))->xyz;
        }
        ResidueCalcAllAtomXYZ(&newResi, pResiTopos, NULL, NULL);

        RotamerCreate(&tempRotamer);
        RotamerSetType(&tempRotamer,"HSE");
        RotamerSetChainName(&tempRotamer,ResidueGetChainName(pCurrentDesignSite->pResidue));
        RotamerSetPosInChain(&tempRotamer,ResidueGetPosInChain(pCurrentDesignSite->pResidue));

        AtomArrayCopy(&tempRotamer.atoms, &pRotamerRepresentative->atoms);
        for(int i = 0; i < RotamerGetAtomCount(&tempRotamer); i++){
          Atom* pAtom = RotamerGetAtom(&tempRotamer, i);
          pAtom->xyz = ResidueGetAtomByName(&newResi, AtomGetName(pAtom))->xyz;
        }
        BondSetCopy(&tempRotamer.bonds,&pRotamerRepresentative->bonds);

        XYZArrayResize(&tempRotamer.xyzs, RotamerGetAtomCount(pRotamerRepresentative));
        for(int i = 0; i < XYZArrayGetLength(&tempRotamer.xyzs); i++){
          XYZArraySet(&tempRotamer.xyzs, i, &AtomArrayGet(&tempRotamer.atoms,i)->xyz);
        }
        DoubleArrayCopy(&tempRotamer.xtorsions,&pDestResidue->xtorsions);
        RotamerSetAdd(&pCurrentDesignSite->rotamers, &tempRotamer);

        RotamerDestroy(&tempRotamer);
        ResidueDestroy(&newResi);
      }
    } // HSD
    else if(strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "HSE") == 0){
      if((pRotamerRepresentative = RotamerSetGetRepresentative(pSetI, "HSD")) != NULL){
        Residue newResi;
        ResidueCreate(&newResi);
        ResidueSetName(&newResi, "HSD");
        AtomArrayCopy(&newResi.atoms, &pRotamerRepresentative->atoms);
        for(int i = 0; i < ResidueGetAtomCount(&newResi); i++){
          Atom* pAtom = ResidueGetAtom(&newResi, i);
          if(pAtom->isBBAtom == FALSE && AtomIsHydrogen(pAtom) == TRUE){
            pAtom->isXyzValid = FALSE;
            continue;
          }
          pAtom->xyz = ResidueGetAtomByName(pCurrentDesignSite->pResidue, AtomGetName(pAtom))->xyz;
        }
        ResidueCalcAllAtomXYZ(&newResi, pResiTopos, NULL, NULL);

        RotamerCreate(&tempRotamer);
        RotamerSetType(&tempRotamer,"HSD");
        RotamerSetChainName(&tempRotamer,ResidueGetChainName(pCurrentDesignSite->pResidue));
        RotamerSetPosInChain(&tempRotamer,ResidueGetPosInChain(pCurrentDesignSite->pResidue));

        AtomArrayCopy(&tempRotamer.atoms, &pRotamerRepresentative->atoms);
        for(int i = 0; i < RotamerGetAtomCount(&tempRotamer); i++){
          Atom* pAtom = RotamerGetAtom(&tempRotamer, i);
          pAtom->xyz = ResidueGetAtomByName(&newResi, AtomGetName(pAtom))->xyz;
        }
        BondSetCopy(&tempRotamer.bonds,&pRotamerRepresentative->bonds);

        XYZArrayResize(&tempRotamer.xyzs, RotamerGetAtomCount(pRotamerRepresentative));
        for(int i = 0; i < XYZArrayGetLength(&tempRotamer.xyzs); i++){
          XYZArraySet(&tempRotamer.xyzs, i, &AtomArrayGet(&tempRotamer.atoms,i)->xyz);
        }
        DoubleArrayCopy(&tempRotamer.xtorsions,&pDestResidue->xtorsions);
        RotamerSetAdd(&pCurrentDesignSite->rotamers, &tempRotamer);
        RotamerDestroy(&tempRotamer);
        ResidueDestroy(&newResi);
      }
    } //HSE
  }

  return Success;
}


//this function is used for build flipped rotamers for ASN/GLN/HIS, note that this function is based on the above function
//ProteinSiteBuildCrystalRotamer(), first build crystal rotamer then flip
int ProteinSiteBuildFlippedCrystalRotamer(Structure* pStructure, int chainIndex, int resiIndex, ResiTopoSet *pResiTopos){
  Chain *pDestChain = StructureGetChain(pStructure,chainIndex);
  Residue* pDesignResi = ChainGetResidue(pDestChain, resiIndex);
  if(pDestChain->type == Type_Chain_Protein){
    if(strcmp(ResidueGetName(pDesignResi),"ASN")!=0 && strcmp(ResidueGetName(pDesignResi),"GLN")!=0 &&
      strcmp(ResidueGetName(pDesignResi),"HSD")!=0 && strcmp(ResidueGetName(pDesignResi),"HSE")!=0){
        return Success;
    }

    DesignSite *pCurrenDesignSite = StructureFindDesignSite(pStructure, chainIndex, resiIndex);
    //step2: flip rotamer
    RotamerSet* pCurrentRotamerSet=DesignSiteGetRotamers(pCurrenDesignSite);
    int rotCount=RotamerSetGetCount(pCurrentRotamerSet);
    for(int i=0; i<rotCount; i++){
      Rotamer* pRotamer=RotamerSetGet(pCurrentRotamerSet,i);
      Rotamer* pRepresentative=RotamerSetGetRepresentative(pCurrentRotamerSet,RotamerGetType(pRotamer));
      Rotamer tempRotamer;
      RotamerCreate(&tempRotamer);
      RotamerCopy(&tempRotamer,pRepresentative);
      if(strcmp(RotamerGetType(pRotamer),"ASN")==0){
        int index1,index2;
        RotamerFindAtom(&tempRotamer,"ND2",&index1);
        RotamerFindAtom(&tempRotamer,"OD1",&index2);
        XYZ tempXYZ1 = RotamerGetAtom(&tempRotamer,index1)->xyz;
        XYZ tempXYZ2 = RotamerGetAtom(&tempRotamer,index2)->xyz;
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ2;
        RotamerGetAtom(&tempRotamer,index2)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs,index2,&tempXYZ1);
        ResidueTopology resiTopo;
        CharmmIC ic;
        ResidueTopologyCreate(&resiTopo);
        CharmmICCreate(&ic);
        ResiTopoSetGet(pResiTopos,"ASN",&resiTopo);
        ResidueTopologyFindCharmmIC(&resiTopo,"HD21",&ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer,"HD21",&index1);
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ1);
        ResidueTopologyFindCharmmIC(&resiTopo,"HD22",&ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer,"HD22",&index1);
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ1);
        ResidueTopologyDestroy(&resiTopo);
        CharmmICDestroy(&ic);
      }
      else if(strcmp(RotamerGetType(pRotamer),"GLN")==0){
        int index1,index2;
        RotamerFindAtom(&tempRotamer,"NE2",&index1);
        RotamerFindAtom(&tempRotamer,"OE1",&index2);
        XYZ tempXYZ1 = RotamerGetAtom(&tempRotamer,index1)->xyz;
        XYZ tempXYZ2 = RotamerGetAtom(&tempRotamer,index2)->xyz;
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ2;
        RotamerGetAtom(&tempRotamer,index2)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs,index2,&tempXYZ1);
        ResidueTopology resiTopo;
        CharmmIC ic;
        ResidueTopologyCreate(&resiTopo);
        CharmmICCreate(&ic);
        ResiTopoSetGet(pResiTopos,"GLN",&resiTopo);
        ResidueTopologyFindCharmmIC(&resiTopo,"HE21",&ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer,"HE21",&index1);
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ1);
        ResidueTopologyFindCharmmIC(&resiTopo,"HE22",&ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer,"HE22",&index1);
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ1);
        ResidueTopologyDestroy(&resiTopo);
        CharmmICDestroy(&ic);
      }
      else if(strcmp(RotamerGetType(pRotamer),"HSD")==0){
        int index1,index2;
        RotamerFindAtom(&tempRotamer,"CD2",&index1);
        RotamerFindAtom(&tempRotamer,"ND1",&index2);
        XYZ tempXYZ1 = RotamerGetAtom(&tempRotamer,index1)->xyz;
        XYZ tempXYZ2 = RotamerGetAtom(&tempRotamer,index2)->xyz;
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ2;
        RotamerGetAtom(&tempRotamer,index2)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs,index2,&tempXYZ1);
        RotamerFindAtom(&tempRotamer,"NE2",&index1);
        RotamerFindAtom(&tempRotamer,"CE1",&index2);
        tempXYZ1 = RotamerGetAtom(&tempRotamer,index1)->xyz;
        tempXYZ2 = RotamerGetAtom(&tempRotamer,index2)->xyz;
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ2;
        RotamerGetAtom(&tempRotamer,index2)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs,index2,&tempXYZ1);
        ResidueTopology resiTopo;
        CharmmIC ic;
        ResidueTopologyCreate(&resiTopo);
        CharmmICCreate(&ic);
        ResiTopoSetGet(pResiTopos,"HSD",&resiTopo);
        ResidueTopologyFindCharmmIC(&resiTopo,"HD1",&ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer,"HD1",&index1);
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ1);
        ResidueTopologyDestroy(&resiTopo);
        CharmmICDestroy(&ic);
      }
      else if(strcmp(RotamerGetType(pRotamer),"HSE")==0){
        int index1,index2;
        RotamerFindAtom(&tempRotamer,"CD2",&index1);
        RotamerFindAtom(&tempRotamer,"ND1",&index2);
        XYZ tempXYZ1 = RotamerGetAtom(&tempRotamer,index1)->xyz;
        XYZ tempXYZ2 = RotamerGetAtom(&tempRotamer,index2)->xyz;
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ2;
        RotamerGetAtom(&tempRotamer,index2)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs,index2,&tempXYZ1);
        RotamerFindAtom(&tempRotamer,"NE2",&index1);
        RotamerFindAtom(&tempRotamer,"CE1",&index2);
        tempXYZ1 = RotamerGetAtom(&tempRotamer,index1)->xyz;
        tempXYZ2 = RotamerGetAtom(&tempRotamer,index2)->xyz;
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ2;
        RotamerGetAtom(&tempRotamer,index2)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs,index2,&tempXYZ1);
        ResidueTopology resiTopo;
        CharmmIC ic;
        ResidueTopologyCreate(&resiTopo);
        CharmmICCreate(&ic);
        ResiTopoSetGet(pResiTopos,"HSE",&resiTopo);
        ResidueTopologyFindCharmmIC(&resiTopo,"HE2",&ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer,"HE2",&index1);
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ1);
        ResidueTopologyDestroy(&resiTopo);
        CharmmICDestroy(&ic);
      }
      tempRotamer.dunbrack=NATIVE_ROTAMER_DUNBRACK;
      RotamerSetAdd(pCurrentRotamerSet,&tempRotamer);
      RotamerDestroy(&tempRotamer);
    }
  }

  return Success;
}


int ProteinSiteExpandHydroxylRotamers(Structure *pStructure, int chainIndex, int resiIndex, ResiTopoSet *pTopos){
  Chain* pChain=StructureGetChain(pStructure,chainIndex);
  if(pChain->type==Type_Chain_Protein){
    DesignSite *pDesignSite = StructureFindDesignSite(pStructure, chainIndex,resiIndex);
    RotamerSet *pRotamerSet = DesignSiteGetRotamers(pDesignSite);
    int rotamerCount = RotamerSetGetCount(pRotamerSet); // the rotamer set will be expanded, we need to record the current count
    ResidueTopology tops;
    CharmmIC ics;
    ResidueTopologyCreate(&tops);
    CharmmICCreate(&ics);
    if(RotamerSetGetRepresentative(pRotamerSet, "SER") != NULL ){
      ResiTopoSetGet(pTopos, "SER", &tops);
      ResidueTopologyFindCharmmIC(&tops, "HG", &ics);
      double icParaX = ics.icParam[2];
      int addedCount=0;
      Rotamer tempRot;
      RotamerCreate(&tempRot);
      for(int j = 0; j < rotamerCount; j++){
        Rotamer *pRotamer = RotamerSetGet(pRotamerSet, j);
        if(strcmp(RotamerGetType(pRotamer), "SER") == 0){
          RotamerRestore(pRotamer, pRotamerSet);
          int atomIndex;
          RotamerFindAtom(pRotamer, "HG", &atomIndex);
          RotamerCopy(&tempRot, pRotamer);
          RotamerExtract(pRotamer);
          RotamerSetAdd(pRotamerSet, &tempRot);
          for(int k = 0; k < EXPANDED_ROT_SER; k++){
            ics.icParam[2] = icParaX + 2.0*PI*(k+1)/(EXPANDED_ROT_SER+1);
            if(ics.icParam[2] > PI) ics.icParam[2] -= 2*PI;
            GetFourthAtom(&RotamerGetAtomByName(&tempRot, "CA")->xyz,&RotamerGetAtomByName(&tempRot, "CB")->xyz,&RotamerGetAtomByName(&tempRot, "OG")->xyz,ics.icParam,&RotamerGetAtomByName(&tempRot, "HG")->xyz);
            XYZArraySet(&tempRot.xyzs, atomIndex, &RotamerGetAtomByName(&tempRot, "HG")->xyz);
            // check if the rotamer should be added
            BOOL expandedRotAccepted = TRUE;
            for(int kk = 0; kk < RotamerGetAtomCount(&tempRot); kk++){
              Atom *pAtomK = RotamerGetAtom(&tempRot, kk);
              if(strcmp(AtomGetName(pAtomK), "HG") != 0) continue;
              for(int ss = 0; ss < RotamerGetAtomCount(&tempRot); ss++){
                Atom *pAtomS = RotamerGetAtom(&tempRot, ss);
                if(strcmp(AtomGetName(pAtomK), AtomGetName(pAtomS)) == 0) continue;
                int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtomK), AtomGetName(pAtomS), RotamerGetBonds(&tempRot));
                if(bondType==14||bondType==15){
                  double distance = XYZDistance(&pAtomK->xyz, &pAtomS->xyz);
                  if(distance < (pAtomK->vdw_radius+pAtomS->vdw_radius)*0.75){
                    expandedRotAccepted = FALSE;
                    break;
                  }
                }
              }
              if(expandedRotAccepted == FALSE) break;
            }
            if(expandedRotAccepted == TRUE){
              RotamerSetAdd(pRotamerSet, &tempRot);
              addedCount++;
            }
          }
        }
      }
      //printf("Design site (%2d, %4d): %d SER rotamers expanded\n", chainIndex, ResidueGetPosInChain(pDesignSite->pResidue), addedCount);
      RotamerDestroy(&tempRot);
    }
    // for thr rotamers
    if(RotamerSetGetRepresentative(pRotamerSet, "THR") != NULL ){
      ResiTopoSetGet(pTopos, "THR", &tops);
      ResidueTopologyFindCharmmIC(&tops, "HG1", &ics);
      double icParaX = ics.icParam[2];
      int addedCount=0;
      Rotamer tempRot;
      RotamerCreate(&tempRot);
      for(int j = 0; j < rotamerCount; j++){
        Rotamer *pRotamer = RotamerSetGet(pRotamerSet, j);
        if(strcmp(RotamerGetType(pRotamer), "THR") == 0){
          RotamerRestore(pRotamer, pRotamerSet);
          int atomIndex;
          RotamerFindAtom(pRotamer, "HG1", &atomIndex);
          RotamerCopy(&tempRot, pRotamer);
          RotamerExtract(pRotamer);
          for(int k = 0; k < EXPANDED_ROT_THR; k++){
            BOOL expandedRotAccepted = TRUE;
            ics.icParam[2] = icParaX + 2.0*PI*(k+1)/(EXPANDED_ROT_THR+1);
            if(ics.icParam[2] > PI) ics.icParam[2] -= 2*PI;
            GetFourthAtom(&RotamerGetAtomByName(&tempRot, "CA")->xyz,&RotamerGetAtomByName(&tempRot, "CB")->xyz,&RotamerGetAtomByName(&tempRot, "OG1")->xyz,ics.icParam,&RotamerGetAtomByName(&tempRot, "HG1")->xyz);
            XYZArraySet(&tempRot.xyzs, atomIndex, &RotamerGetAtomByName(&tempRot, "HG1")->xyz);
            //RotamerSetAdd(&tempSet, &tempRot);
            for(int kk = 0; kk < RotamerGetAtomCount(&tempRot); kk++){
              Atom *pAtomK = RotamerGetAtom(&tempRot, kk);
              if(strcmp(AtomGetName(pAtomK), "HG1") != 0) continue;
              for(int ss = 0; ss < RotamerGetAtomCount(&tempRot); ss++){
                Atom *pAtomS = RotamerGetAtom(&tempRot, ss);
                if(strcmp(AtomGetName(pAtomK), AtomGetName(pAtomS)) == 0) continue;
                int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtomK), AtomGetName(pAtomS), RotamerGetBonds(&tempRot));
                if(bondType==14||bondType==15){
                  double distance = XYZDistance(&pAtomK->xyz, &pAtomS->xyz);
                  if(distance < (pAtomK->vdw_radius+pAtomS->vdw_radius)*0.75){
                    expandedRotAccepted = FALSE;
                    break;
                  }
                }
              }
              if(expandedRotAccepted == FALSE) break;
            }
            if(expandedRotAccepted == TRUE){
              RotamerSetAdd(pRotamerSet, &tempRot);
              addedCount++;
            }
          }

        }
      }
      //printf("Design site (%2d, %4d): %d THR rotamers expanded\n", chainIndex, ResidueGetPosInChain(pDesignSite->pResidue), addedCount);
      RotamerDestroy(&tempRot);
    }
    // for tyr rotamers
    if(RotamerSetGetRepresentative(pRotamerSet, "TYR") != NULL ){
      ResiTopoSetGet(pTopos, "TYR", &tops);
      ResidueTopologyFindCharmmIC(&tops, "HH", &ics);
      double icPara_Tyr = ics.icParam[2];
      int addedCount=0;
      Rotamer tempRot;
      RotamerCreate(&tempRot);
      for(int j = 0; j < rotamerCount; j++){
        Rotamer *pRotamer = RotamerSetGet(pRotamerSet, j);
        if(strcmp(RotamerGetType(pRotamer), "TYR") == 0){
          RotamerRestore(pRotamer, pRotamerSet);
          int atomIndex = -1;
          RotamerFindAtom(pRotamer, "HH", &atomIndex);
          RotamerCopy(&tempRot, pRotamer);
          RotamerExtract(pRotamer);
          for(int k = 0; k < EXPANDED_ROT_TYR; k++){
            ics.icParam[2] = icPara_Tyr + 2.0*PI*(k+1)/(EXPANDED_ROT_TYR+1);
            if(ics.icParam[2] > PI) ics.icParam[2] -= 2*PI;
            GetFourthAtom(&RotamerGetAtomByName(&tempRot, "CE1")->xyz,&RotamerGetAtomByName(&tempRot, "CZ")->xyz,&RotamerGetAtomByName(&tempRot, "OH")->xyz,ics.icParam,&RotamerGetAtomByName(&tempRot, "HH")->xyz);
            XYZArraySet(&tempRot.xyzs, atomIndex, &RotamerGetAtomByName(&tempRot, "HH")->xyz);
            RotamerSetAdd(pRotamerSet, &tempRot);
            addedCount++;
          }
        }

      }
      //printf("Design site (%2d, %4d): %d TYR rotamers expanded\n", chainIndex, ResidueGetPosInChain(pDesignSite->pResidue), addedCount);
      RotamerDestroy(&tempRot);
    }
    ResidueTopologyDestroy(&tops);
    CharmmICDestroy(&ics);
  }
  return Success;
}

int ProteinSiteWriteRotamers(Structure *pStructure, int chainIndex, int resiIndex, char *rotamerFilePath){
  DesignSite *pCurrentDesignSite = StructureFindDesignSite(pStructure, chainIndex, resiIndex);
  RotamerSet *pCurrentRotamerSet  = DesignSiteGetRotamers(pCurrentDesignSite);
  FILE *pOut= fopen(rotamerFilePath, "w");
  for(int i = 0; i < RotamerSetGetCount(pCurrentRotamerSet); i++){
    Rotamer *pRotamer = RotamerSetGet(pCurrentRotamerSet, i);
    RotamerRestore(pRotamer,pCurrentRotamerSet);
    Model(i, pOut);
    RotamerShowInPDBFormat(pRotamer, "ATOM", RotamerGetChainName(pRotamer),1, i, FALSE, pOut);
    EndModel(pOut);
    RotamerExtract(pRotamer);
  }
  fclose(pOut);

  return 0;
}


int StructureGenerateSpecifiedProteinRotamers(Structure* pThis,BBindRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* designsitefile){
  FileReader fr;
  FileReaderCreate(&fr, designsitefile);
  char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
    char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    ExtractFirstStringFromSourceString(keyword, buffer);
    if(strcmp(keyword, "catalytic_sites_begin") == 0){
      StringArray strings;
      StringArrayCreate(&strings);
      while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
        StringArraySplitString(&strings, buffer, ' ');
        if(strcmp(StringArrayGet(&strings, 0), "catalytic_sites_end") == 0){
          break;
        }
        char* chnname = StringArrayGet(&strings, 0);
        int posInChain = atoi(StringArrayGet(&strings, 1));
        Chain* pChain = StructureFindChainByName(pThis, chnname);
        int resiIndex = -1;
        ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
        if(resiIndex != -1){
          ChainGetResidue(pChain, resiIndex)->designSiteType = Type_ResidueDesignType_Catalytic;
        }
      }
      StringArrayDestroy(&strings);
    }
    else if(strcmp(keyword, "mutated_sites_begin") == 0){
      StringArray strings;
      StringArrayCreate(&strings);
      while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
        StringArraySplitString(&strings, buffer, ' ');
        if(strcmp(StringArrayGet(&strings, 0), "mutated_sites_end") == 0){
          break;
        }
        char* chnname = StringArrayGet(&strings, 0);
        int posInChain = atoi(StringArrayGet(&strings, 1));
        Chain* pChain = StructureFindChainByName(pThis, chnname);
        int resiIndex = -1;
        ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
        if(resiIndex != -1){
          ChainGetResidue(pChain, resiIndex)->designSiteType = Type_ResidueDesignType_Mutated;
        }
      }
      StringArrayDestroy(&strings);
    }
    else if(strcmp(keyword, "rotameric_sites_begin") == 0){
      StringArray strings;
      StringArrayCreate(&strings);
      while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
        StringArraySplitString(&strings, buffer, ' ');
        if(strcmp(StringArrayGet(&strings, 0), "rotameric_sites_end") == 0){
          break;
        }
        char* chnname = StringArrayGet(&strings, 0);
        int posInChain = atoi(StringArrayGet(&strings, 1));
        Chain* pChain = StructureFindChainByName(pThis, chnname);
        int resiIndex = -1;
        ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
        if(resiIndex != -1){
          ChainGetResidue(pChain, resiIndex)->designSiteType = Type_ResidueDesignType_Rotameric;
        }
      }
      StringArrayDestroy(&strings);
    }
  }
  FileReaderDestroy(&fr);

  for(int i = 0; i < StructureGetChainCount(pThis); i++){
    Chain *pChainI = StructureGetChain(pThis, i);
    if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
    for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if(pResidue->designSiteType != Type_ResidueDesignType_Fixed){
        //printf("generate rotamers for residue %d\n", j);
        ProteinSiteBuildSpecifiedRotamers(pThis, i, j, rotlib, atomParams, resiTopos, pResidue->designSiteType);
        if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamer(pThis, i, j, resiTopos);
        if(FLAG_EXPAND_HYDROXYL_ROT==TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
      }
    }
  }

  return Success;
}


int StructureGenerateBindingSiteRotamers(Structure* pThis,BBindRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,double range){
  char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  Residue* pSmallMol = NULL;
  StructureFindSmallMol(pThis, &pSmallMol);
  if(pSmallMol==NULL){
    int code = DataNotExistError;
    sprintf(errMsg, "in file %s function %s() line %d, cannot find small molecule", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, code);
    return code;
  }

  for(int i = 0; i < StructureGetChainCount(pThis); i++){
    Chain *pChainI = StructureGetChain(pThis, i);
    if(strcmp(ResidueGetChainName(pSmallMol),ChainGetName(pChainI))==0) continue;
    if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
    for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if(AtomArrayCalcMinDistance(&pSmallMol->atoms, &pResidue->atoms)>range) continue;
      if(pResidue->designSiteType != Type_ResidueDesignType_Fixed) continue;
      ResidueSetDesignSiteFlag(pResidue, Type_ResidueDesignType_Mutated);
      ProteinSiteBuildSpecifiedRotamers(pThis, i, j, rotlib, atomParams, resiTopos, pResidue->designSiteType);
      if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamer(pThis, i, j, resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT==TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
    }
  }

  return Success;
}


int StructureGenerateAllRotamers(Structure* pThis,BBindRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* DESCHNS){
  //first, create rotamers for mutated sites
  for(int i = 0; i < StructureGetChainCount(pThis); i++){
    Chain *pChainI = StructureGetChain(pThis, i);
    if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if(strstr(DESCHNS,ChainGetName(pChainI))==0) continue;
    for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if(pResidue->designSiteType != Type_ResidueDesignType_Fixed) continue;
      ResidueSetDesignSiteFlag(pResidue, Type_ResidueDesignType_Mutated);
      //ResidueSetDesignSiteFlag(pResidue, Type_ResidueDesignType_Rotameric);
      ProteinSiteBuildSpecifiedRotamers(pThis, i, j, rotlib, atomParams, resiTopos, pResidue->designSiteType);
      if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamer(pThis, i, j, resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT==TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
    }
  }
  //second, create rotamers for rotameric sites => only change conformations
  for(int i = 0; i < StructureGetChainCount(pThis); i++){
    Chain *pChainI = StructureGetChain(pThis, i);
    if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if(strstr(DESCHNS,ChainGetName(pChainI))!=0) continue;
    for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if(pResidue->designSiteType != Type_ResidueDesignType_Fixed) continue;
      Atom* pAtomCA1=ResidueGetAtomByName(pResidue,"CA");
      BOOL interResi=FALSE;
      for(int k=0; k<StructureGetChainCount(pThis); k++){
        if(k==i)continue;
        Chain* pChainK=StructureGetChain(pThis,k);
        if(ChainGetType(pChainK)!=Type_Chain_Protein)continue;
        if(strstr(DESCHNS,ChainGetName(pChainK))==0)continue;
        for(int s=0;s<ChainGetResidueCount(pChainK);s++){
          Residue* pResidueKS=ChainGetResidue(pChainK,s);
          Atom* pAtomCA2=ResidueGetAtomByName(pResidueKS,"CA");
          if(XYZDistance(&pAtomCA2->xyz,&pAtomCA1->xyz)>15.0) continue;
          if(AtomArrayCalcMinDistance(&pResidue->atoms,&pResidueKS->atoms)<PPI_DIST_CUTOFF){
            interResi=TRUE;
            break;
          }
        }
        if(interResi==TRUE) break;
      }
      if(interResi==TRUE){
        ResidueSetDesignSiteFlag(pResidue, Type_ResidueDesignType_Rotameric);
        ProteinSiteBuildSpecifiedRotamers(pThis, i, j, rotlib, atomParams, resiTopos, pResidue->designSiteType);
        if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamer(pThis, i, j, resiTopos);
        if(FLAG_EXPAND_HYDROXYL_ROT==TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
      }
    }
  }

  return Success;
}


int StructureGenerateWildtypeRotamers(Structure* pThis,BBindRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* DESCHNS){
  //first, create rotamers for mutated sites
  for(int i = 0; i < StructureGetChainCount(pThis); i++){
    Chain *pChainI = StructureGetChain(pThis, i);
    if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if(strstr(DESCHNS,ChainGetName(pChainI))==0) continue;
    for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if(pResidue->designSiteType != Type_ResidueDesignType_Fixed) continue;
      ResidueSetDesignSiteFlag(pResidue, Type_ResidueDesignType_Rotameric);
      ProteinSiteBuildSpecifiedRotamers(pThis,i,j,rotlib,atomParams,resiTopos,pResidue->designSiteType);
      if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamer(pThis,i,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT==TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
    }
  }
  //second, create rotamers for rotameric sites => only change conformations
  for(int i = 0; i < StructureGetChainCount(pThis); i++){
    Chain *pChainI = StructureGetChain(pThis, i);
    if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if(strstr(DESCHNS,ChainGetName(pChainI))!=0) continue;
    for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if(pResidue->designSiteType != Type_ResidueDesignType_Fixed) continue;
      Atom* pAtomCA1=ResidueGetAtomByName(pResidue,"CA");
      BOOL interResi=FALSE;
      for(int k=0; k<StructureGetChainCount(pThis); k++){
        if(k==i)continue;
        Chain* pChainK=StructureGetChain(pThis,k);
        if(ChainGetType(pChainK)!=Type_Chain_Protein)continue;
        if(strstr(DESCHNS,ChainGetName(pChainK))==0)continue;
        for(int s=0;s<ChainGetResidueCount(pChainK);s++){
          Residue* pResidueKS=ChainGetResidue(pChainK,s);
          Atom* pAtomCA2=ResidueGetAtomByName(pResidueKS,"CA");
          if(XYZDistance(&pAtomCA2->xyz,&pAtomCA1->xyz)>15.0) continue;
          if(AtomArrayCalcMinDistance(&pResidue->atoms,&pResidueKS->atoms)<PPI_DIST_CUTOFF){
            interResi=TRUE;
            break;
          }
        }
        if(interResi==TRUE) break;
      }
      if(interResi==TRUE){
        ResidueSetDesignSiteFlag(pResidue, Type_ResidueDesignType_Rotameric);
        ProteinSiteBuildSpecifiedRotamers(pThis,i,j,rotlib,atomParams,resiTopos,pResidue->designSiteType);
        if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamer(pThis,i,j,resiTopos);
        if(FLAG_EXPAND_HYDROXYL_ROT==TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
      }
    }
  }

  return Success;
}



//#define DEAL_WITH_LIGAND_ROTAMERS
//
//int StructureGetTruncatedBackbone(Structure* pThis, Residue* pSmallMol, double activeSiteRange, BOOL withHydrogen,AtomArray* pBackboneAtoms){
//  int chainIndex;
//  int resiIndex;
//
//  for(chainIndex=0; chainIndex<StructureGetChainCount(pThis); chainIndex++){
//    Chain* pChain = StructureGetChain(pThis,chainIndex);
//    for(resiIndex=0; resiIndex< ChainGetResidueCount(pChain); resiIndex++ ){
//      int atomIndex;
//      double minDist;
//      Residue* pResi = ChainGetResidue(pChain,resiIndex);
//
//      //Make sure it's not the small mol itself
//      if( strcmp(ResidueGetChainName(pResi),ResidueGetChainName(pSmallMol))==0 &&
//        ResidueGetPosInChain(pResi) == ResidueGetPosInChain(pSmallMol) ){
//          continue;
//      }
//
//      minDist = AtomArrayCalcMinDistance( ResidueGetAllAtoms(pResi),ResidueGetAllAtoms(pSmallMol));
//
//      if( minDist > activeSiteRange ){
//        continue;
//      }
//
//      for(atomIndex=0;atomIndex<ResidueGetAtomCount(pResi);atomIndex++){
//        Atom* pAtom = ResidueGetAtom(pResi,atomIndex);
//        if( AtomIsHydrogen(pAtom) && !withHydrogen){
//          continue;
//        }
//        if(pAtom->isBBAtom==FALSE && strcmp(AtomGetName(pAtom), "CB") != 0){
//          continue;
//        }
//        AtomArrayAppend(pBackboneAtoms,pAtom);
//      }
//    }
//  }
//  return Success;
//}
//
//int StructureDeployCataConsSitePair(Structure* pThis, CataConsSitePair* pCataConsSitePair){
//  char errMsg[MAX_LENGTH_ERR_MSG+1];
//  Residue* pSmallMol = NULL;
//
//  RotamerSet* pRotSetOfFirstAndSecondSite[2] = {NULL,NULL};
//
//  RotamerSet tempSmallMolRotamerSet;
//  RotamerSetCreate(&tempSmallMolRotamerSet);
//
//  int result = StructureFindSmallMol(pThis,&pSmallMol);
//  if(FAILED(result)){
//    TraceError("In StructureDeployCataConsSitePair(), Cannot Find Small Molecule",result);
//    return result;
//  }else{
//    Rotamer tempSmallMolRotamer;
//    RotamerCreate(&tempSmallMolRotamer);
//    RotamerSetType(&tempSmallMolRotamer,ResidueGetName(pSmallMol));
//    RotamerAddAtoms(&tempSmallMolRotamer,ResidueGetAllAtoms(pSmallMol));
//    RotamerSetAdd(&tempSmallMolRotamerSet,&tempSmallMolRotamer);
//    RotamerDestroy(&tempSmallMolRotamer);
//  }
//
//  for(int i=0;i<2;i++){
//    int posInChain;
//    char* chainName;
//    char* residueName;
//    if(i==0){
//      posInChain = pCataConsSitePair->firstSitePosInChain;
//      chainName  = pCataConsSitePair->firstSiteChainName;
//      residueName  = pCataConsSitePair->firstSiteResidueName;
//    }
//    else{
//      posInChain = pCataConsSitePair->secondSitePosInChain;
//      chainName  = pCataConsSitePair->secondSiteChainName;
//      residueName  = pCataConsSitePair->secondSiteResidueName;
//    }
//
//    pRotSetOfFirstAndSecondSite[i] = NULL;
//    if( posInChain == ResidueGetPosInChain(pSmallMol) &&
//      strcmp(chainName,ResidueGetChainName(pSmallMol))==0 &&
//      strcmp(residueName,ResidueGetName(pSmallMol))==0 ){
//        pRotSetOfFirstAndSecondSite[i] = &tempSmallMolRotamerSet;
//    }
//    else{
//      int chainIndex=-1, resiIndex = -1;
//      StructureFindChainIndex(pThis, chainName, &chainIndex);
//      ChainFindResidueByPosInChain(StructureGetChain(pThis, chainIndex), posInChain, &resiIndex);
//      DesignSite* pDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
//      pRotSetOfFirstAndSecondSite[i] = DesignSiteGetRotamers(pDesignSite);
//    }
//    if(pRotSetOfFirstAndSecondSite[i]==NULL){
//      result = DataNotExistError;
//      sprintf(errMsg,"In StructureDeployCataConsSitePair(), cannot find site %s %s %d",
//        chainName,residueName,posInChain);
//      TraceError(errMsg,result);
//      return result;
//    }
//  }
//
//  result = CataConsSitePairDeploy(pCataConsSitePair,pRotSetOfFirstAndSecondSite[0],pRotSetOfFirstAndSecondSite[1]);
//
//  if(FAILED(result)){
//    sprintf(errMsg,"In StructureDeployCataConsSitePair(), cannot deploy CataCons "
//      "between Site %s %s %d and %s %s %d",
//      pCataConsSitePair->firstSiteChainName,
//      pCataConsSitePair->firstSiteResidueName,
//      pCataConsSitePair->firstSitePosInChain,
//      pCataConsSitePair->secondSiteChainName,
//      pCataConsSitePair->secondSiteResidueName,
//      pCataConsSitePair->secondSitePosInChain);
//    TraceError(errMsg,result);
//    return result;
//  }
//
//  RotamerSetDestroy(&tempSmallMolRotamerSet);
//  return Success;
//}
//
//int StructurePlaceSmallMol(Structure* pThis, PlacingRule* pPlacingRule, CataConsSitePairArray* pCataConsCollection, int relatedProteinSiteCount,DesignSite** relatedProteinSites, RotamerSet* pSmallMolRotSet){
//  int    result;
//  char errMsg[MAX_LENGTH_ERR_MSG+1];
//  DesignSite* pPlacingStartingSite;
//
//  //Find the placing starting site in pThis->designSites
//  pPlacingStartingSite = NULL;
//  int chainIndex=-1, resiIndex = -1;
//  StructureFindChainIndex(pThis, PlacingRuleGetChainName(pPlacingRule), &chainIndex);
//  ChainFindResidueByPosInChain(StructureGetChain(pThis, chainIndex), PlacingRuleGetPosInChain(pPlacingRule), &resiIndex);
//  pPlacingStartingSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
//  if(pPlacingStartingSite == NULL){
//    result = DataNotExistError;
//    sprintf(errMsg,"In StructurePlaceSmallMol(), the placing rule requires a design site %s %d on "
//      "chain %s. But there is no such design site",
//      PlacingRuleGetResiName(pPlacingRule),
//      PlacingRuleGetPosInChain(pPlacingRule),
//      PlacingRuleGetChainName(pPlacingRule));
//    TraceError(errMsg,result);
//    return result;
//  }
//
//  //All ready, begin calculation
//  result = PlacingRulePlaceSmallMol(pPlacingRule,pCataConsCollection,pPlacingStartingSite, relatedProteinSiteCount,relatedProteinSites,pSmallMolRotSet);
//
//  if(FAILED(result)){
//    return result;
//  }
//
//  return Success;
//}
//
//// this is the main function to generate smallmol rotamers
//int StructureGenerateSmallMolRotamers(Structure* pThis, char* cataConsFileName, char* placingRuleFileName){
//  int    result;
//  Residue* pSmallMol = NULL;
//  AtomArray truncatedBackbone;
//  CataConsSitePairArray cataCons;
//  PlacingRule placingRule; 
//
//  printf("make ligand ensemble for ligand conformational sampling in computational enzyme design\n");
//
//  //1. Create and deploy the catalytic constraints
//  printf("read and deploy catalytic constraints\n");
//  result = CataConsSitePairArrayCreate(&cataCons,cataConsFileName);
//  if(FAILED(result)){
//    return result;
//  }
//  for(int i=0;i<CataConsSitePairArrayGetCount(&cataCons);i++){
//    result = StructureDeployCataConsSitePair(pThis,CataConsSitePairArrayGet(&cataCons,i));
//    if(FAILED(result)){
//      return result;
//    }
//  }
//
//  //2. Find the small molecule
//  printf("find small molecule ligand\n");
//  result = StructureFindSmallMol(pThis,&pSmallMol);
//  if(FAILED(result)){
//    return result;
//  }
//
//  //3. Create and deploy the placing rule, with creation of truncated backbone
//  printf("read and deploy ligand placing process\n");
//  result = PlacingRuleCreate(&placingRule,placingRuleFileName);
//  if(FAILED(result)){
//    return result;
//  }
//
//  printf("get truncated backbone for active site pocket where we want to make the ensemble\n");
//  AtomArrayCreate(&truncatedBackbone);
//  StructureGetTruncatedBackbone(pThis,pSmallMol,PlacingRuleGetTruncatedBackboneRange(&placingRule),TRUE,&truncatedBackbone);
//  //catalytic design sites related to generate smallmole rotamers
//  int relatedProteinSiteCount = 0;
//  DesignSite** ppRelatedDesignSites = NULL;
//  for(int i=0; i<CataConsSitePairArrayGetCount(&cataCons); i++){
//    CataConsSitePair* pSitePair = CataConsSitePairArrayGet(&cataCons, i);
//    for(int k=0; k<2; k++){
//      DesignSite* pSite = NULL;
//      if(k==0){
//        if(ResidueGetPosInChain(pSmallMol)==pSitePair->firstSitePosInChain && strcmp(ResidueGetChainName(pSmallMol), pSitePair->firstSiteChainName)==0){
//          continue;
//        }
//        int chainIndex=-1, resiIndex = -1;
//        StructureFindChainIndex(pThis, pSitePair->firstSiteChainName, &chainIndex);
//        ChainFindResidueByPosInChain(StructureGetChain(pThis, chainIndex), pSitePair->firstSitePosInChain, &resiIndex);
//        pSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
//      }
//      else{
//        if(ResidueGetPosInChain(pSmallMol)==pSitePair->secondSitePosInChain && strcmp(ResidueGetChainName(pSmallMol), pSitePair->secondSiteChainName)==0){
//          continue;
//        }
//        int chainIndex=-1, resiIndex = -1;
//        StructureFindChainIndex(pThis, pSitePair->secondSiteChainName, &chainIndex);
//        ChainFindResidueByPosInChain(StructureGetChain(pThis, chainIndex), pSitePair->secondSitePosInChain, &resiIndex);
//        pSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
//      }
//      //smallmol design site is excluded
//      BOOL designsiteExist = FALSE;
//      for(int j = 0; j < relatedProteinSiteCount; j++){
//        DesignSite* pSite1 = ppRelatedDesignSites[j];
//        //skip the identical sites
//        if(pSite1 == pSite){
//          designsiteExist = TRUE;
//          break;
//        }
//      }
//      if(designsiteExist == FALSE){
//        relatedProteinSiteCount++;
//        ppRelatedDesignSites = (DesignSite**)realloc(ppRelatedDesignSites, sizeof(DesignSite*)*relatedProteinSiteCount);
//        ppRelatedDesignSites[relatedProteinSiteCount-1] = pSite;
//      }
//    }
//  }
//
//  result = PlacingRuleDeploy(&placingRule,pSmallMol,&cataCons,relatedProteinSiteCount, ppRelatedDesignSites,&truncatedBackbone);
//
//  if(FAILED(result)){
//    return result;
//  }
//
//  //4. Prepare the small mol's rotamer set for outputting
//  RotamerSet smallMolRotamerSet;
//  RotamerSetCreate(&smallMolRotamerSet);
//
//  //5. Place small molecule
//  printf("making ligand ensembles\n");
//  result = StructurePlaceSmallMol(pThis,&placingRule,&cataCons,relatedProteinSiteCount, ppRelatedDesignSites,&smallMolRotamerSet);
//  if(FAILED(result)){
//    return result;
//  }
//
//  //6. If the small molecule has already been set as a design site, add rotamers to this design site.
//  printf("add ligand ensemble to a new rotamer set\n");
//  BOOL smallmolRotSetGenerated = FALSE;
//  for(int i=0;i<relatedProteinSiteCount;i++){
//    if(ppRelatedDesignSites[i]->pResidue == pSmallMol){
//      int j;
//      for(j=0;j<RotamerSetGetCount(&smallMolRotamerSet);j++){
//        RotamerSetAdd(DesignSiteGetRotamers(ppRelatedDesignSites[i]), RotamerSetGet(&smallMolRotamerSet,j) );
//      }
//      smallmolRotSetGenerated = TRUE;
//      break;
//    }
//  }
//  //add a new smallmol design site if it doesn't exist
//  if(smallmolRotSetGenerated == FALSE){
//    int chainIndex=-1, resiIndex = -1;
//    StructureFindChainIndex(pThis, ResidueGetChainName(pSmallMol), &chainIndex);
//    ChainFindResidueByPosInChain(StructureGetChain(pThis, chainIndex), ResidueGetPosInChain(pSmallMol), &resiIndex);
//    ProteinSiteAddDesignSite(pThis, chainIndex, resiIndex);
//    RotamerSetCopy(DesignSiteGetRotamers(StructureGetDesignSite(pThis, pThis->designSiteCount-1)), &smallMolRotamerSet);
//  }
//
//  //1.2.3.4. Destroyed    
//  RotamerSetDestroy(&smallMolRotamerSet);
//  PlacingRuleDestroy(&placingRule);
//  CataConsSitePairArrayDestroy(&cataCons);
//  AtomArrayDestroy(&truncatedBackbone);
//  free(ppRelatedDesignSites);
//  ppRelatedDesignSites = NULL;
//
//  return Success;
//}
//
//int StructureReadSmallMolRotamers(Structure* pThis,ResiTopoSet* resiTopos,char* smallMolFileName){
//  char errMsg[MAX_LENGTH_ERR_MSG+1];
//  FileReader smallMolFile;
//  Residue* pSmallMol = NULL;
//
//  int result = StructureFindSmallMol(pThis, &pSmallMol);
//  if(FAILED(result)){
//    sprintf(errMsg,"In StructureReadSmallMolRotamer, cannot file small molecule");
//    TraceError(errMsg,result);
//    return result;
//  }
//
//  int chainIndex=-1, resiIndex = -1;
//  StructureFindChainIndex(pThis, ResidueGetChainName(pSmallMol), &chainIndex);
//  ChainFindResidueByPosInChain(StructureGetChain(pThis, chainIndex), ResidueGetPosInChain(pSmallMol), &resiIndex);
//  DesignSite* pSmallMolSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
//  //the smallmol design site is not defined yet
//  if(pSmallMolSite==NULL){
//    ProteinSiteAddDesignSite(pThis,chainIndex,resiIndex);
//    //pSmallMolSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
//    pSmallMolSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
//  }
//  //set smallmol design type
//  pSmallMol->designSiteType = Type_ResidueDesignType_SmallMol;
//
//
//  result = FileReaderCreate(&smallMolFile,smallMolFileName);
//  if(FAILED(result)){
//    sprintf(errMsg,"In StructureReadSmallMolRotamer, cannot create FileReader for %s",smallMolFileName);
//    TraceError(errMsg,result);
//    return result;
//  }
//  while(!FileReaderEndOfFile(&smallMolFile)){
//    int i;
//    char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
//    char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1];
//    Residue tempResidueForSmallMol;
//
//    FileReaderGetNextLine(&smallMolFile,line);
//    ExtractFirstStringFromSourceString(keyword,line);
//    if(strcmp(keyword,"MODEL")!=0){
//      continue;
//    }
//
//    ResidueCreate(&tempResidueForSmallMol);
//    ResidueCopy(&tempResidueForSmallMol, pSmallMol);
//    for(i=0;i<ResidueGetAtomCount(&tempResidueForSmallMol);i++){
//      ResidueGetAtom(&tempResidueForSmallMol,i)->isXyzValid = FALSE;
//    }
//    if( FAILED(ResidueReadXYZFromPDB(&tempResidueForSmallMol,&smallMolFile))){
//      sprintf(errMsg,"In StructureReadSmallMolRotamer, error when reading PDB file");
//      TraceError(errMsg,result);
//      return result;
//    }
//
//    ResidueCheckAtomCoordinateValidity(&tempResidueForSmallMol);
//
//    result = ResidueCalcAllAtomXYZ(&tempResidueForSmallMol,resiTopos,NULL,NULL);
//
//    if(FAILED(result)){
//      sprintf(errMsg,"In StructureReadSmallMolRotamer, not all atoms' XYZ can be calculated");
//      TraceError(errMsg,result);
//      return result;
//    }else{
//      Rotamer tempRotamer;
//      RotamerCreate(&tempRotamer);
//
//      strcpy(tempRotamer.type, ResidueGetName(&tempResidueForSmallMol));
//      RotamerAddAtoms(&tempRotamer,ResidueGetAllAtoms(&tempResidueForSmallMol));
//      BondSetCopy(RotamerGetBonds(&tempRotamer),ResidueGetBonds(&tempResidueForSmallMol));
//      RotamerSetChainName(&tempRotamer,ResidueGetChainName(&tempResidueForSmallMol));
//      //read the energy from smallmol rotamer file
//      tempRotamer.vdwInternal = tempResidueForSmallMol.internalEnergy;
//      tempRotamer.vdwBackbone = tempResidueForSmallMol.backboneEnergy;
//      //set smallmol internal energy as self energy and recalculate backbone energy later
//      //tempRotamer.selfenergy = tempResidueForSmallMol.internalEnergy;
//      RotamerSetPosInChain(&tempRotamer,ResidueGetPosInChain(&tempResidueForSmallMol));
//      RotamerSetAdd(DesignSiteGetRotamers(pSmallMolSite),&tempRotamer);
//      RotamerDestroy(&tempRotamer);
//    }
//
//    ResidueDestroy(&tempResidueForSmallMol);
//  }
//
//  FileReaderDestroy(&smallMolFile);
//  return Success;
//}
//
//int StructureWriteSmallMolRotamers(Structure* pThis,char* smallMolFile){
//  int i;
//  Residue* pSmallMol = NULL;
//  RotamerSet* pSmallMolRotamers = NULL;
//  FILE* outputFile = NULL;
//  char errMsg[MAX_LENGTH_ERR_MSG+1];
//  int result;
//
//
//  StructureFindSmallMol(pThis,&pSmallMol);
//  int chainIndex=-1, resiIndex = -1;
//  StructureFindChainIndex(pThis, ResidueGetChainName(pSmallMol), &chainIndex);
//  ChainFindResidueByPosInChain(StructureGetChain(pThis, chainIndex), ResidueGetPosInChain(pSmallMol), &resiIndex);
//  DesignSite* pSmallMolSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
//  pSmallMolRotamers = DesignSiteGetRotamers(pSmallMolSite);
//
//  if(pSmallMolRotamers==NULL){
//    result = DataNotExistError;
//    sprintf(errMsg,"In StructureWriteSmallMolRotamers(), cannot find the design site of the small molecule");
//    TraceError(errMsg,result);
//    return result;
//  }
//
//  if(smallMolFile==NULL || strcmp(smallMolFile,"")==0 ){
//    outputFile = NULL;
//  }
//  else{
//    outputFile = fopen(smallMolFile,"w");
//    if(outputFile==NULL){
//      sprintf(errMsg,"In StructureWriteSmallMolRotamers(), cannot create file %s to write",smallMolFile);
//      result = IOError;
//      TraceError(errMsg,result);
//      return result;
//    }
//  }
//
//  for(i=0;i<RotamerSetGetCount(pSmallMolRotamers);i++){
//    Rotamer tempRotamer;
//    RotamerCreate(&tempRotamer);
//    Model(i,outputFile);
//    RotamerCopy(&tempRotamer,RotamerSetGet(pSmallMolRotamers,i));
//    RotamerRestore(&tempRotamer,pSmallMolRotamers);
//    RotamerShowInPDBFormat(&tempRotamer,"ATOM",ResidueGetChainName(pSmallMol),0,i%10000,FALSE,outputFile);
//    RotamerExtract(&tempRotamer);
//    EndModel(outputFile);
//    fprintf(outputFile,"ENERGY INTERNAL: %f BACKBONE: %f\n",tempRotamer.vdwInternal,tempRotamer.vdwBackbone);
//    RotamerDestroy(&tempRotamer);
//  }
//
//  return Success;
//}
//
//
//// this function is used to delete or screen small molecules with correct direct orientation;
//// oriFileName is small molecule library file by combining the tmpPDB files together;
//// newFileName is the reserved library file after screening;
//// screenFileName is the screening rule file;
//// the following information is needed in the screening rule file:
//// two atoms on small molecule: a catalytic atom and a binding atom;
//// a set of residues in protein scaffold;
//// the information is organized in the following format:
//// CATA_ATOM  C15
//// BIND_ATOM  C11
//// BIND_SITE_GROUP
//// BIND_SITE  CHAINB VAL 55
//// BIND_SITE_GROUP
//// BIND_SITE  CHAINB PHE 56
//// BIND_SITE_GROUP
//// BIND_SITE  CHAINB SER 66
//// BIND_SITE_GROUP
//// BIND_SITE  CHAINB TRP 153
//// in each BIND_SITE_GROUP, all the distances between the binding atom and the CA atoms of the residues must be are 
//// larger than those between the catalytic atom and the CA atoms of the residues;
//// and at least one binding site group must satisfy the above relationship.
//int StructureSmallmolOrientationScreen(Structure* pStructure, ResiTopoSet* pResiTopo, char* oriFileName, char* newFileName, char* screenRuleFileName)
//{
//  typedef struct _BindSite{
//    char chainName[MAX_LENGTH_CHAIN_NAME+1];
//    int  posInChain;
//    char resiName[MAX_LENGTH_RESIDUE_NAME+1];
//  } BindSite;
//
//  typedef struct _BindSiteGroup{
//    int       siteNum;
//    BindSite  *pSites;
//  } BindSiteGroup;
//
//  FileReader     file;
//  int            i, j;
//  int            totalRotamerCounter, acceptedRotamerCounter;
//  char           cataAtomName[MAX_LENGTH_ATOM_NAME+1];
//  char           bindAtomName[MAX_LENGTH_ATOM_NAME+1];
//  int            groupCounter;
//  char           line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
//  BindSiteGroup* pGroups = NULL;
//  FILE*          pOutputFile = NULL;
//  FILE*          pInputFile  = NULL;
//  StringArray    buffer;
//  XYZ            xyzCataAtom;
//  XYZ            xyzBindAtom;
//  BOOL           cataAtomExist = FALSE;
//  BOOL           bindAtomExist = FALSE;
//  BOOL           curRotamerAccepted = FALSE;
//
//  FileReaderCreate(&file, screenRuleFileName);
//  groupCounter = 0;
//  while(!FileReaderEndOfFile(&file)){
//    BOOL doneInThisGroup = FALSE;
//    int siteNum;
//    while(!FAILED(FileReaderGetNextLine(&file, line))){
//      StringArray wordsInLine;
//      StringArrayCreate(&wordsInLine);
//      StringArraySplitString(&wordsInLine, line, ' ');
//
//      if(strcmp(StringArrayGet(&wordsInLine, 0), "CATA_ATOM") == 0){
//        strcpy(cataAtomName, StringArrayGet(&wordsInLine, 1));
//      }
//      else if(strcmp(StringArrayGet(&wordsInLine, 0), "BIND_ATOM") == 0){
//        strcpy(bindAtomName, StringArrayGet(&wordsInLine, 1));
//      }
//      else if(strcmp(StringArrayGet(&wordsInLine, 0), "BIND_SITE_GROUP") == 0 && doneInThisGroup == TRUE){
//        FileReaderSetCurrentPos(&file, FileReaderGetCurrentPos(&file)-1);
//        break;
//      }
//      else if(strcmp(StringArrayGet(&wordsInLine, 0), "BIND_SITE_GROUP") == 0 && doneInThisGroup == FALSE){
//        doneInThisGroup = TRUE;
//        groupCounter++;
//        pGroups = (BindSiteGroup*)realloc(pGroups, sizeof(BindSiteGroup)*groupCounter);
//        siteNum=0;
//        pGroups[groupCounter-1].siteNum = 0;
//        pGroups[groupCounter-1].pSites = NULL;
//        continue;
//      }
//      else if(strcmp(StringArrayGet(&wordsInLine, 0), "BIND_SITE") == 0){
//        siteNum++;
//        pGroups[groupCounter-1].pSites = (BindSite*)realloc(pGroups[groupCounter-1].pSites, sizeof(BindSite)*siteNum);
//        strcpy(pGroups[groupCounter-1].pSites[siteNum-1].chainName, StringArrayGet(&wordsInLine, 1));
//        strcpy(pGroups[groupCounter-1].pSites[siteNum-1].resiName, StringArrayGet(&wordsInLine, 2));
//        pGroups[groupCounter-1].pSites[siteNum-1].posInChain = atoi(StringArrayGet(&wordsInLine, 3));
//        pGroups[groupCounter-1].siteNum = siteNum;
//      }
//    }
//  }
//
//  char errMsg[MAX_LENGTH_ERR_MSG+1];
//  printf("reading small-molecule rotamers from file %s\n", oriFileName);
//  pInputFile = fopen(oriFileName, "r");
//  if(pInputFile == NULL){
//    sprintf(errMsg,"in file %s function %s line %d, can not open file %s for reading",oriFileName);
//    TraceError(errMsg,IOError);
//    return IOError;
//  }
//  pOutputFile = fopen(newFileName, "w");
//  if(pOutputFile == NULL){
//    sprintf(errMsg,"in file %s function %s line %d, can not open file %s for reading",newFileName);
//    TraceError(errMsg,IOError);
//    return IOError;
//  }
//
//  totalRotamerCounter = 0;
//  acceptedRotamerCounter = 0;
//  StringArrayCreate(&buffer);
//  while(fgets(line, MAX_LENGTH_ONE_LINE_IN_FILE, pInputFile)){
//    char keyword[10];
//
//    strcpy(keyword,"");
//    ExtractTargetStringFromSourceString(keyword,line,0,4);
//
//    if(strcmp(keyword,"ENDM")==0){
//      if(cataAtomExist == FALSE || bindAtomExist == FALSE){
//        sprintf(errMsg,"in file %s function %s line %d, CataAtom or BindAtom does not exist in file %s",oriFileName);
//        TraceError(errMsg,FormatError);
//        return FormatError;
//      }
//
//      for(i = 0; i < groupCounter; i++){
//        curRotamerAccepted = TRUE;
//        for(j = 0; j < pGroups[i].siteNum; j++){
//          Residue* pResidue = ChainGetResidue(StructureFindChainByName(pStructure, pGroups[i].pSites[j].chainName), pGroups[i].pSites[j].posInChain);
//          XYZ* pXyzCA = &ResidueGetAtomByName(pResidue, "CA")->xyz;
//          if(XYZDistance(&xyzBindAtom, pXyzCA) > XYZDistance(&xyzCataAtom, pXyzCA)){
//            curRotamerAccepted = FALSE;
//            break;
//          }
//        }
//        if(curRotamerAccepted == TRUE){
//          fprintf(pOutputFile,"MODEL     %d\n", acceptedRotamerCounter);
//          for(j=0;j<StringArrayGetCount(&buffer);j++){
//            fprintf(pOutputFile,"%s",StringArrayGet(&buffer,j));
//          }
//          fprintf(pOutputFile,"ENDMDL\n");
//          acceptedRotamerCounter++;
//          break;
//        }
//      }
//
//      printf("%d / %d rotamers have been processed      \r", acceptedRotamerCounter, totalRotamerCounter);
//      //fflush(stdout);
//    }
//    else if(strcmp(keyword,"MODE")==0){
//      StringArrayDestroy(&buffer);
//      StringArrayCreate(&buffer);
//      totalRotamerCounter++;
//    }
//    else if(strcmp(keyword,"ATOM")==0){
//      char atomName[MAX_LENGTH_ATOM_NAME+1];
//      ExtractTargetStringFromSourceString(atomName,line,12,4);
//      if(strcmp(atomName, cataAtomName) == 0){
//        cataAtomExist = TRUE;
//        ExtractTargetStringFromSourceString(keyword,line,31,7);
//        xyzCataAtom.X = atof(keyword);
//        ExtractTargetStringFromSourceString(keyword,line,39,7);
//        xyzCataAtom.Y = atof(keyword);
//        ExtractTargetStringFromSourceString(keyword,line,47,7);
//        xyzCataAtom.Z = atof(keyword);
//      }
//
//      if(strcmp(atomName, bindAtomName) == 0){
//        bindAtomExist = TRUE;
//        ExtractTargetStringFromSourceString(keyword,line,31,7);
//        xyzBindAtom.X = atof(keyword);
//        ExtractTargetStringFromSourceString(keyword,line,39,7);
//        xyzBindAtom.Y = atof(keyword);
//        ExtractTargetStringFromSourceString(keyword,line,47,7);
//        xyzBindAtom.Z = atof(keyword);
//      }
//      StringArrayAppend(&buffer,line);
//    }
//    else if(strcmp(keyword,"ENER")==0 && curRotamerAccepted==TRUE){
//      fprintf(pOutputFile,"%s",line);
//    }
//  }
//  printf("\n");
//  fclose(pInputFile);
//  fclose(pOutputFile);
//
//  FileReaderDestroy(&file);
//  StringArrayDestroy(&buffer);
//  for(j = 0; j < groupCounter; j++){
//    free(pGroups[j].pSites);
//    pGroups[j].siteNum = 0;
//  }
//  free(pGroups);
//
//  return Success;
//}


#define DEAL_WITH_PROTEIN_ROTAMERS_BBDEP
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//build backbone dependent rotamers
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ProteinSiteBuildAllRotamersByBBdepRotLib(Structure* pThis,int chainIndex,int resiIndex,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos){
  int result = Success;
  DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if(pCurrentDesignSite != NULL){
    DesignSiteRemoveRotamers(pCurrentDesignSite);
  }
  else{
    (pThis->designSiteCount)++;
    pThis->designSites=(DesignSite*)realloc(pThis->designSites, sizeof(DesignSite)*pThis->designSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->designSiteCount-1]);
    pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
    Chain* pDestChain = StructureGetChain(pThis, chainIndex);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
    pCurrentDesignSite->pResidue = pDestResidue;
    pCurrentDesignSite->chainIndex = chainIndex;
    pCurrentDesignSite->resiIndex = resiIndex;
  }

  // set design types - 20 AA types;
  StringArray designTypes;
  StringArray patchTypes;
  StringArrayCreate(&designTypes);
  StringArrayCreate(&patchTypes);
  StringArrayAppend(&designTypes, "ALA"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ARG"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ASN"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ASP"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "CYS"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "GLN"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "GLU"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "GLY"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "HSD"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "HSE"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ILE"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "LEU"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "LYS"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "MET"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "PHE"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "PRO"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "SER"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "THR"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "TRP"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "TYR"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "VAL"); StringArrayAppend(&patchTypes, "");
  result = RotamerSetOfProteinGenerateByBBdepRotLib(DesignSiteGetRotamers(pCurrentDesignSite),pCurrentDesignSite->pResidue,&designTypes,&patchTypes,rotlib,atomParams,resiTopos);
  //ResidueSetDesignSiteFlag(pCurrentDesignSite->pResidue, Type_ResidueDesignType_Mutated);
  StringArrayDestroy(&patchTypes);
  StringArrayDestroy(&designTypes);

  return result;
}


int ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(Structure* pThis, int chainIndex, int resiIndex,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos, Type_ResidueDesignType type){
  int result = Success;
  DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if(pCurrentDesignSite != NULL){
    DesignSiteRemoveRotamers(pCurrentDesignSite);
  }
  else{
    (pThis->designSiteCount)++;
    pThis->designSites=(DesignSite*)realloc(pThis->designSites, sizeof(DesignSite)*pThis->designSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->designSiteCount-1]);
    pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
    Chain* pDestChain = StructureGetChain(pThis, chainIndex);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
    pCurrentDesignSite->pResidue = pDestResidue;
    pCurrentDesignSite->chainIndex = chainIndex;
    pCurrentDesignSite->resiIndex = resiIndex;
  }
  // set design types - 20 AA types;
  StringArray designTypes;
  StringArray patchTypes;
  StringArrayCreate(&designTypes);
  StringArrayCreate(&patchTypes);
  if(type == Type_ResidueDesignType_Mutated){
    StringArrayAppend(&designTypes, "ALA"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "ARG"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "ASN"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "ASP"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "CYS"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "GLN"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "GLU"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "GLY"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "HSD"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "HSE"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "ILE"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "LEU"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "LYS"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "MET"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "PHE"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "PRO"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "SER"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "THR"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "TRP"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "TYR"); StringArrayAppend(&patchTypes, "");
    StringArrayAppend(&designTypes, "VAL"); StringArrayAppend(&patchTypes, "");
  }
  else if(type == Type_ResidueDesignType_Rotameric || type == Type_ResidueDesignType_Catalytic){
    StringArrayAppend(&designTypes, ResidueGetName(pCurrentDesignSite->pResidue)); StringArrayAppend(&patchTypes, "");
  }
  result = RotamerSetOfProteinGenerateByBBdepRotLib(DesignSiteGetRotamers(pCurrentDesignSite),pCurrentDesignSite->pResidue,&designTypes,&patchTypes,rotlib,atomParams,resiTopos);
  StringArrayDestroy(&patchTypes);
  StringArrayDestroy(&designTypes);

  return result;
}




int ProteinSiteBuildMutatedRotamersByBBdepRotLib(Structure* pThis, int chainIndex, int resiIndex, BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos, StringArray *pDesignTypes, StringArray *pPatchTypes){
  int result = Success;
  DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if(pCurrentDesignSite != NULL){
    DesignSiteRemoveRotamers(pCurrentDesignSite);
  }
  else{
    (pThis->designSiteCount)++;
    pThis->designSites=(DesignSite*)realloc(pThis->designSites, sizeof(DesignSite)*pThis->designSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->designSiteCount-1]);
    pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
    Chain* pDestChain = StructureGetChain(pThis, chainIndex);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
    pCurrentDesignSite->pResidue = pDestResidue;
    pCurrentDesignSite->chainIndex = chainIndex;
    pCurrentDesignSite->resiIndex = resiIndex;
  }

  result = RotamerSetOfProteinGenerateByBBdepRotLib(DesignSiteGetRotamers(pCurrentDesignSite),pCurrentDesignSite->pResidue,pDesignTypes,pPatchTypes,rotlib,atomParams,resiTopos);
  ResidueSetDesignSiteFlag(pCurrentDesignSite->pResidue,Type_ResidueDesignType_Mutated);
  return result;
}

int ProteinSiteBuildWildtypeRotamersByBBdepRotLib(Structure* pThis,int chainIndex,int resiIndex,BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos){
  int result = Success;
  DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if(pCurrentDesignSite != NULL){
    DesignSiteRemoveRotamers(pCurrentDesignSite);
  }
  else{
    (pThis->designSiteCount)++;
    pThis->designSites=(DesignSite*)realloc(pThis->designSites, sizeof(DesignSite)*pThis->designSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->designSiteCount-1]);
    pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
    Chain* pDestChain = StructureGetChain(pThis, chainIndex);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
    pCurrentDesignSite->pResidue = pDestResidue;
    pCurrentDesignSite->chainIndex = chainIndex;
    pCurrentDesignSite->resiIndex = resiIndex;
  }

  // set native design types;
  StringArray designTypes;
  StringArray patchTypes;
  StringArrayCreate(&designTypes);
  StringArrayCreate(&patchTypes);
  StringArrayAppend(&designTypes, ResidueGetName(pCurrentDesignSite->pResidue));
  StringArrayAppend(&patchTypes, "");
  if(strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "HSD") == 0){
    StringArrayAppend(&designTypes, "HSE");
    StringArrayAppend(&patchTypes, "");
  }
  else if(strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "HSE") == 0){
    StringArrayAppend(&designTypes, "HSD");
    StringArrayAppend(&patchTypes, "");
  }

  result = RotamerSetOfProteinGenerateByBBdepRotLib(DesignSiteGetRotamers(pCurrentDesignSite),pCurrentDesignSite->pResidue,&designTypes,&patchTypes,rotlib,atomParams,resiTopos);
  ResidueSetDesignSiteFlag(pCurrentDesignSite->pResidue, Type_ResidueDesignType_Rotameric);

  StringArrayDestroy(&patchTypes);
  StringArrayDestroy(&designTypes);

  return result;
}


// this function can be used to build crystal rotamers for every amino acid type
int ProteinSiteAddCrystalRotamerByBBdepRotLib(Structure* pThis,int chainIndex,int resiIndex,ResiTopoSet *pResiTopos,BBdepRotamerLib* pBBdepRotLib){
  Chain* pChain = StructureGetChain(pThis, chainIndex);
  Residue* pResidue = ChainGetResidue(pChain, resiIndex);
  if(pChain->type == Type_Chain_Protein){
    DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
    if(pCurrentDesignSite == NULL){
      (pThis->designSiteCount)++;
      pThis->designSites=(DesignSite*)realloc(pThis->designSites, sizeof(DesignSite)*pThis->designSiteCount);
      DesignSiteCreate(&pThis->designSites[pThis->designSiteCount-1]);
      pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
      pCurrentDesignSite->pResidue = pResidue;
      pCurrentDesignSite->chainIndex = chainIndex;
      pCurrentDesignSite->resiIndex = resiIndex;
    }

    //do not add rotamer for residue ala and gly
    if(strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "ALA") == 0 || strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "GLY") == 0){
      return Success;
    }

    //step2: get the phi&psi angles for the current residue
    double DUNBRACK_ENERGY=0.0;
    int binindex=((int)(pResidue->phipsi[0]+180)/10)*36+(int)(pResidue->phipsi[1]+180)/10;
    RotLibPhiPsi* pRotLibPhiPsi=&pBBdepRotLib->rotlibphipsis[binindex];
    int rotTypeIndex=-1;
    StringArrayFind(&pRotLibPhiPsi->rotTypes,ResidueGetName(pResidue),&rotTypeIndex);
    DoubleArray* pTorsionsArrayForTypeI = pRotLibPhiPsi->torsions[rotTypeIndex];
    DoubleArray* pDeviationsArrayForTypeI = pRotLibPhiPsi->deviations[rotTypeIndex];
    int matchIndex=-1;
    for(int i=0;i<IntArrayGet(&pRotLibPhiPsi->rotamerCounts,rotTypeIndex);i++){
      DoubleArray* pTorsions=&pTorsionsArrayForTypeI[i];
      DoubleArray* pDeviations=&pDeviationsArrayForTypeI[i];
      BOOL match=TRUE;
      for(int j=0;j<DoubleArrayGetLength(&pResidue->xtorsions);j++){
        double min=DoubleArrayGet(pTorsions,j)-DegToRad(30);
        double max=DoubleArrayGet(pTorsions,j)+DegToRad(30);
        double torsion=DoubleArrayGet(&pResidue->xtorsions,j);
        double torsionm2pi=torsion-2*PI;
        double torsionp2pi=torsion+2*PI;
        double torsion2=torsion;
        if((strcmp(ResidueGetName(pResidue),"PHE")==0 && j==1)||
          (strcmp(ResidueGetName(pResidue),"TYR")==0 && j==1)||
          (strcmp(ResidueGetName(pResidue),"ASP")==0 && j==1)||
          strcmp(ResidueGetName(pResidue),"GLU")==0 && j==2){
            torsion2=torsion+PI;
            torsion2=torsion>0?torsion-PI:torsion2;
        }
        double torsion2m2pi=torsion2-2*PI;
        double torsion2p2pi=torsion2+2*PI;
        if(!(
          (torsion    <=max && torsion>=min) ||
          (torsionm2pi<=max && torsionm2pi>=min) ||
          (torsionp2pi<=max && torsionp2pi>=min) ||
          (torsion2    <=max && torsion2>=min) ||
          (torsion2m2pi<=max && torsion2m2pi>=min) ||
          (torsion2p2pi<=max && torsion2p2pi>=min)
          )){
            match=FALSE;
            break;
        }
      }
      if(match==TRUE){
        matchIndex=i;
        break;
      }
    }
    double delta_prob=1e-7;
    if(matchIndex!=-1){
      DUNBRACK_ENERGY=-1.0*log(DoubleArrayGet(&pRotLibPhiPsi->probability[rotTypeIndex],matchIndex)+delta_prob);
    }
    else{
      //printf("cannot find similar rotamer for residue %s %d %s in bbdep rotlib\n",ResidueGetChainName(pDestResidue),ResidueGetPosInChain(pDestResidue),ResidueGetName(pDestResidue));
      DUNBRACK_ENERGY=-1.0*log(DoubleArrayGet(&pRotLibPhiPsi->probability[rotTypeIndex],DoubleArrayGetLength(&pRotLibPhiPsi->probability[rotTypeIndex])-1)+delta_prob);
    }

    //else add a crystal rotamer for the residue
    RotamerSet* pSetI = DesignSiteGetRotamers(pCurrentDesignSite);
    Rotamer tempRotamer;
    Rotamer* pRotamerRepresentative;
    RotamerCreate(&tempRotamer);
    RotamerSetType(&tempRotamer,ResidueGetName(pCurrentDesignSite->pResidue));
    RotamerSetChainName(&tempRotamer,ResidueGetChainName(pCurrentDesignSite->pResidue));
    RotamerSetPosInChain(&tempRotamer,ResidueGetPosInChain(pCurrentDesignSite->pResidue));
    pRotamerRepresentative = RotamerSetGetRepresentative(pSetI, RotamerGetType(&tempRotamer));
    if(pRotamerRepresentative != NULL){
      AtomArrayCopy(&tempRotamer.atoms, &pRotamerRepresentative->atoms);
      for(int i = 0; i < RotamerGetAtomCount(&tempRotamer); i++){
        Atom* pAtom = RotamerGetAtom(&tempRotamer, i);
        pAtom->xyz = ResidueGetAtomByName(pCurrentDesignSite->pResidue, AtomGetName(pAtom))->xyz;
      }
      BondSetCopy(&tempRotamer.bonds,&pRotamerRepresentative->bonds);
      XYZArrayResize(&tempRotamer.xyzs, RotamerGetAtomCount(pRotamerRepresentative));
      for(int i = 0; i < XYZArrayGetLength(&tempRotamer.xyzs); i++){
        XYZArraySet(&tempRotamer.xyzs, i, &AtomArrayGet(&tempRotamer.atoms,i)->xyz);
      }
      tempRotamer.dunbrack=DUNBRACK_ENERGY;
      DoubleArrayCopy(&tempRotamer.xtorsions,&pResidue->xtorsions);
      RotamerSetAdd(&pCurrentDesignSite->rotamers, &tempRotamer);
    }
    else{
      AtomArrayCopy(&tempRotamer.atoms, &pResidue->atoms);
      BondSetCopy(&tempRotamer.bonds,&pResidue->bonds);
      XYZArrayResize(&tempRotamer.xyzs, AtomArrayGetCount(&tempRotamer.atoms));
      for(int i = 0; i < XYZArrayGetLength(&tempRotamer.xyzs); i++){
        XYZArraySet(&tempRotamer.xyzs, i, &AtomArrayGet(&tempRotamer.atoms,i)->xyz);
      }
      tempRotamer.dunbrack=DUNBRACK_ENERGY;
      DoubleArrayCopy(&tempRotamer.xtorsions,&pResidue->xtorsions);
      RotamerSetAdd(&pCurrentDesignSite->rotamers, &tempRotamer);
    }
    RotamerDestroy(&tempRotamer);
    //if residue is histidine, add a flipped rotamer
    if(strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "HSD") == 0){
      if((pRotamerRepresentative = RotamerSetGetRepresentative(pSetI, "HSE")) != NULL){
        Residue newResi;
        ResidueCreate(&newResi);
        ResidueSetName(&newResi, "HSE");
        AtomArrayCopy(&newResi.atoms, &pRotamerRepresentative->atoms);
        for(int i = 0; i < ResidueGetAtomCount(&newResi); i++){
          Atom* pAtom = ResidueGetAtom(&newResi, i);
          if(pAtom->isBBAtom == FALSE && AtomIsHydrogen(pAtom) == TRUE){
            pAtom->isXyzValid = FALSE;
            continue;
          }
          pAtom->xyz = ResidueGetAtomByName(pCurrentDesignSite->pResidue, AtomGetName(pAtom))->xyz;
        }
        ResidueCalcAllAtomXYZ(&newResi, pResiTopos, NULL, NULL);

        RotamerCreate(&tempRotamer);
        RotamerSetType(&tempRotamer,"HSE");
        RotamerSetChainName(&tempRotamer,ResidueGetChainName(pCurrentDesignSite->pResidue));
        RotamerSetPosInChain(&tempRotamer,ResidueGetPosInChain(pCurrentDesignSite->pResidue));

        AtomArrayCopy(&tempRotamer.atoms, &pRotamerRepresentative->atoms);
        for(int i = 0; i < RotamerGetAtomCount(&tempRotamer); i++){
          Atom* pAtom = RotamerGetAtom(&tempRotamer, i);
          pAtom->xyz = ResidueGetAtomByName(&newResi, AtomGetName(pAtom))->xyz;
        }
        BondSetCopy(&tempRotamer.bonds,&pRotamerRepresentative->bonds);

        XYZArrayResize(&tempRotamer.xyzs, RotamerGetAtomCount(pRotamerRepresentative));
        for(int i = 0; i < XYZArrayGetLength(&tempRotamer.xyzs); i++){
          XYZArraySet(&tempRotamer.xyzs, i, &AtomArrayGet(&tempRotamer.atoms,i)->xyz);
        }
        tempRotamer.dunbrack=DUNBRACK_ENERGY;
        DoubleArrayCopy(&tempRotamer.xtorsions,&pResidue->xtorsions);
        RotamerSetAdd(&pCurrentDesignSite->rotamers, &tempRotamer);

        RotamerDestroy(&tempRotamer);
        ResidueDestroy(&newResi);
      }
    } // HSD
    else if(strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "HSE") == 0){
      if((pRotamerRepresentative = RotamerSetGetRepresentative(pSetI, "HSD")) != NULL){
        Residue newResi;
        ResidueCreate(&newResi);
        ResidueSetName(&newResi, "HSD");
        AtomArrayCopy(&newResi.atoms, &pRotamerRepresentative->atoms);
        for(int i = 0; i < ResidueGetAtomCount(&newResi); i++){
          Atom* pAtom = ResidueGetAtom(&newResi, i);
          if(pAtom->isBBAtom == FALSE && AtomIsHydrogen(pAtom) == TRUE){
            pAtom->isXyzValid = FALSE;
            continue;
          }
          pAtom->xyz = ResidueGetAtomByName(pCurrentDesignSite->pResidue, AtomGetName(pAtom))->xyz;
        }
        ResidueCalcAllAtomXYZ(&newResi, pResiTopos, NULL, NULL);

        RotamerCreate(&tempRotamer);
        RotamerSetType(&tempRotamer,"HSD");
        RotamerSetChainName(&tempRotamer,ResidueGetChainName(pCurrentDesignSite->pResidue));
        RotamerSetPosInChain(&tempRotamer,ResidueGetPosInChain(pCurrentDesignSite->pResidue));

        AtomArrayCopy(&tempRotamer.atoms, &pRotamerRepresentative->atoms);
        for(int i = 0; i < RotamerGetAtomCount(&tempRotamer); i++){
          Atom* pAtom = RotamerGetAtom(&tempRotamer, i);
          pAtom->xyz = ResidueGetAtomByName(&newResi, AtomGetName(pAtom))->xyz;
        }
        BondSetCopy(&tempRotamer.bonds,&pRotamerRepresentative->bonds);

        XYZArrayResize(&tempRotamer.xyzs, RotamerGetAtomCount(pRotamerRepresentative));
        for(int i = 0; i < XYZArrayGetLength(&tempRotamer.xyzs); i++){
          XYZArraySet(&tempRotamer.xyzs, i, &AtomArrayGet(&tempRotamer.atoms,i)->xyz);
        }
        tempRotamer.dunbrack=DUNBRACK_ENERGY;
        DoubleArrayCopy(&tempRotamer.xtorsions,&pResidue->xtorsions);
        RotamerSetAdd(&pCurrentDesignSite->rotamers, &tempRotamer);
        RotamerDestroy(&tempRotamer);
        ResidueDestroy(&newResi);
      }
    } //HSE
  }

  return Success;
}




int StructureGenerateSpecifiedProteinRotamersByBBdepRotLib(Structure* pThis,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* designsitefile){
  FileReader fr;
  FileReaderCreate(&fr, designsitefile);
  char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
    char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    ExtractFirstStringFromSourceString(keyword, buffer);
    if(strcmp(keyword, "catalytic_sites_begin") == 0){
      StringArray strings;
      StringArrayCreate(&strings);
      while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
        StringArraySplitString(&strings, buffer, ' ');
        if(strcmp(StringArrayGet(&strings, 0), "catalytic_sites_end") == 0){
          break;
        }
        char* chnname = StringArrayGet(&strings, 0);
        int posInChain = atoi(StringArrayGet(&strings, 1));
        Chain* pChain = StructureFindChainByName(pThis, chnname);
        int resiIndex = -1;
        ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
        if(resiIndex != -1){
          ChainGetResidue(pChain, resiIndex)->designSiteType = Type_ResidueDesignType_Catalytic;
        }
      }
      StringArrayDestroy(&strings);
    }
    else if(strcmp(keyword, "mutated_sites_begin") == 0){
      StringArray strings;
      StringArrayCreate(&strings);
      while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
        StringArraySplitString(&strings, buffer, ' ');
        if(strcmp(StringArrayGet(&strings, 0), "mutated_sites_end") == 0){
          break;
        }
        char* chnname = StringArrayGet(&strings, 0);
        int posInChain = atoi(StringArrayGet(&strings, 1));
        Chain* pChain = StructureFindChainByName(pThis, chnname);
        int resiIndex = -1;
        ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
        if(resiIndex != -1){
          ChainGetResidue(pChain, resiIndex)->designSiteType = Type_ResidueDesignType_Mutated;
        }
      }
      StringArrayDestroy(&strings);
    }
    else if(strcmp(keyword, "rotameric_sites_begin") == 0){
      StringArray strings;
      StringArrayCreate(&strings);
      while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
        StringArraySplitString(&strings, buffer, ' ');
        if(strcmp(StringArrayGet(&strings, 0), "rotameric_sites_end") == 0){
          break;
        }
        char* chnname = StringArrayGet(&strings, 0);
        int posInChain = atoi(StringArrayGet(&strings, 1));
        Chain* pChain = StructureFindChainByName(pThis, chnname);
        int resiIndex = -1;
        ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
        if(resiIndex != -1){
          ChainGetResidue(pChain, resiIndex)->designSiteType = Type_ResidueDesignType_Rotameric;
        }
      }
      StringArrayDestroy(&strings);
    }
  }
  FileReaderDestroy(&fr);

  for(int i = 0; i < StructureGetChainCount(pThis); i++){
    Chain *pChainI = StructureGetChain(pThis, i);
    if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
    for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if(pResidue->designSiteType != Type_ResidueDesignType_Fixed){
        //printf("generate rotamers for residue %d\n", j);
        ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pThis,i,j,rotlib,atomParams,resiTopos,pResidue->designSiteType);
        if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamer(pThis,i,j,resiTopos);
        if(FLAG_EXPAND_HYDROXYL_ROT=TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
      }
    }
  }

  return Success;
}


int StructureGenerateCatalyticSiteProteinRotamersByBBdepRotLib(Structure* pThis,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* catasitefile){
  FileReader fr;
  if(!FAILED(FileReaderCreate(&fr, catasitefile))){
    BOOL flagCataSite=FALSE;
    char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
      char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1];
      ExtractFirstStringFromSourceString(keyword, buffer);
      if(strcmp(keyword, "catalytic_sites_begin") == 0){
        flagCataSite=TRUE;
        StringArray strings;
        StringArrayCreate(&strings);
        while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
          StringArraySplitString(&strings, buffer, ' ');
          if(strcmp(StringArrayGet(&strings, 0), "catalytic_sites_end") == 0){
            break;
          }
          char* chnname = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pThis, chnname);
          int resiIndex = -1;
          ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
          if(resiIndex != -1){
            ChainGetResidue(pChain, resiIndex)->designSiteType = Type_ResidueDesignType_Catalytic;
          }
        }
        StringArrayDestroy(&strings);
      }
    }
    FileReaderDestroy(&fr);
    if(flagCataSite==FALSE){
      char errMsg[MAX_LENGTH_ERR_MSG+1];
      sprintf(errMsg,"in file %s function %s line %d, cannot find design sites. User should specify catalytic sites using keyword 'catalytic_sites_begein' "
        "and 'catalytic_sites_end' in the catalytic site file",__FILE__,__FUNCTION__,__LINE__);
      TraceError(errMsg,FormatError);
      return FormatError;
    }

    for(int i = 0; i < StructureGetChainCount(pThis); i++){
      Chain *pChainI = StructureGetChain(pThis, i);
      if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
      for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
        Residue* pResidue = ChainGetResidue(pChainI, j);
        if(pResidue->designSiteType != Type_ResidueDesignType_Fixed){
          //printf("generate rotamers for residue %d\n", j);
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pThis,i,j,rotlib,atomParams,resiTopos,pResidue->designSiteType);
          if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamer(pThis,i,j,resiTopos);
          if(FLAG_EXPAND_HYDROXYL_ROT=TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
        }
      }
    }
  }

  return Success;
}


int StructureGenerateFirstShellProteinRotamersByBBdepRotLib(Structure* pThis,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* shell1sitefile,double shell1cutoff){
  FileReader fr;
  if(!FAILED(FileReaderCreate(&fr, shell1sitefile))){
    char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
      char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1];
      ExtractFirstStringFromSourceString(keyword, buffer);
      if(strcmp(keyword, "mutated_sites_begin") == 0){
        StringArray strings;
        StringArrayCreate(&strings);
        while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
          StringArraySplitString(&strings, buffer, ' ');
          if(strcmp(StringArrayGet(&strings, 0), "mutated_sites_end") == 0){
            break;
          }
          char* chnname = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pThis, chnname);
          int resiIndex = -1;
          ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
          if(resiIndex != -1){
            ChainGetResidue(pChain, resiIndex)->designSiteType = Type_ResidueDesignType_Mutated;
          }
        }
        StringArrayDestroy(&strings);
      }
    }
    FileReaderDestroy(&fr);

    for(int i = 0; i < StructureGetChainCount(pThis); i++){
      Chain *pChainI = StructureGetChain(pThis, i);
      if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
      for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
        Residue* pResidue = ChainGetResidue(pChainI, j);
        if(pResidue->designSiteType == Type_ResidueDesignType_Mutated){
          //printf("generate rotamers for residue %d\n", j);
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pThis,i,j,rotlib,atomParams,resiTopos,pResidue->designSiteType);
          if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamer(pThis,i,j,resiTopos);
          if(FLAG_EXPAND_HYDROXYL_ROT=TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
        }
      }
    }
  }
  else{
    char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    Residue* pSmallMol = NULL;
    StructureFindSmallMol(pThis, &pSmallMol);
    if(pSmallMol==NULL){
      int code = DataNotExistError;
      sprintf(errMsg, "in file %s function %s() line %d, cannot find small molecule", __FILE__, __FUNCTION__, __LINE__);
      TraceError(errMsg, code);
      return code;
    }

    for(int i = 0; i < StructureGetChainCount(pThis); i++){
      Chain *pChainI = StructureGetChain(pThis, i);
      if(strcmp(ResidueGetChainName(pSmallMol),ChainGetName(pChainI))==0) continue;
      if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
      for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
        Residue* pResidue = ChainGetResidue(pChainI, j);
        double minDist=AtomArrayCalcMinDistance(&pSmallMol->atoms, &pResidue->atoms);
        if(minDist>shell1cutoff) continue;
        if(pResidue->designSiteType != Type_ResidueDesignType_Fixed) continue;
        ResidueSetDesignSiteFlag(pResidue, Type_ResidueDesignType_Mutated);
        ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pThis,i,j,rotlib,atomParams,resiTopos,pResidue->designSiteType);
        if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamer(pThis,i,j,resiTopos);
        if(FLAG_EXPAND_HYDROXYL_ROT=TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
      }
    }
  }


  return Success;
}


int StructureGenerateSecondShellProteinRotamersByBBdepRotLib(Structure* pThis,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* shell2sitefile,double shell2cutoff){
  FileReader fr;
  if(!FAILED(FileReaderCreate(&fr, shell2sitefile))){
    char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
      char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1];
      ExtractFirstStringFromSourceString(keyword, buffer);
      if(strcmp(keyword, "rotameric_sites_begin") == 0){
        StringArray strings;
        StringArrayCreate(&strings);
        while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
          StringArraySplitString(&strings, buffer, ' ');
          if(strcmp(StringArrayGet(&strings, 0), "rotameric_sites_end") == 0){
            break;
          }
          char* chnname = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pThis, chnname);
          int resiIndex = -1;
          ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
          if(resiIndex != -1){
            ChainGetResidue(pChain, resiIndex)->designSiteType = Type_ResidueDesignType_Rotameric;
          }
        }
        StringArrayDestroy(&strings);
      }
    }
    FileReaderDestroy(&fr);

    for(int i = 0; i < StructureGetChainCount(pThis); i++){
      Chain *pChainI = StructureGetChain(pThis, i);
      if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
      for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
        Residue* pResidue = ChainGetResidue(pChainI, j);
        if(pResidue->designSiteType == Type_ResidueDesignType_Rotameric){
          //printf("generate rotamers for residue %d\n", j);
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pThis,i,j,rotlib,atomParams,resiTopos,pResidue->designSiteType);
          if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamer(pThis,i,j,resiTopos);
          if(FLAG_EXPAND_HYDROXYL_ROT=TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
        }
      }
    }
  }
  else{
    char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    Residue* pSmallMol = NULL;
    StructureFindSmallMol(pThis, &pSmallMol);
    if(pSmallMol==NULL){
      int code = DataNotExistError;
      sprintf(errMsg, "in file %s function %s() line %d, cannot find small molecule", __FILE__, __FUNCTION__, __LINE__);
      TraceError(errMsg, code);
      return code;
    }

    for(int i = 0; i < StructureGetChainCount(pThis); i++){
      Chain *pChainI = StructureGetChain(pThis, i);
      if(strcmp(ResidueGetChainName(pSmallMol),ChainGetName(pChainI))==0) continue;
      if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
      for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
        Residue* pResidue = ChainGetResidue(pChainI, j);
        double minDist=AtomArrayCalcMinDistance(&pSmallMol->atoms, &pResidue->atoms);
        if(minDist>shell2cutoff) continue;
        if(pResidue->designSiteType != Type_ResidueDesignType_Fixed) continue;
        ResidueSetDesignSiteFlag(pResidue, Type_ResidueDesignType_Rotameric);
        ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pThis,i,j,rotlib,atomParams,resiTopos,pResidue->designSiteType);
        if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamer(pThis,i,j,resiTopos);
        if(FLAG_EXPAND_HYDROXYL_ROT=TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
      }
    }
  }


  return Success;
}


int StructureGenerateBindingSiteRotamersByBBdepRotLib(Structure* pThis,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,double range){
  char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  Residue* pSmallMol = NULL;
  StructureFindSmallMol(pThis, &pSmallMol);
  if(pSmallMol==NULL){
    int code = DataNotExistError;
    sprintf(errMsg, "in file %s function %s() line %d, cannot find small molecule", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, code);
    return code;
  }

  for(int i = 0; i < StructureGetChainCount(pThis); i++){
    Chain *pChainI = StructureGetChain(pThis, i);
    if(strcmp(ResidueGetChainName(pSmallMol),ChainGetName(pChainI))==0) continue;
    if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
    for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if(AtomArrayCalcMinDistance(&pSmallMol->atoms, &pResidue->atoms)>range) continue;
      if(pResidue->designSiteType != Type_ResidueDesignType_Fixed) continue;
      ResidueSetDesignSiteFlag(pResidue, Type_ResidueDesignType_Mutated);
      ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pThis,i,j,rotlib,atomParams,resiTopos,pResidue->designSiteType);
      if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamer(pThis,i,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT=TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
    }
  }

  return Success;
}


int StructureGenerateAllRotamersByBBdepRotLib(Structure* pThis,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* DESCHNS){
  //first, create rotamers for mutated sites
  for(int i = 0; i < StructureGetChainCount(pThis); i++){
    Chain *pChainI = StructureGetChain(pThis, i);
    if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if(strstr(DESCHNS,ChainGetName(pChainI))==0) continue;
    for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if(pResidue->designSiteType != Type_ResidueDesignType_Fixed) continue;
      ResidueSetDesignSiteFlag(pResidue, Type_ResidueDesignType_Mutated);
      ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pThis,i,j,rotlib,atomParams,resiTopos,pResidue->designSiteType);
      if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamer(pThis,i,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT=TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
    }
  }
  //second, create rotamers for rotameric sites => only change conformations
  for(int i = 0; i < StructureGetChainCount(pThis); i++){
    Chain *pChainI = StructureGetChain(pThis, i);
    if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if(strstr(DESCHNS,ChainGetName(pChainI))!=0) continue;
    for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if(pResidue->designSiteType != Type_ResidueDesignType_Fixed) continue;
      Atom* pAtomCA1=ResidueGetAtomByName(pResidue,"CA");
      BOOL interResi=FALSE;
      for(int k=0; k<StructureGetChainCount(pThis); k++){
        if(k==i)continue;
        Chain* pChainK=StructureGetChain(pThis,k);
        if(ChainGetType(pChainK)!=Type_Chain_Protein)continue;
        if(strstr(DESCHNS,ChainGetName(pChainK))==0)continue;
        for(int s=0;s<ChainGetResidueCount(pChainK);s++){
          Residue* pResidueKS=ChainGetResidue(pChainK,s);
          Atom* pAtomCA2=ResidueGetAtomByName(pResidueKS,"CA");
          if(XYZDistance(&pAtomCA2->xyz,&pAtomCA1->xyz)>15.0) continue;
          if(AtomArrayCalcMinDistance(&pResidue->atoms,&pResidueKS->atoms)<PPI_DIST_CUTOFF){
            interResi=TRUE;
            break;
          }
        }
        if(interResi==TRUE) break;
      }
      if(interResi==TRUE){
        ResidueSetDesignSiteFlag(pResidue, Type_ResidueDesignType_Rotameric);
        ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pThis,i,j,rotlib,atomParams,resiTopos,pResidue->designSiteType);
        if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamer(pThis,i,j,resiTopos);
        if(FLAG_EXPAND_HYDROXYL_ROT=TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
      }
    }
  }

  return Success;
}


int StructureGenerateWildtypeRotamersByBBdepRotLib(Structure* pThis,BBdepRotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* DESCHNS){
  //first, create rotamers for mutated sites
  for(int i = 0; i < StructureGetChainCount(pThis); i++){
    Chain *pChainI = StructureGetChain(pThis, i);
    if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if(strstr(DESCHNS,ChainGetName(pChainI))==0) continue;
    for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if(pResidue->designSiteType != Type_ResidueDesignType_Fixed) continue;
      ResidueSetDesignSiteFlag(pResidue, Type_ResidueDesignType_Rotameric);
      ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pThis,i,j,rotlib,atomParams,resiTopos,pResidue->designSiteType);
      if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamerByBBdepRotLib(pThis,i,j,resiTopos,rotlib);
      if(FLAG_EXPAND_HYDROXYL_ROT=TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
    }
  }
  //second, create rotamers for rotameric sites => only change conformations
  for(int i = 0; i < StructureGetChainCount(pThis); i++){
    Chain *pChainI = StructureGetChain(pThis, i);
    if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if(strstr(DESCHNS,ChainGetName(pChainI))!=0) continue;
    for(int j = 0; j < ChainGetResidueCount(pChainI); j++){
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if(pResidue->designSiteType != Type_ResidueDesignType_Fixed) continue;
      Atom* pAtomCA1=ResidueGetAtomByName(pResidue,"CA");
      BOOL interResi=FALSE;
      for(int k=0; k<StructureGetChainCount(pThis); k++){
        if(k==i)continue;
        Chain* pChainK=StructureGetChain(pThis,k);
        if(ChainGetType(pChainK)!=Type_Chain_Protein)continue;
        if(strstr(DESCHNS,ChainGetName(pChainK))==0)continue;
        for(int s=0;s<ChainGetResidueCount(pChainK);s++){
          Residue* pResidueKS=ChainGetResidue(pChainK,s);
          Atom* pAtomCA2=ResidueGetAtomByName(pResidueKS,"CA");
          if(XYZDistance(&pAtomCA2->xyz,&pAtomCA1->xyz)>15.0) continue;
          if(AtomArrayCalcMinDistance(&pResidue->atoms,&pResidueKS->atoms)<PPI_DIST_CUTOFF){
            interResi=TRUE;
            break;
          }
        }
        if(interResi==TRUE) break;
      }
      if(interResi==TRUE){
        ResidueSetDesignSiteFlag(pResidue, Type_ResidueDesignType_Rotameric);
        ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pThis,i,j,rotlib,atomParams,resiTopos,pResidue->designSiteType);
        if(FLAG_ADD_CRYSTAL_ROT==TRUE) ProteinSiteAddCrystalRotamerByBBdepRotLib(pThis,i,j,resiTopos,rotlib);
        if(FLAG_EXPAND_HYDROXYL_ROT=TRUE) ProteinSiteExpandHydroxylRotamers(pThis,i,j,resiTopos);
      }
    }
  }

  return Success;
}



#define CHECK_ROTAMER_CONF_IN_ROTLIB

BOOL ProteinSiteCheckCrystalRotamerInBBdepRotLib(Structure* pThis,int chainIndex,int resiIndex,ResiTopoSet *pResiTopos,BBdepRotamerLib* pBBdepRotLib,double torsionStd){
  char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  Chain* pDestChain = StructureGetChain(pThis, chainIndex);
  Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
  if(pDestChain->type == Type_Chain_Protein){
    //do not add rotamer for residue ala and gly
    if(strcmp(ResidueGetName(pDestResidue), "ALA") == 0 || strcmp(ResidueGetName(pDestResidue), "GLY") == 0){
      return TRUE;
    }

    //step1: calculate the X angles for the current residue
    double DUNBRACK_ENERGY=0.0;
    int torsionCount=0;
    char resname[MAX_LENGTH_RESIDUE_NAME+1];
    strcpy(resname,ResidueGetName(pDestResidue));
    if(!strcmp(resname,"ALA"))     { torsionCount=0;}
    else if(!strcmp(resname,"ARG")){ torsionCount=4;}
    else if(!strcmp(resname,"ASN")){ torsionCount=2;}
    else if(!strcmp(resname,"ASP")){ torsionCount=2;}
    else if(!strcmp(resname,"CYS")){ torsionCount=1;}
    else if(!strcmp(resname,"GLN")){ torsionCount=3;}
    else if(!strcmp(resname,"GLU")){ torsionCount=3;}
    else if(!strcmp(resname,"GLY")){ torsionCount=0;}
    else if(!strcmp(resname,"HSD")){ torsionCount=2;}
    else if(!strcmp(resname,"HSE")){ torsionCount=2;}
    else if(!strcmp(resname,"ILE")){ torsionCount=2;}
    else if(!strcmp(resname,"LEU")){ torsionCount=2;}
    else if(!strcmp(resname,"LYS")){ torsionCount=4;}
    else if(!strcmp(resname,"MET")){ torsionCount=3;}
    else if(!strcmp(resname,"PHE")){ torsionCount=2;}
    else if(!strcmp(resname,"PRO")){ torsionCount=2;}
    else if(!strcmp(resname,"SER")){ torsionCount=1;}
    else if(!strcmp(resname,"THR")){ torsionCount=1;}
    else if(!strcmp(resname,"TRP")){ torsionCount=2;}
    else if(!strcmp(resname,"TYR")){ torsionCount=2;}
    else if(!strcmp(resname,"VAL")){ torsionCount=1;}

    DoubleArray xangles;
    DoubleArrayCreate(&xangles,0);
    ResidueTopology resiTop;
    ResidueTopologyCreate(&resiTop);
    ResiTopoSetGet(pResiTopos,resname,&resiTop);
    for(int torsionIndex=0;torsionIndex<torsionCount;torsionIndex++){
      Type_ProteinAtomOrder desiredAtomBOrder = Type_ProteinAtomOrder_FromInt(torsionIndex);
      Type_ProteinAtomOrder desiredAtomCOrder = Type_ProteinAtomOrder_FromInt(torsionIndex+1);
      CharmmIC icOfCurrentTorsion;
      CharmmICCreate(&icOfCurrentTorsion);
      BOOL icFound=FALSE;
      for(int icIndex=0; icIndex<ResidueTopologyGetCharmmICCount(&resiTop); icIndex++){
        Type_ProteinAtomOrder atomBOrder;
        Type_ProteinAtomOrder atomCOrder;
        ResidueTopologyGetCharmmIC(&resiTop,icIndex,&icOfCurrentTorsion);
        atomBOrder = Type_ProteinAtomOrder_JudgedByAtomName(CharmmICGetAtomB(&icOfCurrentTorsion));
        atomCOrder = Type_ProteinAtomOrder_JudgedByAtomName(CharmmICGetAtomC(&icOfCurrentTorsion));
        if(desiredAtomBOrder==atomBOrder && desiredAtomCOrder==atomCOrder){
          icFound = TRUE;
          break;
        }
      }

      if( !icFound ){
        sprintf(errMsg,"In file %s function %s() line %d, cannot find the CharmmIC for the %dth torsion, residue name: %s",__FILE__, __FUNCTION__, __LINE__,torsionIndex+1, resname);
        TraceError(errMsg,DataNotExistError);
        CharmmICDestroy(&icOfCurrentTorsion);
        return DataNotExistError;
      }

      Atom* pAtomA=ResidueGetAtomByName(pDestResidue,CharmmICGetAtomA(&icOfCurrentTorsion));
      Atom* pAtomB=ResidueGetAtomByName(pDestResidue,CharmmICGetAtomB(&icOfCurrentTorsion));
      Atom* pAtomC=ResidueGetAtomByName(pDestResidue,CharmmICGetAtomC(&icOfCurrentTorsion));
      Atom* pAtomD=ResidueGetAtomByName(pDestResidue,CharmmICGetAtomD(&icOfCurrentTorsion));
      double torsion =GetTorsionAngle(&pAtomA->xyz,&pAtomB->xyz,&pAtomC->xyz,&pAtomD->xyz);
      DoubleArrayAppend(&xangles,torsion);

      CharmmICDestroy(&icOfCurrentTorsion);
    }
    ResidueTopologyDestroy(&resiTop);

    //step2: get the phi&psi angles for the current residue
    int binindex=((int)(pDestResidue->phipsi[0]+180)/10)*36+(int)(pDestResidue->phipsi[1]+180)/10;
    RotLibPhiPsi* pRotLibPhiPsi=&pBBdepRotLib->rotlibphipsis[binindex];
    int rotTypeIndex=-1;
    StringArrayFind(&pRotLibPhiPsi->rotTypes,ResidueGetName(pDestResidue),&rotTypeIndex);
    DoubleArray* pTorsionsArrayForTypeI = pRotLibPhiPsi->torsions[rotTypeIndex];
    DoubleArray* pDeviationsArrayForTypeI = pRotLibPhiPsi->deviations[rotTypeIndex];
    int matchIndex=-1;
    for(int i=0;i<IntArrayGet(&pRotLibPhiPsi->rotamerCounts,rotTypeIndex);i++){
      DoubleArray* pTorsions=&pTorsionsArrayForTypeI[i];
      DoubleArray* pDeviations=&pDeviationsArrayForTypeI[i];
      BOOL match=TRUE;
      for(int j=0;j<DoubleArrayGetLength(&xangles);j++){
        double min=DoubleArrayGet(pTorsions,j)-DegToRad(torsionStd);
        double max=DoubleArrayGet(pTorsions,j)+DegToRad(torsionStd);
        double torsion=DoubleArrayGet(&xangles,j);
        double torsionm2pi=torsion-2*PI;
        double torsionp2pi=torsion+2*PI;
        double torsion2=torsion;
        if((strcmp(resname,"PHE")==0 && j==1)||
          (strcmp(resname,"TYR")==0 && j==1)||
          (strcmp(resname,"ASP")==0 && j==1)||
          strcmp(resname,"GLU")==0 && j==2){
            torsion2=torsion+PI;
            torsion2=torsion>0?torsion-PI:torsion2;
        }
        double torsion2m2pi=torsion2-2*PI;
        double torsion2p2pi=torsion2+2*PI;
        if(!(
          (torsion    <=max && torsion>=min) ||
          (torsionm2pi<=max && torsionm2pi>=min) ||
          (torsionp2pi<=max && torsionp2pi>=min) ||
          (torsion2    <=max && torsion2>=min) ||
          (torsion2m2pi<=max && torsion2m2pi>=min) ||
          (torsion2p2pi<=max && torsion2p2pi>=min)
          )){
            match=FALSE;
            break;
        }
      }
      if(match==TRUE){
        matchIndex=i;
        break;
      }
    }
    DoubleArrayDestroy(&xangles);

    if(matchIndex!=-1){
      return TRUE;
    }
    else{
      return FALSE;
    }

  }

  return Success;
}

BOOL ProteinSiteCheckCrystalRotamerInBBindRotLib(Structure* pThis,int chainIndex,int resiIndex,ResiTopoSet *pResiTopos,BBindRotamerLib* pBBindRotLib,double torsionStd){
  char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  Chain* pDestChain = StructureGetChain(pThis, chainIndex);
  Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
  if(pDestChain->type == Type_Chain_Protein){
    //do not add rotamer for residue ala and gly
    if(strcmp(ResidueGetName(pDestResidue), "ALA") == 0 || strcmp(ResidueGetName(pDestResidue), "GLY") == 0){
      return TRUE;
    }

    //step1: calculate the X angles for the current residue
    double DUNBRACK_ENERGY=0.0;
    int torsionCount=0;
    char resname[MAX_LENGTH_RESIDUE_NAME+1];
    strcpy(resname,ResidueGetName(pDestResidue));
    if(!strcmp(resname,"ALA"))     { torsionCount=0;}
    else if(!strcmp(resname,"ARG")){ torsionCount=4;}
    else if(!strcmp(resname,"ASN")){ torsionCount=2;}
    else if(!strcmp(resname,"ASP")){ torsionCount=2;}
    else if(!strcmp(resname,"CYS")){ torsionCount=1;}
    else if(!strcmp(resname,"GLN")){ torsionCount=3;}
    else if(!strcmp(resname,"GLU")){ torsionCount=3;}
    else if(!strcmp(resname,"GLY")){ torsionCount=0;}
    else if(!strcmp(resname,"HSD")){ torsionCount=2;}
    else if(!strcmp(resname,"HSE")){ torsionCount=2;}
    else if(!strcmp(resname,"ILE")){ torsionCount=2;}
    else if(!strcmp(resname,"LEU")){ torsionCount=2;}
    else if(!strcmp(resname,"LYS")){ torsionCount=4;}
    else if(!strcmp(resname,"MET")){ torsionCount=3;}
    else if(!strcmp(resname,"PHE")){ torsionCount=2;}
    else if(!strcmp(resname,"PRO")){ torsionCount=2;}
    else if(!strcmp(resname,"SER")){ torsionCount=1;}
    else if(!strcmp(resname,"THR")){ torsionCount=1;}
    else if(!strcmp(resname,"TRP")){ torsionCount=2;}
    else if(!strcmp(resname,"TYR")){ torsionCount=2;}
    else if(!strcmp(resname,"VAL")){ torsionCount=1;}

    DoubleArray xangles;
    DoubleArrayCreate(&xangles,0);
    ResidueTopology resiTop;
    ResidueTopologyCreate(&resiTop);
    ResiTopoSetGet(pResiTopos,resname,&resiTop);
    for(int torsionIndex=0;torsionIndex<torsionCount;torsionIndex++){
      Type_ProteinAtomOrder desiredAtomBOrder = Type_ProteinAtomOrder_FromInt(torsionIndex);
      Type_ProteinAtomOrder desiredAtomCOrder = Type_ProteinAtomOrder_FromInt(torsionIndex+1);
      CharmmIC icOfCurrentTorsion;
      CharmmICCreate(&icOfCurrentTorsion);
      BOOL icFound=FALSE;
      for(int icIndex=0; icIndex<ResidueTopologyGetCharmmICCount(&resiTop); icIndex++){
        Type_ProteinAtomOrder atomBOrder;
        Type_ProteinAtomOrder atomCOrder;
        ResidueTopologyGetCharmmIC(&resiTop,icIndex,&icOfCurrentTorsion);
        atomBOrder = Type_ProteinAtomOrder_JudgedByAtomName(CharmmICGetAtomB(&icOfCurrentTorsion));
        atomCOrder = Type_ProteinAtomOrder_JudgedByAtomName(CharmmICGetAtomC(&icOfCurrentTorsion));
        if(desiredAtomBOrder==atomBOrder && desiredAtomCOrder==atomCOrder){
          icFound = TRUE;
          break;
        }
      }

      if( !icFound ){
        sprintf(errMsg,"In file %s function %s() line %d, cannot find the CharmmIC for the %dth torsion, residue name: %s",__FILE__, __FUNCTION__, __LINE__,torsionIndex+1, resname);
        TraceError(errMsg,DataNotExistError);
        CharmmICDestroy(&icOfCurrentTorsion);
        return DataNotExistError;
      }

      Atom* pAtomA=ResidueGetAtomByName(pDestResidue,CharmmICGetAtomA(&icOfCurrentTorsion));
      Atom* pAtomB=ResidueGetAtomByName(pDestResidue,CharmmICGetAtomB(&icOfCurrentTorsion));
      Atom* pAtomC=ResidueGetAtomByName(pDestResidue,CharmmICGetAtomC(&icOfCurrentTorsion));
      Atom* pAtomD=ResidueGetAtomByName(pDestResidue,CharmmICGetAtomD(&icOfCurrentTorsion));
      double torsion =GetTorsionAngle(&pAtomA->xyz,&pAtomB->xyz,&pAtomC->xyz,&pAtomD->xyz);
      DoubleArrayAppend(&xangles,torsion);

      CharmmICDestroy(&icOfCurrentTorsion);
    }
    ResidueTopologyDestroy(&resiTop);


    //step2: check the torsions in rotamer library
    int rotTypeIndex=-1,matchIndex=-1;
    StringArrayFind(&pBBindRotLib->residueTypeNames,ResidueGetName(pDestResidue),&rotTypeIndex);
    int rotCount=IntArrayGet(&pBBindRotLib->rotamerCounts,rotTypeIndex);
    DoubleArray torsions;
    DoubleArrayCreate(&torsions,0);
    for(int i=0;i<rotCount;i++){
      BBindRotamerLibGet(pBBindRotLib,ResidueGetName(pDestResidue),i,&torsions);
      BOOL match=TRUE;
      for(int j=0;j<DoubleArrayGetLength(&xangles);j++){
        double min=DoubleArrayGet(&torsions,j)-DegToRad(torsionStd);
        double max=DoubleArrayGet(&torsions,j)+DegToRad(torsionStd);
        double torsion=DoubleArrayGet(&xangles,j);
        double torsionm2pi=torsion-2*PI;
        double torsionp2pi=torsion+2*PI;
        double torsion2=torsion;
        if((strcmp(resname,"PHE")==0 && j==1)||
          (strcmp(resname,"TYR")==0 && j==1)||
          (strcmp(resname,"ASP")==0 && j==1)||
          strcmp(resname,"GLU")==0 && j==2){
            torsion2=torsion+PI;
            torsion2=torsion>0?torsion-PI:torsion2;
        }
        double torsion2m2pi=torsion2-2*PI;
        double torsion2p2pi=torsion2+2*PI;
        if(!(
          (torsion    <=max && torsion>=min) ||
          (torsionm2pi<=max && torsionm2pi>=min) ||
          (torsionp2pi<=max && torsionp2pi>=min) ||
          (torsion2    <=max && torsion2>=min) ||
          (torsion2m2pi<=max && torsion2m2pi>=min) ||
          (torsion2p2pi<=max && torsion2p2pi>=min)
          )){
            match=FALSE;
            break;
        }
      }
      if(match==TRUE){
        matchIndex=i;
        break;
      }
    }
    DoubleArrayDestroy(&torsions);
    DoubleArrayDestroy(&xangles);

    if(matchIndex!=-1){
      return TRUE;
    }
    else{
      return FALSE;
    }

  }

  return FALSE;
}
