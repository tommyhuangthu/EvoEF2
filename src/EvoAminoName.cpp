///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include "EvoAminoName.h"

//  FUNCTION    :   Protein::charToAmino()
//  ARGUMENT    :   ---> name: one-letter representation of an amino acid type.
//  RETURN VALUE:   integer representation of that type.
//  COMPLEXITY  :   constant.
CharSeq::AminoName CharSeq::charToAmino(const char name) {
  switch(name){
    case 'A': return ALA; 
    case 'V': return VAL;
    case 'F': return PHE;
    case 'P': return PRO;
    case 'M': return MET;
    case 'I': return ILE;
    case 'L': return LEU;
    case 'D': return ASP;
    case 'E': return GLU;
    case 'K': return LYS;
    case 'R': return ARG;
    case 'S': return SER;
    case 'T': return THR;
    case 'Y': return TYR;
    case 'H': return HIS;
    case 'C': return CYS;
    case 'N': return ASN;
    case 'Q': return GLN;
    case 'W': return TRP;
    case 'G': return GLY;
    default : return ALA;
  }
} 
