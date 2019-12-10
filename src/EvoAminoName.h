///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

/******************************************************************************/
/*                                                                            */
/*                                AminoName.h                                 */
/*                                                                            */
/******************************************************************************/

#ifndef AMINO_NAME_H
#define AMINO_NAME_H

namespace CharSeq {
  const int AMINOS = 20;
  enum AminoName {GLY, PRO, ASP, GLU, LYS, ARG, HIS, SER, THR, ASN,
                  GLN, ALA, MET, TYR, TRP, VAL, ILE, LEU, PHE, CYS};
  // conversion between one-letter and integer representations of aa types
  AminoName charToAmino(const char name);
}

#endif
