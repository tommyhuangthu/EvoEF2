///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#ifndef EVOLUTION_H
#define EVOLUTION_H

//the first parameter is the path of the current executable program used to compute evolution score
//this path is very important, because it is used to determine the location of feature parameters required by the program
float EvolutionScore(char* program_path, char* prffile, char* ssfile, char* safile, char* phipsifile, char* fasseqfile);
float EvolutionScore2(char* programpath, char* prffile, char* ssfile, char* safile, char* phipsifile, char* fasseq);
float EvolutionScore3(char* programpath, char* prffile, char* fasseq);

#endif
