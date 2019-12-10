///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#ifndef _GETPHIPSI_H_
#define _GETPHIPSI_H_

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "EvoUtility.h"

namespace PhiPsiPrediction {
  class getPhiPsi {
    //char pwd[500];
    public:
      getPhiPsi(SequenceData dsInfo,char evolutionpath[500]);
      void systemCall(char* evolutionpath);
  };
}


#endif
