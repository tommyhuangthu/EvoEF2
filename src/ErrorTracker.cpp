///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>
//////////////////////////////////////////////////////////////////////////////////////

#include "ErrorTracker.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

BOOL FAILED(int errorCode){
  if(errorCode == Success || errorCode == Warning){
    return FALSE;
  }
  return TRUE;
}

int TraceError(char* errMsg, int errorCode){
  char defaultMsg[MAX_LENGTH_ERR_MSG+1];

  if(errorCode == Success){
    return errorCode;
  }
  else if (errorCode == Warning){
    strcpy(defaultMsg, "Warning:");
    printf("%s %s\n", defaultMsg, errMsg);
    return errorCode;
  }

  switch(errorCode){
    case IOError:
      strcpy(defaultMsg, "IOError:"); break;
    case FormatError:
      strcpy(defaultMsg, "FormatError:"); break;
    case IndexError:
      strcpy(defaultMsg, "IndexError:"); break;
    case ValueError:
      strcpy(defaultMsg, "ValueError:"); break;
    case ZeroDivisonError:
      strcpy(defaultMsg, "ZeroDivisonError:"); break;
    case DataNotExistError:
      strcpy(defaultMsg, "DataNotExistError:"); break;
    case NameError:
      strcpy(defaultMsg, "NameError:"); break;
    default:
      strcpy(defaultMsg, "OtherError:");
  }

  printf("%s %s\n", defaultMsg, errMsg);
  return errorCode;
}
