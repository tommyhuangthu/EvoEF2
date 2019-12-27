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
