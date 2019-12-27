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

#ifndef ERROR_HANDLING_H
#define ERROR_HANDLING_H

#define _CRT_SECURE_NO_WARNINGS

#define MAX_LENGTH_ERR_MSG 1024
#define BOOL               char
#define TRUE               1
#define FALSE              0

typedef enum _ErrorCode{
  Success           =   0, 
  Warning           =   1, 
  NameError         =  -1, 
  IOError           =  -2, 
  FormatError       =  -3, 
  IndexError        =  -4, 
  ValueError        =  -5, 
  ZeroDivisonError  =  -6, 
  DataNotExistError =  -7, 
  AssertionError    =  -8,
  Exception         = -10, 
} ErrorCode;

BOOL FAILED(int errorCode);
int TraceError(char* errMsg, int errorCode);

#endif // ERROR_HANDLING_H
