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

#include "Utility.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <errno.h>


int StringArrayCreate(StringArray* pThis){
  pThis->stringCount = 0;
  pThis->capacity = 1;
  pThis->strings = (char**)malloc(sizeof(char*)*(pThis->capacity));
  return Success;
}

int StringArrayDestroy(StringArray* pThis){
  for(int i=0;i<pThis->stringCount;i++)  {
    free(pThis->strings[i]);
  }
  free(pThis->strings);
  pThis->stringCount = pThis->capacity = 0;
  pThis->strings = NULL;
  return Success;
}

int StringArrayCopy(StringArray* pThis, StringArray* pOther){
  StringArrayDestroy(pThis);
  pThis->capacity = pOther->capacity;
  pThis->stringCount = pOther->stringCount;
  pThis->strings = (char**)malloc(sizeof(char*)*pThis->capacity);
  for(int i=0;i<pThis->stringCount;i++){
    pThis->strings[i] = (char*)malloc( sizeof(char) * (strlen(pOther->strings[i])+1) );
    strcpy(pThis->strings[i], pOther->strings[i]);
  }
  return Success;
}

int StringArrayGetCount(StringArray* pThis){
  return pThis->stringCount;
}

int StringArraySet(StringArray* pThis, int index, char* srcString){
  if(index<0 || index>=pThis->stringCount){
    return IndexError;
  }
  if(strlen(srcString) > strlen(pThis->strings[index])){
    pThis->strings[index] = (char*)realloc(pThis->strings[index], 
      sizeof(char)*(strlen(srcString)+1));
  }
  strcpy(pThis->strings[index], srcString);
  return Success;
}

char* StringArrayGet(StringArray* pThis, int index){
  if(index<0 || index>=pThis->stringCount){
    return NULL;
  }
  return pThis->strings[index];
}

int StringArrayFind(StringArray* pThis, char* srcString, int* pos){
  int i;
  for(i=0;i<pThis->stringCount;i++){
    if(strcmp(pThis->strings[i], srcString)==0)
      break;
  }
  if(i==pThis->stringCount){
    return DataNotExistError;
  }
  *pos = i;
  return Success;
}

int StringArrayInsert(StringArray* pThis, int pos, char* srcString){
  if(pos<0 || pos>pThis->stringCount){
    return IndexError;
  }
  
  if(pThis->stringCount == pThis->capacity){
    int newCap = pThis->capacity*2;
    pThis->strings = (char**)realloc(pThis->strings, sizeof(char*) * newCap);
    pThis->capacity = newCap;
  }
  
  char* newString = (char*)malloc(sizeof(char)*( strlen(srcString)+1 ));
  strcpy(newString, srcString);

  for(int i=pThis->stringCount; i>pos; i--){
    pThis->strings[i] = pThis->strings[i-1];
  }
  pThis->strings[pos] = newString;
  pThis->stringCount++;
  return Success;
}

int StringArrayRemove(StringArray* pThis, int pos){
  if(pos<0 || pos>=pThis->stringCount){
    return IndexError;
  }
  free(pThis->strings[pos]);
  for(int i=pos;i < pThis->stringCount-1;i++){
    pThis->strings[i] = pThis->strings[i+1];
  }
  pThis->stringCount--;
  return Success;  
}

int StringArrayAppend(StringArray* pThis, char* srcString){
  return StringArrayInsert(pThis, StringArrayGetCount(pThis), srcString);
}

int StringArraySplitString(StringArray* pThis, char* srcStr, char splitter){
  StringArrayDestroy(pThis);
  StringArrayCreate(pThis);
  int length = (int)strlen(srcStr);
  char* buffer = (char*)malloc(sizeof(char)*(length+1));

  int beg = 0;
  while(beg<length){
    int end = beg;
    while(end<length){
      char endChar = srcStr[end];
      if(splitter==endChar || (isspace(splitter)&&isspace(endChar))){
        break;
      }
      else{
        end++;
      }
    }
    if(end>beg){
      strcpy(buffer, srcStr+beg);
      buffer[end-beg] = '\0';
      StringArrayAppend(pThis, buffer);
    }
    beg = end+1;
  }
  free(buffer);
  return Success;
}

int StringArrayShow(StringArray* pThis){
  printf("[");
  for(int i=0;i<pThis->stringCount;i++){
    printf("\"%s\", ", pThis->strings[i]);
  }
  printf("]");
  return Success;
}

int FileReaderCreate(FileReader* pThis, char* path){
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  StringArrayCreate(&pThis->lines);
  FILE* pFile = fopen(path, "r");
  if(pFile==NULL){
    int    result = IOError;
    sprintf(errMsg, "in file %s function %s() line %d, when opening:\n%s", __FILE__, __FUNCTION__, __LINE__, path);
    TraceError(errMsg, result);
    return result;
  }
  char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(fgets(buffer, MAX_LENGTH_ONE_LINE_IN_FILE, pFile)){
    int length = 0;
    while(length< (int)(strlen(buffer)) && buffer[length] != '!')
      length++;
    while(length>0 && buffer[length-1]>0 && isspace(buffer[length-1]))
      length--;
    if(length==0)
      continue;
    buffer[length]='\0';
    for(int i=0;i<length;i++){
      if(buffer[i] < 0){
        int result = FormatError;
        sprintf(errMsg, "in file %s function %s() line %d, only ASCII characters are allowed in non-comment Area : \n%s", __FILE__, __FUNCTION__, __LINE__, buffer);
        TraceError(errMsg, result);
        return result;
      }
    }
    StringArrayAppend(&pThis->lines, buffer);
  }
  pThis->position = 0;

  fclose(pFile);
  return Success;
}

int FileReaderDestroy(FileReader* pThis){
  StringArrayDestroy(&pThis->lines);
  return Success;
}

int FileReaderGetLineCount(FileReader* pThis){
  return StringArrayGetCount(&pThis->lines);
}

int FileReaderGetLine(FileReader* pThis, int index, char* dest){
  char* line = StringArrayGet(&pThis->lines, index);
  if(line==NULL){
    return IndexError;
  }
  else{
    strcpy(dest, line);
    return Success;
  }
}

int FileReaderGetCurrentPos(FileReader* pThis){
  return pThis->position;
}

int FileReaderSetCurrentPos(FileReader* pThis, int index){
  if(index<0 || index>=FileReaderGetLineCount(pThis)){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = IndexError;
    sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  pThis->position = index;
  return Success;
}

int FileReaderGetNextLine(FileReader* pThis, char* dest){
  int result = FileReaderGetLine(pThis, pThis->position, dest);
  if(FAILED(result)){
    return result;
  }
  else{
    (pThis->position)++;
     return Success;
  }
}

BOOL FileReaderEndOfFile(FileReader* pThis){
  if(pThis->position == FileReaderGetLineCount(pThis))
    return TRUE;
  else
    return FALSE;
}

Type_CoordinateFile FileReaderRecognizeCoordinateFileType(FileReader* pThis){
  for(int i=0;i<FileReaderGetLineCount(pThis);i++){
    char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    FileReaderGetLine(pThis, i, line);
    if(strlen(line)>9){
      line[9] = '\0';
      if(strcmp(line, "@<TRIPOS>")==0)
        return Type_CoordinateFile_MOL2;
    }
    if(strlen(line)>6){
      line[6] = '\0';
      if(strcmp(line, "ATOM  ")==0 || strcmp(line, "HETATM")==0 )
        return Type_CoordinateFile_PDB;
    }
  }
  return Type_CoordinateFile_Unrecognized;
}

int IntArrayCreate(IntArray* pThis, int length){
  if(length<0)
    return IndexError;
  pThis->content = (int*)calloc(length, sizeof(int));
  pThis->length = length;
  pThis->capacity = length;
  return Success;
}

int IntArrayDestroy(IntArray* pThis){
  free(pThis->content);
  pThis->content = NULL;
  pThis->length = 0;
  pThis->capacity = 0;
  return Success;
}

int IntArrayCopy(IntArray* pThis, IntArray* pOther){
  IntArrayDestroy(pThis);
  pThis->content = (int*)calloc(pOther->length, sizeof(int));
  pThis->length = pOther->length;
  pThis->capacity = pThis->length;
  IntArraySetAll(pThis, IntArrayGetAll(pOther));
  return Success;
}

int IntArrayResize(IntArray* pThis, int newLength){
  if(newLength<0)
    return IndexError;
  pThis->content = (int*)realloc(pThis->content, sizeof(int)*newLength);
  for(int i=pThis->length;i<newLength;i++){
    pThis->content[i] = 0;
  }
  pThis->length = newLength;
  pThis->capacity = pThis->length;
  return Success;
}

int IntArrayGetLength(IntArray* pThis){
  return pThis->length;
}

int IntArrayGet(IntArray* pThis, int index){
  if(index<0 || index>=pThis->length){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = IndexError;
    sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  return pThis->content[index];
}
int IntArraySet(IntArray* pThis, int index, int newValue){
  if(index<0 || index>=pThis->length)
    return IndexError;
  pThis->content[index] = newValue;
  return Success;
}

int* IntArrayGetAll(IntArray* pThis){
  return pThis->content;
}

int IntArraySetAll(IntArray* pThis, int* pNewContent){
  if(pNewContent==NULL)
    return ValueError;
  for(int i=0;i<pThis->length;i++){
    pThis->content[i] = pNewContent[i];
  }
  return Success;
}

int IntArrayInsert(IntArray* pThis, int index, int newValue){
  if(index<0 || index>pThis->length)
    return IndexError;
  if(pThis->capacity == pThis->length){
    pThis->capacity = pThis->capacity*2 + 1;
    pThis->content = (int*)realloc(pThis->content, sizeof(int)*pThis->capacity);
  }
  for(int i=pThis->length;i>index;i--){
    pThis->content[i] = pThis->content[i-1];
  }
  pThis->content[index] = newValue;

  (pThis->length)++;
  return Success;
}
  
int IntArrayRemove(IntArray* pThis, int index){
  if(index<0 || index>pThis->length){
    return IndexError;
  }
  for(int i=index;i<pThis->length-1;i++){
    pThis->content[i] = pThis->content[i+1];
  }
  (pThis->length)--;
  return Success;
}

int IntArrayAppend(IntArray* pThis, int newValue){
  return IntArrayInsert(pThis, IntArrayGetLength(pThis), newValue);
}

int IntArrayShow(IntArray* pThis){
  printf("[");
  for(int i=0;i<pThis->length;i++){
    printf("%d, ", pThis->content[i]);
  }
  printf("]");
  return Success;
}

int IntArrayFind(IntArray *pThis, int num){
  for(int i = 0;i < pThis->length; i++){
    if(pThis->content[i] == num){
      return i;
    }
  }
  return -1;
}


int DoubleArrayCreate(DoubleArray* pThis, int length){
  if(length<0)
    return IndexError;

  pThis->content = (double*)calloc(length, sizeof(double));
  pThis->length = length;
  pThis->capacity = length;
  for(int i=0;i<pThis->length;i++){
    pThis->content[i] = 0.0;
  }
  return Success;
}

int DoubleArrayDestroy(DoubleArray* pThis){
  free(pThis->content);
  pThis->content = NULL;
  pThis->length = 0;
  pThis->capacity = 0;
  return Success;
}

int DoubleArrayCopy(DoubleArray* pThis, DoubleArray* pOther){
  if(pThis->capacity >= pOther->length){
    DoubleArraySetAll(pThis, DoubleArrayGetAll(pOther));
    return Success;
  }
  DoubleArrayDestroy(pThis);
  pThis->content = (double*)calloc(pOther->length, sizeof(double));
  pThis->length = pOther->length;
  pThis->capacity = pThis->length;
  DoubleArraySetAll(pThis, DoubleArrayGetAll(pOther));
  return Success;
}

int DoubleArrayResize(DoubleArray* pThis, int newLength){
  if(newLength<0)
    return IndexError;

  pThis->content = (double*)realloc(pThis->content, sizeof(double)*newLength);
  for(int i=pThis->length;i<newLength;i++){
    pThis->content[i] = 0.0;
  }
  pThis->length = newLength;
  pThis->capacity = pThis->length;
  return Success;
}

int DoubleArrayGetLength(DoubleArray* pThis){
  return pThis->length;
}

double DoubleArrayGet(DoubleArray* pThis, int index){
  if(index<0 || index>=pThis->length){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = IndexError;
    sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  return pThis->content[index];
}

int DoubleArraySet(DoubleArray* pThis, int index, double newValue){
  if(index<0 || index>=pThis->length)
    return IndexError;
  pThis->content[index] = newValue;
  return Success;
}

double* DoubleArrayGetAll(DoubleArray* pThis){
  return pThis->content;
}

int DoubleArraySetAll(DoubleArray* pThis, double* pNewContent){
  if(pNewContent==NULL)
    return ValueError;
  for(int i=0;i<pThis->length;i++){
    pThis->content[i] = pNewContent[i];
  }
  return Success;
}

int DoubleArrayInsert(DoubleArray* pThis, int index, double newValue){
  if(index<0 || index>pThis->length)
    return IndexError;
  if(pThis->capacity == pThis->length){
    pThis->capacity = pThis->capacity*2 + 1;
    pThis->content = (double*)realloc(pThis->content, sizeof(double)*pThis->capacity);
  }
  for(int i=pThis->length;i>index;i--){
    pThis->content[i] = pThis->content[i-1];
  }
  pThis->content[index] = newValue;
  (pThis->length)++;
  return Success;
}

int DoubleArrayRemove(DoubleArray* pThis, int index){
  if(index<0 || index>pThis->length)
    return IndexError;

  for(int i=index;i<pThis->length-1;i++){
    pThis->content[i] = pThis->content[i+1];
  }
  (pThis->length)--;
  return Success;
}

int DoubleArrayAppend(DoubleArray* pThis, double newValue){
  return DoubleArrayInsert(pThis, DoubleArrayGetLength(pThis), newValue);
}


double DoubleArrayInnerProduct(DoubleArray *pThis, DoubleArray *pOther){
  double product = 0.0;
  if(pThis->length == pOther->length){
    for(int i = 0; i < pThis->length; i++){
      product += pThis->content[i]*pOther->content[i];
    }
    return product;
  }
  printf("in file %s function %s() line %, two vectors have different dimensions\n",__FILE__, __FUNCTION__, __LINE__);
  return ValueError;
}

int DoubleArrayScale(DoubleArray *pThis, double scale){
  for(int i = 0; i < pThis->length; i++){
    pThis->content[i] *= scale;
  }
  return Success;
}

double DoubleArrayNorm(DoubleArray *pThis){
  double norm = 0.0;
  for(int i = 0; i < pThis->length; i++){
    norm += pThis->content[i]*pThis->content[i];
  }
  return sqrt(norm);
}

int  DoubleArrayMinus(DoubleArray *pThis, DoubleArray *pOther){
  if(pThis->length != pOther->length){
    printf("in file %s function %s() line %, two vectors have different dimensions\n",__FILE__, __FUNCTION__, __LINE__);
    return ValueError;
  }
  for(int i = 0; i < pThis->length; i++){
    pThis->content[i] -= pOther->content[i];
  }
  return Success;
}

int ExtractTargetStringFromSourceString(char* dest, char* src, int start, int length){
  // the function strncpy() from <string.h> has the following format;
  // char *strncpy(char *strDest, char *strSource, size_t count);
  // the function copies characters with a maximum number of 'count' from string 'strSource' to string 'strDest', 
  // if 'count' <= strlen(strSource), '\0' will not be added to the end of 'strDest'; else added.
  // given the parameters 'beg' and 'length';
  // the following operation:
  // strncpy(strDest, strSource+beg, length);
  // strDest[length] = '\0';
  // equals to
  // copying a string with a length of 'length' from the 'beg' character of 'strSource' to 'strDest'.
  char* temp = (char*)malloc(sizeof(char)*(length+1));
  strncpy(temp, src+start, length);
  temp[length]='\0';
  // discard redundant space characters;
  sscanf(temp, "%s", dest);
  free(temp);
  return Success;
}

int ExtractFirstStringFromSourceString(char* dest, char* src){
  if(sscanf(src, "%s", dest)==EOF)
    return ValueError;
  // the number of space characters before the first string;
  int from = 0;
  while(isspace(src[from]))
    from++;
  // add the length of the first string;
  from += (int)strlen(dest);
  // add the number of space characters after the first string;
  while(isspace(src[from]))
    from++;
  int to = strlen(src);
  // the left string is still stored in string 'src';
  for(int i=from;i<=to;i++){
    src[i-from] = src[i];
  }
  return Success;
}

int Model(int i, FILE* pFile){
  if(pFile==NULL){
    pFile = stdout;
  }
  fprintf(pFile, "MODEL     %d\n", i);
  return Success;
}

int EndModel(FILE* pFile){
  if(pFile==NULL){
    pFile = stdout;
  }
  fprintf(pFile, "ENDMDL\n");
  return Success;
}

int ShowProgress(int width, double percentage){
  if(width<0){
    width=0;
  }
  int complete = (int)(width*percentage/100.0);
  printf("|");
  for(int i=0;i<complete;i++){
    printf(">");
  }
  for(int i=complete;i<width;i++){
    printf("=");
  }
  printf("|");
  printf("%5.2f%% ", percentage);
  return Success;
}

int SpentTimeShow(time_t ts, time_t te){
  printf("Time spent: %f\n", (double)(te-ts)/CLOCKS_PER_SEC);
  return Success;
}

int OneLetterAAToThreeLetterAA(char name1, char name3[MAX_LENGTH_RESIDUE_NAME]){
  switch(name1){
    case 'A':
      strcpy(name3, "ALA");
      break;
    case 'C':
      strcpy(name3, "CYS");
      break;
    case 'D':
      strcpy(name3, "ASP");
      break;
    case 'E':
      strcpy(name3, "GLU");
      break;
    case 'F':
      strcpy(name3, "PHE");
      break;
    case 'G':
      strcpy(name3, "GLY");
      break;
    case 'H':
      strcpy(name3, "HSD");
      break;
    case 'I':
      strcpy(name3, "ILE");
      break;
    case 'K':
      strcpy(name3, "LYS");
      break;
    case 'L':
      strcpy(name3, "LEU");
      break;
    case 'M':
      strcpy(name3, "MET");
      break;
    case 'N':
      strcpy(name3, "ASN");
      break;
    case 'P':
      strcpy(name3, "PRO");
      break;
    case 'Q':
      strcpy(name3, "GLN");
      break;
    case 'R':
      strcpy(name3, "ARG");
      break;
    case 'S':
      strcpy(name3, "SER");
      break;
    case 'T':
      strcpy(name3, "THR");
      break;
    case 'V':
      strcpy(name3, "VAL");
      break;
    case 'W':
      strcpy(name3, "TRP");
      break;
    case 'Y':
      strcpy(name3, "TYR");
      break;
    default:
      strcpy(name3, "XXX");
      break;
  }
  return Success;
}


char ThreeLetterAAToOneLetterAA(char name3[MAX_LENGTH_RESIDUE_NAME]){
  if(!strcmp(name3,"ALA")) return 'A';
  else if(!strcmp(name3,"CYS")) return 'C';
  else if(!strcmp(name3,"ASP")) return 'D';
  else if(!strcmp(name3,"GLU")) return 'E';
  else if(!strcmp(name3,"PHE")) return 'F';
  else if(!strcmp(name3,"GLY")) return 'G';
  else if(!strcmp(name3,"HSD")) return 'H';
  else if(!strcmp(name3,"HSE")) return 'H';
  else if(!strcmp(name3,"HSP")) return 'H';
  else if(!strcmp(name3,"HIS")) return 'H';
  else if(!strcmp(name3,"ILE")) return 'I';
  else if(!strcmp(name3,"LYS")) return 'K';
  else if(!strcmp(name3,"LEU")) return 'L';
  else if(!strcmp(name3,"MET")) return 'M';
  else if(!strcmp(name3,"ASN")) return 'N';
  else if(!strcmp(name3,"PRO")) return 'P';
  else if(!strcmp(name3,"GLN")) return 'Q';
  else if(!strcmp(name3,"ARG")) return 'R';
  else if(!strcmp(name3,"SER")) return 'S';
  else if(!strcmp(name3,"THR")) return 'T';
  else if(!strcmp(name3,"VAL")) return 'V';
  else if(!strcmp(name3,"TRP")) return 'W';
  else if(!strcmp(name3,"TYR")) return 'Y';
  else if(!strcmp(name3,"XXX")) return 'X'; //ligand
  else return 'X'; //ligand
}


BOOL AminoAcidIsPolar(char* resiName){
  char oneLetterName = ThreeLetterAAToOneLetterAA(resiName);
  switch(oneLetterName){
    case 'A':
      return FALSE;
    case 'C':
      return FALSE;
    case 'D':
      return TRUE;
    case 'E':
      return TRUE;
    case 'F':
      return FALSE;
    case 'G':
      return FALSE;
    case 'H':
      return TRUE;
    case 'I':
      return FALSE;
    case 'K':
      return TRUE;
    case 'L':
      return FALSE;
    case 'M':
      return FALSE;
    case 'N':
      return TRUE;
    case 'P':
      return FALSE;
    case 'Q':
      return TRUE;
    case 'R':
      return TRUE;
    case 'S':
      return TRUE;
    case 'T':
      return TRUE;
    case 'V':
      return FALSE;
    case 'W':
      return FALSE;
    case 'Y':
      return FALSE;
    default:
      return FALSE;
  }
  return FALSE;
}

int AminoAcidGetIndex(char name1){
  switch(name1){
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'D':
      return 2;
    case 'E':
      return 3;
    case 'F':
      return 4;
    case 'G':
      return 5;
    case 'H':
      return 6;
    case 'I':
      return 7;
    case 'K':
      return 8;
    case 'L':
      return 9;
    case 'M':
      return 10;
    case 'N':
      return 11;
    case 'P':
      return 12;
    case 'Q':
      return 13;
    case 'R':
      return 14;
    case 'S':
      return 15;
    case 'T':
      return 16;
    case 'V':
      return 17;
    case 'W':
      return 18;
    case 'Y':
      return 19;
    default:
      return 20;
  }
}


int ThreeLetterAAGetIndex(char name3[MAX_LENGTH_RESIDUE_NAME]){
  if(!strcmp(name3,"ALA")) return 0;
  else if(!strcmp(name3,"CYS")) return 1;
  else if(!strcmp(name3,"ASP")) return 2;
  else if(!strcmp(name3,"GLU")) return 3;
  else if(!strcmp(name3,"PHE")) return 4;
  else if(!strcmp(name3,"GLY")) return 5;
  else if(!strcmp(name3,"HSD")) return 6;
  else if(!strcmp(name3,"HSE")) return 6;
  else if(!strcmp(name3,"HSP")) return 6;
  else if(!strcmp(name3,"HIS")) return 6;
  else if(!strcmp(name3,"ILE")) return 7;
  else if(!strcmp(name3,"LYS")) return 8;
  else if(!strcmp(name3,"LEU")) return 9;
  else if(!strcmp(name3,"MET")) return 10;
  else if(!strcmp(name3,"ASN")) return 11;
  else if(!strcmp(name3,"PRO")) return 12;
  else if(!strcmp(name3,"GLN")) return 13;
  else if(!strcmp(name3,"ARG")) return 14;
  else if(!strcmp(name3,"SER")) return 15;
  else if(!strcmp(name3,"THR")) return 16;
  else if(!strcmp(name3,"VAL")) return 17;
  else if(!strcmp(name3,"TRP")) return 18;
  else if(!strcmp(name3,"TYR")) return 19;
  else if(!strcmp(name3,"XXX")) return 20; //ligand
  else return 20; //ligand
}