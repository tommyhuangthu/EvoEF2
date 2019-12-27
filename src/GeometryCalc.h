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

#ifndef GEO_CALC_H
#define GEO_CALC_H

#include "ErrorTracker.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define PI 3.1415926535898
#define MIN_ZERO_TOLERANCE 1e-7

typedef struct _XYZ{
  double X, Y, Z;
} XYZ;

int XYZShow(XYZ* pThis);
int XYZScale(XYZ* pThis, double ratio);
int XYZAdd(XYZ* pThis, XYZ* pOther);
int XYZMinus(XYZ* pThis, XYZ* pOther);
XYZ XYZSum(XYZ* pThis, XYZ* pOther);
XYZ XYZDifference(XYZ* pThis, XYZ* pOther);
double XYZNormalization(XYZ* pThis);
double XYZDistance(XYZ* pThis, XYZ* pOther);
double XYZDotProduct(XYZ* pThis, XYZ* pOther);
XYZ XYZCrossProduct(XYZ* pThis, XYZ* pOther);
double XYZAngle(XYZ* pThis, XYZ* pOther);
XYZ XYZRotateAround(XYZ* pThis, XYZ* axisFrom, XYZ* axisTo, double angle);
int XYZRandomlyGenerate(XYZ* pThis, double range);


typedef struct _XYZArray{
  XYZ* xyzs;
  int xyzCount;
} XYZArray;

int XYZArrayCreate(XYZArray* pThis, int length);
int XYZArrayDestroy(XYZArray* pThis);
int XYZArrayCopy(XYZArray* pThis, XYZArray* pOther);
int XYZArrayResize(XYZArray* pThis, int newLength);
int XYZArrayGetLength(XYZArray* pThis);
XYZ* XYZArrayGet(XYZArray* pThis, int index);
int XYZArraySet(XYZArray* pThis, int index, XYZ* newXYZ);
XYZ* XYZArrayGetAll(XYZArray* pThis);
int XYZArrayShow(XYZArray* pThis);
double XYZArrayRMSD(XYZArray* pThis, XYZArray* pOther);


typedef struct _FourXYZsGroup{
  XYZ atomA, atomB, atomC, atomD;
} FourXYZsGroup;

int FourXYZsGroupCreate(FourXYZsGroup* pThis, XYZ* pAtomA, XYZ* pAtomB, XYZ* pAtomC, XYZ* pAtomD);
int FourXYZsGroupGetTorsionAngle(FourXYZsGroup* pThis,double* angle);
int FourXYZsGroupGetFourthAtom(FourXYZsGroup* pThis, double* icParam, XYZ* pAtomD);
int FourXYZsGroupGetFourthAtomNew(FourXYZsGroup* pThis, double* icParam, XYZ* pAtomD);
int FourXYZsGroupGetICParam(FourXYZsGroup* pThis, int torsionProperFlag, double* icParam);
int GetFourthAtom(XYZ* pAtomA, XYZ* pAtomB, XYZ* pAtomC, double* icParam, XYZ* pAtomD);
double GetTorsionAngle(XYZ* pAtomA, XYZ* pAtomB, XYZ* pAtomC, XYZ* pAtomD);


double RadToDeg(double rad);
double DegToRad(double degree);
int Matrix4By4TimesVector(double result[4], double matrix[4][4], double v[4]);
double SafeArccos(double cosValue);
double RandomDouble(double low, double high);
BOOL RadInRange(double value, double low, double high);

#endif // GEO_CALC_H
