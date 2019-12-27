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

#pragma warning(disable:4305)
#pragma warning(disable:4244)
#include "GeometryCalc.h"
#include "ErrorTracker.h"
#include <time.h>

int XYZShow(XYZ* pThis){
  printf("[%.2f, %.2f, %.2f]", pThis->X, pThis->Y, pThis->Z);
  return Success;
}

int XYZScale(XYZ* pThis, double ratio){
  pThis->X *= ratio;
  pThis->Y *= ratio;
  pThis->Z *= ratio;
  return Success;
}

int XYZAdd(XYZ* pThis, XYZ* pOther){
  pThis->X += pOther->X;
  pThis->Y += pOther->Y;
  pThis->Z += pOther->Z;
  return Success;
}

int XYZMinus(XYZ* pThis, XYZ* pOther){
  pThis->X -= pOther->X;
  pThis->Y -= pOther->Y;
  pThis->Z -= pOther->Z;
  return Success;
}

XYZ XYZSum(XYZ* pThis, XYZ* pOther){
  XYZ xyz;
  xyz.X = pThis->X + pOther->X;
  xyz.Y = pThis->Y + pOther->Y;
  xyz.Z = pThis->Z + pOther->Z;
  return xyz;
}

XYZ XYZDifference(XYZ* pThis, XYZ* pOther){
  XYZ xyz;
  xyz.X = pOther->X - pThis->X;
  xyz.Y = pOther->Y - pThis->Y;
  xyz.Z = pOther->Z - pThis->Z;
  return xyz;
}

double XYZNormalization(XYZ* pThis){
  return sqrt(XYZDotProduct(pThis, pThis));
}

double XYZDistance(XYZ* pThis, XYZ* pOther){
  double distX = pThis->X - pOther->X;
  double distY = pThis->Y - pOther->Y;
  double distZ = pThis->Z - pOther->Z;
  return sqrt(distX * distX + distY * distY + distZ * distZ);
}

double XYZDotProduct(XYZ* pThis, XYZ* pOther){
  return pThis->X*pOther->X + pThis->Y*pOther->Y + pThis->Z*pOther->Z;
}

XYZ XYZCrossProduct(XYZ* pThis, XYZ* pOther){
  XYZ product;
  product.X = pThis->Y*pOther->Z - pOther->Y*pThis->Z;
  product.Y = -pThis->X*pOther->Z + pOther->X*pThis->Z;
  product.Z = pThis->X*pOther->Y - pOther->X*pThis->Y;
  return product;
}

double XYZAngle(XYZ* pThis, XYZ* pOther){
  double cosValue;
  double this2 = pThis->X * pThis->X + pThis->Y * pThis->Y + pThis->Z * pThis->Z;
  double other2 = pOther->X * pOther->X + pOther->Y * pOther->Y + pOther->Z * pOther->Z;
  double norm = sqrt(this2 * other2);
  if(norm < MIN_ZERO_TOLERANCE)  {
    printf("pThis: [%f, %f, %f]\n",pThis->X,pThis->Y,pThis->Z);
    printf("pOther: [%f, %f, %f]\n",pOther->X,pOther->Y,pOther->Z);
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ZeroDivisonError;
    sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return 1000.0;
  }
  cosValue =  XYZDotProduct(pThis, pOther) / norm;
  return SafeArccos(cosValue);
}

XYZ XYZRotateAround(XYZ* pThis, XYZ* axisFrom, XYZ* axisTo, double angle){
  double s = sin(angle);
  double c = cos(angle);
  XYZ result = *pThis;
  XYZ n = XYZDifference(axisFrom, axisTo);
  double normOfAxis = XYZNormalization(&n);
  n.X /= normOfAxis;
  n.Y /= normOfAxis;
  n.Z /= normOfAxis;
  XYZMinus(&result, axisFrom);

  double A[3][3];
  A[0][0] = n.X*n.X + (1-n.X*n.X)*c;
  A[0][1] = n.X*n.Y*(1-c) - n.Z*s;
  A[0][2] = n.X*n.Z*(1-c) + n.Y*s;
  A[1][0] = n.X*n.Y*(1-c) + n.Z*s;
  A[1][1] = n.Y*n.Y + (1-n.Y*n.Y)*c;
  A[1][2] = n.Y*n.Z*(1-c) - n.X*s;
  A[2][0] = n.X*n.Z*(1-c) - n.Y*s;
  A[2][1] = n.Y*n.Z*(1-c) + n.X*s;
  A[2][2] = n.Z*n.Z + (1-n.Z*n.Z)*c;

  double beforeMultiplication[3];
  beforeMultiplication[0] = result.X;
  beforeMultiplication[1] = result.Y;
  beforeMultiplication[2] = result.Z;
  double afterMultiplication[3];
  for(int i=0;i<3;i++){
    afterMultiplication[i] = 0.0;
    for(int j=0;j<3;j++){
      afterMultiplication[i] += beforeMultiplication[j] * A[i][j];
    }
  }
  result.X = afterMultiplication[0];
  result.Y = afterMultiplication[1];
  result.Z = afterMultiplication[2];
  XYZAdd(&result, axisFrom);
  return result;
}

int XYZRandomlyGenerate(XYZ* pThis, double range){
  pThis->X = RandomDouble(-range, range);
  pThis->Y = RandomDouble(-range, range);
  pThis->Z = RandomDouble(-range, range);
  return Success;
}


int XYZArrayCreate(XYZArray* pThis, int length){
  pThis->xyzCount = length;
  pThis->xyzs = (XYZ*)calloc(length, sizeof(XYZ));
  return Success;
}

int XYZArrayDestroy(XYZArray* pThis){
  free(pThis->xyzs);
  pThis->xyzs = NULL;
  pThis->xyzCount = 0;
  return Success;
}

int XYZArrayCopy(XYZArray* pThis, XYZArray* pOther){
  XYZArrayDestroy(pThis);
  XYZArrayCreate(pThis, pOther->xyzCount);
  for(int i=0;i<pThis->xyzCount;i++){
    pThis->xyzs[i] = pOther->xyzs[i];
  }
  return Success;
}
int XYZArrayResize(XYZArray* pThis, int newLength){
  pThis->xyzs = (XYZ*)realloc(pThis->xyzs, sizeof(XYZ)*newLength);
  for(int i=pThis->xyzCount; i<newLength; i++){
    pThis->xyzs[i].X = pThis->xyzs[i].Y = pThis->xyzs[i].Z = 0;
  }
  pThis->xyzCount = newLength;
  return Success;
}

int XYZArrayGetLength(XYZArray* pThis){
  return pThis->xyzCount;
}

XYZ* XYZArrayGet(XYZArray* pThis, int index){
  if(index<0 || index>=pThis->xyzCount) return NULL;
  return &pThis->xyzs[index];
}

int XYZArraySet(XYZArray* pThis, int index, XYZ* newXYZ){
  if(index<0 || index>=pThis->xyzCount) return IndexError;
  pThis->xyzs[index] = *newXYZ;
  return Success;
}

XYZ* XYZArrayGetAll(XYZArray* pThis){
  return pThis->xyzs;
}

int XYZArrayShow(XYZArray* pThis){
  for(int i=0;i<pThis->xyzCount;i++){
    printf("%d : ", i);
    XYZShow(&pThis->xyzs[i]);
    printf("\n");
  }
  return Success;
}
double XYZArrayRMSD(XYZArray* pThis, XYZArray* pOther){
  double sum = 0.0;
  for(int i=0;i<pThis->xyzCount;i++){
    double distance = XYZDistance(&pThis->xyzs[i], &pOther->xyzs[i]);
    sum += distance * distance;
  }
  return sqrt(sum/ pThis->xyzCount);
}

int FourXYZsGroupCreate(FourXYZsGroup* pThis, XYZ* pAtomA, XYZ* pAtomB, XYZ* pAtomC, XYZ* pAtomD){
  if(pAtomA==NULL || pAtomB==NULL || pAtomC==NULL) return ValueError;
  pThis->atomA = *pAtomA;
  pThis->atomB = *pAtomB;
  pThis->atomC = *pAtomC;
  if(pAtomD != NULL){ pThis->atomD = *pAtomD; }
  return Success;
}

int FourXYZsGroupGetTorsionAngle(FourXYZsGroup* pThis,double* angle){
  XYZ rAB = XYZDifference(&pThis->atomA, &pThis->atomB);
  XYZ rBC = XYZDifference(&pThis->atomB, &pThis->atomC);
  XYZ rCD = XYZDifference(&pThis->atomC, &pThis->atomD);
  XYZ rABXrBC = XYZCrossProduct(&rAB, &rBC);
  XYZ rBCXrCD = XYZCrossProduct(&rBC, &rCD);
  XYZ rCDXrAB = XYZCrossProduct(&rCD, &rAB);

  if(XYZNormalization(&rABXrBC)<MIN_ZERO_TOLERANCE || XYZNormalization(&rBCXrCD)<MIN_ZERO_TOLERANCE ){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ZeroDivisonError;
    sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  double cosValue = XYZDotProduct(&rABXrBC, &rBCXrCD)/XYZNormalization(&rABXrBC)/XYZNormalization(&rBCXrCD);
  double sinValue = XYZDotProduct(&rBC, &rCDXrAB);
  if(sinValue<0) *angle=-1*SafeArccos(cosValue);
  else *angle=SafeArccos(cosValue);
  
  return Success;
}

int FourXYZsGroupGetFourthAtom(FourXYZsGroup* pThis, double* icParam, XYZ* pAtomD){
  // calculate the fourth atom from three given atoms, please refer to the proda reference;
  double result[4] = {0.0, 0.0, 0.0, 1.0}; 
  XYZ atoms[7] = {
    {0.0, 1.0, 0.0}, 
    {0.0, 0.0, 0.0}, 
    {-1.0, 0.0, 0.0}, 
    {-1.0, 1.0, 0.0}, 
    //atoms[4], atoms[5], atoms[6] here is left uninitialized, they will be initialized below
  };
  atoms[4] = pThis->atomA;
  atoms[5] = pThis->atomB;
  atoms[6] = pThis->atomC;

  // multiply from B7, B6, ..., to B3;
  for(int i=7;i>=3;i--){
    double sa, sb, ca, cb;
    double B[4][4];
    double temp[4];
    double bond, angle, torsion;
    if(i==7){
      bond    = icParam[4]; 
      angle   = icParam[3]; 
      torsion = icParam[2]; 
    }
    else{
      XYZ vectorI, vectorJ;
      vectorI = XYZDifference(&atoms[i-1], &atoms[i-2]);
      vectorJ = XYZDifference(&atoms[i-1], &atoms[i]);
      angle = XYZAngle(&vectorI, &vectorJ);
      bond = XYZDistance(&atoms[i-1], &atoms[i]);
      torsion = GetTorsionAngle(&atoms[i-3], &atoms[i-2], &atoms[i-1], &atoms[i]);
    }
    ca = cos(angle);
    sa = sin(angle);
    cb = cos(torsion);
    sb = sin(torsion);
    B[0][0] = -ca;   B[0][1] = -sa;    B[0][2] = 0.0;  B[0][3] = -bond*ca;
    B[1][0] = sa*cb; B[1][1] = -ca*cb; B[1][2] = -sb;  B[1][3] = bond*sa*cb;
    B[2][0] = sa*sb; B[2][1] = -ca*sb; B[2][2] = cb;   B[2][3] = bond*sa*sb;
    B[3][0] = 0.0;   B[3][1] = 0.0;    B[3][2] = 0.0;  B[3][3] = 1.0;

    Matrix4By4TimesVector(temp, B, result);
    for(int j=0;j<4;j++){
      result[j] = temp[j];
    }
  }

  // the last step, equals to multiplying B2;
  pAtomD->X = -result[0] - 1.0; 
  pAtomD->Y =  result[1];
  pAtomD->Z = -result[2];
  return Success;
}


int FourXYZsGroupGetFourthAtomNew(FourXYZsGroup* pThis, double* icParam, XYZ* pAtomD){
  // another method for calculating the fourth atom from three given atoms, different from the above method
  // much faster than the above method
  XYZ ba = XYZDifference(&pThis->atomB, &pThis->atomA);
  XYZ bc = XYZDifference(&pThis->atomB, &pThis->atomC);
  XYZ baXbc = XYZCrossProduct(&ba, &bc);
  double angleABC = XYZAngle(&ba, &bc);

  if(XYZNormalization(&baXbc)<MIN_ZERO_TOLERANCE){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ZeroDivisonError;
    sprintf(errMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }

  XYZAdd(&baXbc, &pThis->atomB);
  *pAtomD = XYZRotateAround(&pThis->atomA, &pThis->atomB, &baXbc, icParam[3]-(PI-angleABC));
  *pAtomD = XYZRotateAround(pAtomD, &pThis->atomB, &pThis->atomC, icParam[2]);
  XYZMinus(pAtomD, &pThis->atomB);
  XYZScale(pAtomD, icParam[4]/XYZNormalization(pAtomD));
  XYZAdd(pAtomD, &pThis->atomC);

  return Success;
}

int FourXYZsGroupGetICParam(FourXYZsGroup* pThis, int torsionProperFlag, double* icParam){
  XYZ vectorI, vectorJ;
  // the first IC parameter; if the dihedral angle is proper, this parameter is R(AB); otherwise R(AC) instead;
  if(torsionProperFlag){
    icParam[0] = XYZDistance(&pThis->atomA, &pThis->atomB);
  }
  else{
    icParam[0] = XYZDistance(&pThis->atomA, &pThis->atomC);   
  }
  // the second IC parameter; if the dihedral angle is proper, this parameter is Theta(ABC); otherwise Theta(ACB) instead;
  if(torsionProperFlag){
    vectorI = XYZDifference(&pThis->atomA, &pThis->atomB);
    vectorJ = XYZDifference(&pThis->atomC, &pThis->atomB);
  }
  else{
    vectorI = XYZDifference(&pThis->atomA, &pThis->atomC);
    vectorJ = XYZDifference(&pThis->atomB, &pThis->atomC);  
  }
  icParam[1] = XYZAngle(&vectorI, &vectorJ);
  // the third parameter; Phi(ABCD);
  FourXYZsGroupGetTorsionAngle(pThis,&icParam[2]);
  // the fourth parameter, Theta(BCD);
  vectorI = XYZDifference(&pThis->atomB, &pThis->atomC);
  vectorJ = XYZDifference(&pThis->atomD, &pThis->atomC);
  icParam[3] = XYZAngle(&vectorI, &vectorJ);
  // the fifth parameter, R(CD);
  icParam[4] = XYZDistance(&pThis->atomC, &pThis->atomD);
  return Success;
}


int GetFourthAtom(XYZ* pAtomA, XYZ* pAtomB, XYZ* pAtomC, double* icParam, XYZ* pAtomD){
  FourXYZsGroup group;
  FourXYZsGroupCreate(&group, pAtomA, pAtomB, pAtomC, NULL);
  int result = FourXYZsGroupGetFourthAtomNew(&group, icParam, pAtomD);
  return result;
}

double GetTorsionAngle(XYZ* pAtomA, XYZ* pAtomB, XYZ* pAtomC, XYZ* pAtomD){
  FourXYZsGroup group;
  FourXYZsGroupCreate(&group, pAtomA, pAtomB, pAtomC, pAtomD);
  double torsion;
  FourXYZsGroupGetTorsionAngle(&group,&torsion);
  return torsion;
}


double RadToDeg(double rad){
    return rad*180.0/PI;
}
double DegToRad(double degree){
    return degree*PI/180.0;
}
double SafeArccos(double cosValue){
  // in header file math.h, the function acos() must have a parameter lying in [-1.0, +1.0] strictly;
  // however, in proda, the parameter may have little deviation, i.e. cos(pi) get -1.000000001, and acos( -1.000000001)
  // may not run successfully.
  // SafeArccos() is used to adjust the parameter for acos(); if the prameter is less than -1.0, alter it into -1.0;
  // if the paramter is larger than 1.0, alter it into 1.0; then call for function acos();
  if(cosValue >  1.0) cosValue = 1.0;
  else if(cosValue < -1.0) cosValue = -1.0;
  return acos(cosValue);
}

int Matrix4By4TimesVector(double result[4], double matrix[4][4], double v[4]){
  for(int i=0;i<4;i++){
    result[i] = 0.0;
    for(int j=0;j<4;j++){
      result[i] += matrix[i][j] * v[j];
    }
  }
  return Success;
}

double RandomDouble(double low, double high){
  double precision = 0.01;
  int range = (int)(fabs(high-low)/precision);
  return (rand()%range)*precision + low;
}

BOOL RadInRange(double value, double low, double high){
  while(value<low){value += (PI*2);}
  while(value>high){value -= (PI*2);}
  if(value>low) return TRUE;
  else return FALSE;
}
