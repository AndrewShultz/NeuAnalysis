#ifndef ANALYTIC_SPHERE_METHOD_H
#define ANALYTIC_SPHERE_METHOD_H

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "assert.h"

#include "TMath.h"

using namespace std;

class AnalyticSphereMethod {

 public:

  double aveX;
  double aveY;
  double aveZ;

  double stationCenterX;
  double stationCenterY;
  double stationCenterZ;

  double theta;
  double phi;
  double r;

  AnalyticSphereMethod();

  inline double GetX(){return _x;};
  inline double GetY(){return _y;};
  inline double GetZ(){return _z;};
  inline double GetT(){return _t;};
  inline double Dist3d(double x1, double y1, double z1,double x2, double y2, double z2){return sqrt(SQ(x2-x1) + SQ(y2-y1) + SQ(z2-z1));}
  inline void SetIndexOfRefractionParameters(double a, double b, double c) { _nofzA = a; _nofzB = b; _nofzC = c; }
  inline double SQ(double a){return(a*a);}

  void ConvertCartesian2Spherical(double origx, double origy, double origz, double vecx, double vecy, double vecz, double &Theta, double &Phi, double &R);
  void Orderer(double *values, int *order, int n);
  double Nofz(double z);
  double nTimeExpected(double z1, double z2, double dtot);
  bool Reconstruct(vector< vector<double> > antLocation, vector<bool> hits, vector<double> hitTimes);
  void SetAnt(double x, double y, double z, double t);
  bool Solve();
  void QUAD(double A, double B, double C, double &Pos, double &Neg);

 private:

  double   _AntInfo[4][4];
  double   _CIceP1;
  double   _CIceP2;
  double   _CIceP3;
  double   _x;
  double   _y;
  double   _z;
  double   _t;
  int      _SetCnt;

  double _nofzA;
  double _nofzB;
  double _nofzC;

};

#endif
