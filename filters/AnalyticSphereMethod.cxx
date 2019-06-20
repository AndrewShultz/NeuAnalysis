#include "AnalyticSphereMethod.h"

AnalyticSphereMethod::AnalyticSphereMethod() {

  _x           = -9996;
  _y           = -9997;
  _z           = -9998;
  _t           = -9999;
  _SetCnt      = 0;

  _nofzA = 1.78;
  _nofzB = -.24 * _nofzA;
  _nofzC = 0.016;

}

//
//
//
//
//
//
//
//

void AnalyticSphereMethod::ConvertCartesian2Spherical(double origx, double origy, double origz, double vecx, double vecy, double vecz, double &Theta, double &Phi, double &R) {

  R = Dist3d(origx,origy,origz,vecx,vecy,vecz);

  if(R==0) {
    Phi   = nan("");
    Theta = nan("");
    return;
  }

  else {
    Theta = acos((vecz-origz)/R);

    double top = (vecy - origy);
    double bottom = (vecx - origx);

    Phi = atan2(top,bottom);

    if(Phi<0)
      Phi += 2*TMath::Pi();

  }

}

//
//
//
//
//
//
//
//

void AnalyticSphereMethod::Orderer(double *values, int *order, int n) {

  double valord[n];
  double high=0;

  for (int i=0; i<n; i++) {
    valord[i] = values[i];
    high += fabs(values[i]);
  }

  for(int i = 0;i<n;i++) {
    double checker = high;
    for (int j = 0;j<n;j++) {
      if (valord[j]<=checker) {
	checker = values[j];
	order[i] = j;
      }
    }
    valord[order[i]] = high+1;
  }

}

//
//
//
//
//
//
//
//

double AnalyticSphereMethod::Nofz(double z) {

  if(z>=0) return 1.0;

  double n = _nofzA + _nofzB * exp(_nofzC*z);

  return n;
}

//
//
//
//
//
//
//
//

double AnalyticSphereMethod::nTimeExpected(double z1, double z2, double dtot) {

  const double cVac=0.2998;
  double ctheta;
  double ftimexp=dtot/cVac;

  double zmax=z2;
  double zmin=z1;
  if(zmin>zmax){
    zmin=z2;
    zmax=z1;
  }

  if(dtot>0) ctheta=fabs(z2-z1)/dtot;
  else return 0;

  if(ctheta!=0){

    if(zmax<0. && zmin<0.)
      ftimexp=(_nofzA*(zmax-zmin)+(_nofzB/_nofzC)*(exp(_nofzC*zmax)-exp(_nofzC*zmin)))/(cVac * ctheta);

    if(zmax>=0. && zmin<0.)
      ftimexp=(_nofzA*(0.-zmin)+(_nofzB/_nofzC)*(exp(_nofzC*0.)-exp(_nofzC*zmin))+zmax)/(cVac * ctheta);

    if(zmax>=0. && zmin>=0.)
      ftimexp=(zmax-zmin)/(cVac * ctheta);

  }

  return ftimexp;

}

//
//
//
//
//
//
//
//

bool AnalyticSphereMethod::Reconstruct(vector< vector<double> > antLocation, vector<bool> hits, vector<double> hitTimes) {

  stationCenterX = 0.0;
  stationCenterY = 0.0;
  stationCenterZ = 0.0;

  int nAnt = antLocation.size();

  for(int iant=0; iant<nAnt; iant++) {
    stationCenterX += antLocation[iant][0];
    stationCenterY += antLocation[iant][1];
    stationCenterZ += antLocation[iant][2];
  }

  stationCenterX /= nAnt;
  stationCenterY /= nAnt;
  stationCenterZ /= nAnt;

  aveX = 0.0;
  aveY = 0.0;
  aveZ = 0.0;
  int solcnt = 0;

  for(int iant1=0; iant1<nAnt; iant1++) {

    if(!hits[iant1]) continue;
    _AntInfo[0][0] = antLocation[iant1][0];
    _AntInfo[0][1] = antLocation[iant1][1];
    _AntInfo[0][2] = antLocation[iant1][2];
    _AntInfo[0][3] = hitTimes[iant1];


    for(int iant2=iant1+1; iant2<nAnt; iant2++) {

      if(!hits[iant2]) continue;
      _AntInfo[1][0] = antLocation[iant2][0];
      _AntInfo[1][1] = antLocation[iant2][1];
      _AntInfo[1][2] = antLocation[iant2][2];
      _AntInfo[1][3] = hitTimes[iant2];

      for(int iant3=iant2+1; iant3<nAnt; iant3++) {

	if(!hits[iant3]) continue;
	_AntInfo[2][0] = antLocation[iant3][0];
	_AntInfo[2][1] = antLocation[iant3][1];
	_AntInfo[2][2] = antLocation[iant3][2];
	_AntInfo[2][3] = hitTimes[iant3];

	for(int iant4=iant3+1; iant4<nAnt; iant4++) {

	  if(!hits[iant4]) continue;
	  _AntInfo[3][0] = antLocation[iant4][0];
	  _AntInfo[3][1] = antLocation[iant4][1];
	  _AntInfo[3][2] = antLocation[iant4][2];
	  _AntInfo[3][3] = hitTimes[iant4];

	  _SetCnt = 4;

	  if(Solve()) {
	    aveX += _x;
	    aveY += _y;
	    aveZ += _z;
	    solcnt++;
	  }

	}
      }
    }
  }

  if(solcnt>0) {
    aveX /= solcnt;
    aveY /= solcnt;
    aveZ /= solcnt;
    ConvertCartesian2Spherical(stationCenterX,stationCenterY,stationCenterZ,aveX,aveY,aveZ,theta,phi,r);
    return true;
  }

  return false;

}

//
//
//
//
//
//
//
//

void AnalyticSphereMethod::SetAnt( double x, double y, double z, double t) {

  if(_SetCnt>=4){
    printf("AnalyticSphereMethod:  Cannot set anymore antennas right now... ASM Ready to solve.\n");
    return;
  }

  else {
    _AntInfo[_SetCnt][0] = x;
    _AntInfo[_SetCnt][1] = y;
    _AntInfo[_SetCnt][2] = z;
    _AntInfo[_SetCnt][3] = t;
    _SetCnt++;
  }

}

//
//
//
//
//
//
//
//

bool AnalyticSphereMethod::Solve() {

  if(_SetCnt!=4){
    printf("AnalyticSphereMethod:  Need 4 antennas set... Criterion not met, exiting... \n");
    assert(0);
  }

  _SetCnt = 0;

  // Find Order of Antenna Hits

  double tempT[4];
  tempT[0] = _AntInfo[0][3];
  tempT[1] = _AntInfo[1][3];
  tempT[2] = _AntInfo[2][3];
  tempT[3] = _AntInfo[3][3];

  Int_t order[4];

  Orderer(tempT,order,4);

  // Parameters

  double x1, y1, z1;
  double x2, y2, z2;
  double x3, y3, z3;
  double x4, y4, z4;

  x1 = _AntInfo[order[0]][0];
  y1 = _AntInfo[order[0]][1];
  z1 = _AntInfo[order[0]][2];

  x2 = _AntInfo[order[1]][0];
  y2 = _AntInfo[order[1]][1];
  z2 = _AntInfo[order[1]][2];

  x3 = _AntInfo[order[2]][0];
  y3 = _AntInfo[order[2]][1];
  z3 = _AntInfo[order[2]][2];

  x4 = _AntInfo[order[3]][0];
  y4 = _AntInfo[order[3]][1];
  z4 = _AntInfo[order[3]][2];

  // Find Index of Refraction

  _CIceP1 = (TMath::C()/(1e9))/Nofz((z1+z2)/2);
  _CIceP2 = (TMath::C()/(1e9))/Nofz((z1+z3)/2);
  _CIceP3 = (TMath::C()/(1e9))/Nofz((z1+z4)/2);
  //printf(" 1-2: (%.3f,%.3f): %.3f \t 1-3: (%.3f,%.3f): %.3f \t 1-4: (%.3f,%.3f): %.3f \n",z1,z2,NofZ((z1+z2)/2),z1,z3,NofZ((z1+z3)/2),z1,z4,NofZ((z1+z4)/2));
  /*
    _CIceP1 = (TMath::C()/(1e9))/1.45;
    _CIceP2 = (TMath::C()/(1e9))/1.45;
    _CIceP3 = (TMath::C()/(1e9))/1.45;
  */
  // Derived Parameters

  double deltaD12, deltaD13, deltaD14;

  deltaD12 = (_AntInfo[order[1]][3]-_AntInfo[order[0]][3])*_CIceP1;
  deltaD13 = (_AntInfo[order[2]][3]-_AntInfo[order[0]][3])*_CIceP2;
  deltaD14 = (_AntInfo[order[3]][3]-_AntInfo[order[0]][3])*_CIceP3;

  double RD12, RD13, RD14;

  RD12 = SQ(x2)+SQ(y2)+SQ(z2) - (SQ(x1)+SQ(y1)+SQ(z1)) - SQ(deltaD12);
  RD13 = SQ(x3)+SQ(y3)+SQ(z3) - (SQ(x1)+SQ(y1)+SQ(z1)) - SQ(deltaD13);
  RD14 = SQ(x4)+SQ(y4)+SQ(z4) - (SQ(x1)+SQ(y1)+SQ(z1)) - SQ(deltaD14);


  // Z(d) Terms

  double alpha0, alpha1, alpha2, alpha3;

  alpha0 = ( y1-y4 - (y1-y2)*((x1-x4)/(x1-x2)) ) / (y1-y3-(y1-y2)*((x1-x3)/(x1-x2)));

  alpha1 =  alpha0 * (RD12*((x1-x3)/(x1-x2)) - RD13) + RD14 - RD12*((x1-x4)/(x1-x2));

  alpha2 = 2 * ( alpha0 * (deltaD13 - deltaD12*((x1-x3)/(x1-x2))) + deltaD12*((x1-x4)/(x1-x2)) - deltaD14);

  alpha3 = 2 * ( alpha0 * ( (z1-z2)*((x1-x3)/(x1-x2)) -z1+z3) + z1 - z4 - (z1-z2)*((x1-x4)/(x1-x2)));


  // X(d) Terms

  double beta0, beta1, beta2, beta3;

  beta0 = ((z1-z4)*((y1-y3)/(y1-y4)) -z1+z3) / (z1-z2-(z1-z4)*((y1-y2)/(y1-y4)));

  beta1 = beta0 * (RD12-RD14*((y1-y2)/(y1-y4))) + RD13 - RD14*((y1-y3)/(y1-y4));

  beta2 = 2 * ( beta0 * (deltaD14*((y1-y2)/(y1-y4)) - deltaD12) + deltaD14*((y1-y3)/(y1-y4)) - deltaD13);

  beta3 = 2 * ( beta0 * (x1-x2-(x1-x4)*((y1-y2)/(y1-y4))) + x1 - x3 - (x1-x4)*((y1-y3)/(y1-y4)) );

  // Y(d) Terms

  double gamma0, gamma1, gamma2, gamma3;

  gamma0 = ((x1-x3)*((z1-z2)/(z1-z3)) -x1+x2) / (x1-x4 - (x1-x3)*((z1-z4)/(z1-z3)));

  gamma1 = gamma0 * (RD14-RD13*((z1-z4)/(z1-z3))) + RD12 - RD13*((z1-z2)/(z1-z3));

  gamma2 = 2 * ( gamma0 * (deltaD13*((z1-z4)/(z1-z3)) - deltaD14) + deltaD13*((z1-z2)/(z1-z3)) - deltaD12);

  gamma3 = 2 * ( gamma0 * (y1-y4-(y1-y3)*((z1-z4)/(z1-z3))) + y1 - y2 - (y1-y3)*((z1-z2)/(z1-z3)));

  // Distance Terms

  double QuadTerms, LinTerms, Cons;

  QuadTerms = SQ(beta2/beta3) + SQ(gamma2/gamma3) + SQ(alpha2/alpha3) - 1;
  LinTerms = 2*( (beta1*beta2)/SQ(beta3) + (gamma1*gamma2)/SQ(gamma3) + (alpha1*alpha2)/SQ(alpha3) + x1*(beta2/beta3) + y1*(gamma2/gamma3) + z1*(alpha2/alpha3) );
  Cons = SQ(beta1/beta3) + SQ(gamma1/gamma3) + SQ(alpha1/alpha3) + 2*(x1*(beta1/beta3) + y1*(gamma1/gamma3) + z1*(alpha1/alpha3)) +SQ(x1)+SQ(y1)+SQ(z1);

  // Solve for d with quadratic equation

  double timepos, timeneg;

  QUAD(QuadTerms,LinTerms,Cons,timepos,timeneg);

  // Positive solution is taken
  //printf("TimePos: %.6f \t TimeNeg: %.6f \n",timepos,timeneg);
  if(timepos>0 && timeneg<0) {
    _x = -(beta2*timepos + beta1)/beta3;
    _y = -(gamma2*timepos + gamma1)/gamma3;
    _z = -(alpha2*timepos + alpha1)/alpha3;
  }

  else if(timepos<0 && timeneg>0) {
    _x = -(beta2*timeneg + beta1)/beta3;
    _y = -(gamma2*timeneg + gamma1)/gamma3;
    _z = -(alpha2*timeneg + alpha1)/alpha3;
  }

  // If both negative or positive returns false

  else return false;

  // Calculate time to antenna 1 based on coordinates found

  _t = nTimeExpected(z1,_z,Dist3d(x1,y1,z1,_x,_y,_z));

  return true;

}

//
//
//
//
//
//
//
//

void AnalyticSphereMethod::QUAD(double A, double B, double C, double &Pos, double &Neg) {
  double determinant = B*B - 4*A*C;
  Pos = (-B + sqrt(determinant))/(2*A);
  Neg = (-B - sqrt(determinant))/(2*A);
}
