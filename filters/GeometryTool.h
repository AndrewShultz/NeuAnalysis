#ifndef GEOMETRYTOOL_H
#define GEOMETRYTOOL_H

#include <fstream>
#include <vector>

#include "AraStationInfo.h"
#include "AraGeomTool.h"

using namespace std;

void LoadAntennasAra02_v3_13(int nchan, vector< vector<double> > &antLocation) {

  int station = 2;

  AraGeomTool *fGeomTool = AraGeomTool::Instance();
  double antDelayArray[4][4]={{0}};
  double tbzCorr[16] = {0};

  ifstream ind;
  char antDelayFile[200];
  sprintf(antDelayFile, "/home/aschultz/Ferrum/geometryResultsARA%dE.txt",station);
  ind.open(antDelayFile);
  double corrections;
  if(ind.good()){
    for(int i=0;i<4;i++){
      for(int j=0; j<4;j++){
	ind >> corrections;
	antDelayArray[j][i] = corrections;
      }
    }
    for(int i=0;i<5;i++) ind >> corrections;
  }
  else printf("couldn't read position and delay correction file!! \n");
  ind.close();

  antLocation.clear();
  antLocation.resize(nchan);

  AraStationInfo *SInfo=new AraStationInfo(station);
  for(int ant=0; ant<nchan; ant++) {
    Double_t *Loc = SInfo->getAntennaInfo(ant)->getLocationXYZ();
    if(station==2 && ant==0) Loc[2] += 1.68;
    antLocation[ant].push_back(Loc[0]+antDelayArray[ant%4][0]);
    antLocation[ant].push_back(Loc[1]+antDelayArray[ant%4][1]);
    antLocation[ant].push_back(Loc[2]+antDelayArray[ant%4][2]+tbzCorr[ant]);
  }

  printf("Loaded Station %i /n",station);
  for(int ant=0; ant<nchan; ant++) printf("Antenna %i: (%.3f,%.3f,%.3f)  \n",ant,antLocation[ant][0],antLocation[ant][1],antLocation[ant][2]);

}

void LoadPulserAra02_v3_13(int pulserID, double &x, double &y, double &z) {

  int station = 2;
  AraGeomTool *fGeomTool = AraGeomTool::Instance();
  AraStationInfo *SInfo = fGeomTool->getStationInfo(station);
  AraCalAntennaInfo *CInfo=SInfo->getCalAntennaInfo(pulserID);
  x = CInfo->getLocationXYZ()[0];
  y = CInfo->getLocationXYZ()[1];
  z = CInfo->getLocationXYZ()[2];
  
  if(station==2 && (pulserID==0 || pulserID==1)) { x+=-4.6202589895470378; y+=-0.016119269379674783; z+=0.28036008549941882; }
  if(station==2 && (pulserID==2 || pulserID==3)) z+=-5.3502286411149811;

}

#endif
