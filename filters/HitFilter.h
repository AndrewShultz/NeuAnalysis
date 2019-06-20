
#ifndef HITFILTER_H
#define HITFILTER_H

#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "AraGeomTool.h"
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"

#include "PedestalFinder.h"

using namespace std;

class HitFilter {

 public:

  HitFilter();

  //variables
  int runNum;
  int stationID;
  int nChannels;
  vector<bool> hits;
  vector<double> hitTimes;
  vector<double> noiseRMS;
  double RMSMultiplier;
  bool scanSuccess;
  int nsums;

  //main utility functions
  void ScanAtriFileForRMS(TFile *f);
  int ScanForHits(vector<TGraph*> &gras);
  double GetNeededRMSMultiplier(int channel, TGraph *gra);

  //getters
  double GetRMS(int channel);
  double GetRMSMultiplier();
  double GetThreshold(int channel);

  //setters
  void SetRMS(vector<double> &RMS);
  double SetRMSMultiplier(double mult);

 private:

};



#endif


