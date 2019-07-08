// ROOT headers
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
// ARA headers
#include "RawAtriStationEvent.h"
#include "RawIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "UsefulIcrrStationEvent.h"
#include "FFTtools.h"
// AS Filters headers
#include "loadBar.h"
#include "GeometryTool.h"
#include "PedestalFinder.h"
#include "HitFilter.h"
#include "TimeSequenceFilter.h"
#include "AnalyticSphereMethod.h"

// Constants
const double rad2deg = 57.2958;

// Filter Thresholds
int hitThresh = 3;
double timeSeqThresh = 1.5;
double recoAngleThresh = 38.5;


// For Ara02 only
double timeCalibThomas[16] = {-4.000,-4.421,-10.107,-7.136,-12.000,-12.421,-18.107,-15.136,-0.000,-0.421,-6.107,-3.136,-8.000,-8.421,-14.107,-11.136}; //CALIBRATION


// Main method
void filterReal(char *filename, char *outputDirectory){

  RawAtriStationEvent *rawAtriEvtPtr=0;

  HitFilter *HF = new HitFilter();
  TimeSequenceFilter *TSF = new TimeSequenceFilter();
  AnalyticSphereMethod *ASM = new AnalyticSphereMethod();

  int runNum;
  int stationID;
  int passingEvents = 0;
  string pedestalFilename;
  vector<TGraph*> gras;
  vector< vector<double> > antLocations;
  double Qarray[5];

  TFile *f = new TFile(filename);
  HF->ScanAtriFileForRMS(f);
  HF->SetRMSMultiplier(2.05);

  TTree *t = (TTree*)f->Get("eventTree");
  t->SetBranchAddress("event",&rawAtriEvtPtr);
  t->SetBranchAddress("run",&runNum);
  t->GetEntry(0);

  stationID = rawAtriEvtPtr->getStationId();
  pedestalFilename = findPedestalFile(stationID,runNum);

  TFile *outFilteredFile = new TFile(Form("%s/filteredEvents-ARA0%i-%i.root",outputDirectory,stationID,runNum),"recreate");
  TTree *outFilteredTree = new TTree("eventTree","eventTree");
  outFilteredTree->Branch("event",&rawAtriEvtPtr);

  AraEventCalibrator *araEventCalibrator = AraEventCalibrator::Instance();
  araEventCalibrator->loadAtriCalib(stationID);
  if(pedestalFilename!="") araEventCalibrator->setAtriPedFile(&pedestalFilename[0],stationID);
  
  AraGeomTool *geometryTool = AraGeomTool::Instance();
  AraStationInfo *antennaInfo = new AraStationInfo(stationID);
  AraCalAntennaInfo *pulserInfo = antennaInfo->getCalAntennaInfo(3); // 3 is D6V for ARA02, the pulser that is usually on

  //double pulserLocation[3] = {pulserInfo->getLocationXYZ()[0],pulserInfo->getLocationXYZ()[1],pulserInfo->getLocationXYZ()[2]};
  double pulserLocation[3];
  LoadPulserAra02_v3_13(3,pulserLocation[0],pulserLocation[1],pulserLocation[2]); //CALIBRATION

  int nChannels = antennaInfo->getNumAntennasByPol(AraAntPol::kVertical)+antennaInfo->getNumAntennasByPol(AraAntPol::kHorizontal);

  gras.resize(nChannels);

  /*
  for(int ich=0; ich<nChannels; ich++) {
  antLocations.push_back(vector<double>());
  antLocations[ich].push_back(antennaInfo->getAntennaInfo(ich)->getLocationXYZ()[0]);
    antLocations[ich].push_back(antennaInfo->getAntennaInfo(ich)->getLocationXYZ()[1]);
    antLocations[ich].push_back(antennaInfo->getAntennaInfo(ich)->getLocationXYZ()[2]);
    printf("Channel %i location (X,Y,Z): (%.3f,%.3f,%.3f) \n",ich,antLocations[ich][0],antLocations[ich][1],antLocations[ich][2]);
  }
  */  

  LoadAntennasAra02_v3_13(nChannels,antLocations); //CALIBRATION

  int nev = t->GetEntries();

  for(int iev=0; iev<nev; iev++) {

    loadBar(iev,nev);
    t->GetEntry(iev);

    // Skip calibration pulser events
    if(rawAtriEvtPtr->isCalpulserEvent()) continue;
    // Skip forced trigger events (noise)
    if(rawAtriEvtPtr->isSoftwareTrigger()) continue;

    UsefulAtriStationEvent *realAtriEvtPtr = new UsefulAtriStationEvent(rawAtriEvtPtr, AraCalType::kLatestCalib);

    for(int ich=0; ich<nChannels; ich++) gras[ich] = realAtriEvtPtr->getGraphFromRFChan(ich);

    int nhits = HF->ScanForHits(gras);
    //for(int ich=0; ich<nChannels; ich++) printf("Channel %i needed multiplier = %.3f \n",ich,HF->GetNeededRMSMultiplier(ich,gras[ich]));
    
    for(int ich=0; ich<nChannels; ich++) HF->hitTimes[ich] -= timeCalibThomas[ich]; //CALIBRATION

    if(nhits>hitThresh) {

      double QP = TSF->getQualityParameter(gras,antLocations,stationID,&Qarray[0]);
      if(QP>timeSeqThresh) {

	bool reconstructed = ASM->Reconstruct(antLocations,HF->hits,HF->hitTimes);
	if(reconstructed && rad2deg*ASM->theta>recoAngleThresh) {

	  //loadBar(iev,nev,Form("Position (X,Y,Z): (%.3f,%.3f,%.3f) \t Pulser: (%.3f,%.3f,%.3f) \n",ASM->aveX,ASM->aveY,ASM->aveZ,pulserLocation[0],pulserLocation[1],pulserLocation[2]));
	  passingEvents++;
	  outFilteredTree->Fill();

	}
      }
    }

    for(int ich=0; ich<nChannels; ich++) delete gras[ich];

    delete realAtriEvtPtr;

  }

  loadBarEnd();

  f->Close();

  outFilteredFile->cd();
  outFilteredTree->Write();
  outFilteredFile->Close();

  printf("Events passing: %i \n",passingEvents);

} // End of main program
