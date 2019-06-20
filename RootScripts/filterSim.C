
//C++ includes
#include "stdlib.h"
#include "stdio.h"
#include <vector>

//ROOT includes
#include "TFile.h"
#include "TTree.h"

// ARAROOT includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"

// ARASIM includes
#include "EarthModel.h"
#include "Report.h"
#include "IceModel.h"
#include "Detector.h"
#include "Event.h"
#include "Settings.h"

// AS Filters headers
#include "loadBar.h"
#include "PedestalFinder.h"
#include "HitFilter.h"
#include "TimeSequenceFilter.h"
#include "AnalyticSphereMethod.h"

using namespace std;

// Constants
const double rad2deg = 57.2958;

// Filter Thresholds
int hitThresh = 3;
double timeSeqThresh = 1.5;
double recoAngleThresh = 38.5;

// Global variables
Detector *detector=0;
Settings *settings=0;
IceModel *icemodel=0;
Report *report=0;
Event *event=0;

UsefulAtriStationEvent *realAtriEvtPtr=0;

int runNum;
int stationID;
vector< vector<double> > antLocations;
vector<int> polarization;
vector<int> arbChNum;
int nChannels = 0;

bool geometryLoaded = false;
bool useAraTree2 = false;
bool useEventTree = false;

double stationX;
double stationY;
double stationZ;

TTree *inTree = 0;
TTree *inTree2 = 0;
TTree *outFilteredTree=0;

// Functions to set TTree and get waveforms
void DetermineAndSetTTree(TFile *inputFile);
void GetEventWaveforms(vector<TGraph*> &gras);

// Main method
void filterSim(char *filename, char *outputDirectory){

  HitFilter *HF = new HitFilter();
  TimeSequenceFilter *TSF = new TimeSequenceFilter();
  AnalyticSphereMethod *ASM = new AnalyticSphereMethod();

  int passingEvents = 0;
  vector<TGraph*> gras;
  vector<double> rms;
  double Qarray[5];

  string tempFilename = filename;
  std::string::size_type pos1 = tempFilename.find("run");
  std::string::size_type pos2 = tempFilename.find(".root");
  if(pos1==std::string::npos || pos2==std::string::npos) runNum = 0;
  else runNum = atoi(tempFilename.substr(pos1+3,pos2-pos1-3).c_str());
  printf("Run Number = %i \n",runNum);

  TFile *inputFile = new TFile(filename);

  // This will set the TTree to either AraTree2 (preferred) or eventTree, whichever is available.
  // Will also set station geometry (antLocations).
  DetermineAndSetTTree(inputFile);

  if( !(useAraTree2 || useEventTree) || !geometryLoaded ) {
    printf("File unusable...\n");
    if(!(useAraTree2 || useEventTree)) printf("No usable event TTree\n");
    if(!geometryLoaded) printf("Cannot load station geometry\n");
    printf("Exiting...\n");
    exit(0);
  }

  for(int ich=0; ich<nChannels; ich++) gras.push_back(0);
  for(int ich=0; ich<nChannels; ich++) rms.push_back(43.0);
  HF->SetRMS(rms);

  int nev = inTree->GetEntries();

  for(int iev=0; iev<nev; iev++) {

    inTree->GetEntry(iev);
    loadBar(iev,nev);

    if(useAraTree2)
      if(report->stations[0].Global_Pass<=0) continue; // Did the event trigger? Need to check when using AraTree2.

    GetEventWaveforms(gras);

    int nhits = HF->ScanForHits(gras);
    //for(int ich=0; ich<nChannels; ich++) printf("Channel %i needed multiplier = %.3f \n",ich,HF->GetNeededRMSMultiplier(ich,gras[ich]));
    if(nhits>hitThresh) {

      double QP = TSF->getQualityParameter(gras,antLocations,stationID,&Qarray[0]);
      if(QP>timeSeqThresh) {

	bool reconstructed = ASM->Reconstruct(antLocations,HF->hits,HF->hitTimes);
	if(reconstructed && rad2deg*ASM->theta>recoAngleThresh) {

	  if(useAraTree2) loadBar(iev,nev,Form("Neutrino Position (X,Y,Z): (%.3f,%.3f,%.3f) \n",event->Nu_Interaction[0].posnu.GetX()-stationX,event->Nu_Interaction[0].posnu.GetY()-stationY,event->Nu_Interaction[0].posnu.GetZ()-stationZ));
	  //printf("Position (X,Y,Z): (%.3f,%.3f,%.3f) \t Pulser: (%.3f,%.3f,%.3f) \n",ASM->aveX,ASM->aveY,ASM->aveZ);

	  passingEvents++;
	  //outFilteredTree->Fill();

	}
      }
    }

  }

  inputFile->Close();

  printf("Events passing: %i \n",passingEvents);

}

//
//
//
//
//
//
//
//

void DetermineAndSetTTree(TFile *inputFile) {

  if(inputFile->FindKey("AraTree")) {

    inTree2 = (TTree*)inputFile->Get("AraTree");
    inTree2->SetBranchAddress("detector",&detector);
    inTree2->SetBranchAddress("icemodel",&icemodel);
    inTree2->SetBranchAddress("settings",&settings);

    inTree2->GetEntry(0);

    int nStations = detector->params.number_of_stations;
    int nStringsPerStation = detector->params.number_of_strings_station;
    int nAntennasPerString = detector->params.number_of_antennas_string;
    nChannels = nStringsPerStation*nAntennasPerString;

    stationX = detector->stations[0].GetX();
    stationY = detector->stations[0].GetY();
    stationZ = detector->stations[0].GetZ();

    stationID = detector->stations[0].StationID;

    int ichann = 0;

    antLocations.resize(nChannels);
    polarization.resize(nChannels);
    arbChNum.resize(nChannels);

    for(int istring=0; istring<nStringsPerStation; istring++) {
      for(int iant=0; iant<nAntennasPerString; iant++) {
	int channel = detector->GetChannelfromStringAntenna(stationID,istring,iant,settings);
	antLocations[channel].push_back(detector->stations[0].strings[istring].antennas[iant].GetX()-stationX);
	antLocations[channel].push_back(detector->stations[0].strings[istring].antennas[iant].GetY()-stationY);
	antLocations[channel].push_back(detector->stations[0].strings[istring].antennas[iant].GetZ()-stationZ);
	polarization[channel] = detector->stations[0].strings[istring].antennas[iant].type;
	arbChNum[channel] = ichann;
	ichann++;
      }
    }

    for(int ich=0; ich<nChannels; ich++) printf("Channel %i location (X,Y,Z): (%.3f,%.3f,%.3f) \t Polarization: %i \n",ich,antLocations[ich][0],antLocations[ich][1],antLocations[ich][2],polarization[ich]);

    geometryLoaded = true;

  }

  if(inputFile->FindKey("AraTree2")) {
    inTree = (TTree*)inputFile->Get("AraTree2");
    inTree->SetBranchAddress("report",&report);
    inTree->SetBranchAddress("event",&event);
    if(inTree->GetEntries()>0) useAraTree2 = true;
  }

  if(!useAraTree2 && inputFile->FindKey("eventTree")) {

    inTree = (TTree*)inputFile->Get("eventTree");
    inTree->SetBranchAddress("UsefulAtriStationEvent",&realAtriEvtPtr);
    if(inTree->GetEntries()>0) useEventTree = true;
    else { printf("No usable TTree's or no events, exiting...\n"); exit(0); }

    if(!geometryLoaded) {

      inTree->GetEntry(0);

      stationID = realAtriEvtPtr->getStationId();

      AraGeomTool *geometryTool = AraGeomTool::Instance();
      AraStationInfo *antennaInfo=new AraStationInfo(stationID);

      nChannels = antennaInfo->getNumAntennasByPol(AraAntPol::kVertical)+antennaInfo->getNumAntennasByPol(AraAntPol::kHorizontal);

      for(int ich=0; ich<nChannels; ich++) {
	antLocations[ich].push_back(antennaInfo->getAntennaInfo(ich)->getLocationXYZ()[0]);
	antLocations[ich].push_back(antennaInfo->getAntennaInfo(ich)->getLocationXYZ()[1]);
	antLocations[ich].push_back(antennaInfo->getAntennaInfo(ich)->getLocationXYZ()[2]);
	printf("Channel %i location (X,Y,Z): (%.3f,%.3f,%.3f) \n",ich,antLocations[ich][0],antLocations[ich][1],antLocations[ich][2]);
      }

      geometryLoaded = true;

    }

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

void GetEventWaveforms(vector<TGraph*> &gras) {

  if(useAraTree2) {
    for(int ich=0; ich<nChannels; ich++) {
      if(gras[ich]) delete gras[ich];
      gras[ich] = report->getWaveform(detector,arbChNum[ich]);
    }
  }

  else if(useEventTree) {
    for(int ich=0; ich<nChannels; ich++) {
      if(gras[ich]) delete gras[ich];
      gras[ich] = realAtriEvtPtr->getGraphFromRFChan(ich);
    }
  }

}
