
#include "HitFilter.h"

HitFilter::HitFilter() {

  RMSMultiplier = 4.5;
  scanSuccess = false;
  nsums = 10;

}

//
//
//
//
//
//
//
//

void HitFilter::ScanAtriFileForRMS(TFile *f) {

  RawAtriStationEvent *rawAtriEvtPtr=0;
  TTree *t = (TTree*) f->Get("eventTree");
  t->SetBranchAddress("event",&rawAtriEvtPtr);
  t->SetBranchAddress("run",&runNum);
  int nevt = t->GetEntries();

  int cnt = 0;

  t->GetEntry(0);
  stationID = rawAtriEvtPtr->getStationId();
  string pedestalFilename = findPedestalFile(stationID,runNum);

  AraEventCalibrator *araEventCalibrator = AraEventCalibrator::Instance();
  araEventCalibrator->loadAtriCalib(stationID);
  if(pedestalFilename!="") araEventCalibrator->setAtriPedFile(&pedestalFilename[0],stationID);

  AraGeomTool *geometryTool = AraGeomTool::Instance();
  AraStationInfo *antennaInfo=new AraStationInfo(stationID);
  nChannels = antennaInfo->getNumAntennasByPol(AraAntPol::kVertical)+antennaInfo->getNumAntennasByPol(AraAntPol::kHorizontal);

  noiseRMS.resize(nChannels);

  for(int j=0; j<nevt; j++) {

    t->GetEntry(j);

    if(rawAtriEvtPtr->isSoftwareTrigger()) {

      UsefulAtriStationEvent *realAtriEvtPtr = new UsefulAtriStationEvent(rawAtriEvtPtr, AraCalType::kLatestCalib);

      for(int ichl=0; ichl<nChannels; ichl++) {
	TGraph *gra = realAtriEvtPtr->getGraphFromRFChan(ichl);
	noiseRMS[ichl] += gra->GetRMS(2);
	delete gra;
      }

      cnt++;
      delete realAtriEvtPtr;

    }

  }

  if(cnt>0) {
    for(int ichl=0; ichl<nChannels; ichl++) {
      noiseRMS[ichl] /= cnt;
      printf("Channel %i RMS: %.3f \n",ichl,noiseRMS[ichl]);
    }
    scanSuccess = true;
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

int HitFilter::ScanForHits(vector<TGraph*> &gras) {

  if(!scanSuccess) { printf("Need to set RMS, either by hand or by ScanAtriFileForRMS() \n"); return 0; }
  if(noiseRMS.size()!=gras.size()) { printf("Input TGraph vector size does not match Noise RMS vector size!!! \n"); return 0; }

  int nhits = 0;
  hits.resize(noiseRMS.size());
  hitTimes.resize(noiseRMS.size());

  for(int igra=0; igra<gras.size(); igra++) {

    bool hit = false;
    hits[igra] = false;

    double val = 0.0;

    for(int ipnt=0; ipnt<nsums && ipnt<gras[igra]->GetN(); ipnt++) val += fabs(gras[igra]->GetY()[ipnt]);

    if(val/nsums>RMSMultiplier*noiseRMS[igra]) {
      hit = true;
      nhits++;
      hits[igra] = true;
      hitTimes[igra] = gras[igra]->GetX()[0];
    }

    for(int ipnt=0; ipnt+nsums<gras[igra]->GetN() && !hit; ipnt++) {

      val -= fabs(gras[igra]->GetY()[ipnt]);
      val += fabs(gras[igra]->GetY()[ipnt+nsums]);

      if(val/nsums>RMSMultiplier*noiseRMS[igra]) {
	hit = true;
	nhits++;
	hits[igra] = true;
	hitTimes[igra] = gras[igra]->GetX()[ipnt];
      }

    }

  }

  return nhits;

}

//
//
//
//
//
//
//
//

double HitFilter::GetNeededRMSMultiplier(int channel, TGraph *gra) {

  if(!scanSuccess) { printf("Need to set RMS, either by hand or by ScanAtriFileForRMS() \n"); return 0; }
  if(channel>noiseRMS.size()) { printf("Requested channel has larger index than Noise RMS vector size!!! \n"); return 0; }

  double max = 0.0;
  double val = 0.0;

  for(int ipnt=0; ipnt<nsums && ipnt<gra->GetN(); ipnt++) val += fabs(gra->GetY()[ipnt]);
  if(max<val/nsums) max = val/nsums;

  for(int ipnt=0; ipnt+nsums<gra->GetN(); ipnt++) {
    val -= fabs(gra->GetY()[ipnt]);
    val += fabs(gra->GetY()[ipnt+nsums]);
    if(max<val/nsums) max = val/nsums;
  }

  return (max/noiseRMS[channel]);

}

//
//
//
//
//
//
//
//

double HitFilter::GetRMS(int channel) {
  if(!scanSuccess) { printf("Need to set RMS, either by hand or by ScanAtriFileForRMS() \n"); return 0; }
  else return noiseRMS[channel];
}

//
//
//
//
//
//
//
//

double HitFilter::GetRMSMultiplier() {
  if(!scanSuccess) { printf("Need to set RMS, either by hand or by ScanAtriFileForRMS() \n"); return 0; }
  else return RMSMultiplier;
}

//
//
//
//
//
//
//
//

double HitFilter::GetThreshold(int channel) {
  if(!scanSuccess) { printf("Need to set RMS, either by hand or by ScanAtriFileForRMS() \n"); return 0; }
  else return RMSMultiplier*noiseRMS[channel];
}

//
//
//
//
//
//
//
//

void HitFilter::SetRMS(vector<double> &RMS) {
  noiseRMS.resize(RMS.size());
  scanSuccess = true;
  for(int ichl; ichl<RMS.size(); ichl++) noiseRMS[ichl] = RMS[ichl];
}

//
//
//
//
//
//
//
//

double HitFilter::SetRMSMultiplier(double mult) {
  RMSMultiplier = mult;
}
