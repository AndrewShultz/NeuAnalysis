
#ifndef TIMESEQUENCEFILTER_H
#define TIMESEQUENCEFILTER_H

#include <iostream>
#include <libgen.h>     
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <vector>

//ROOT Includes
#include "TTree.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TMath.h"
#include "TH2D.h"

//AraEvent includes
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraStationInfo.h"
#include "AraGeomTool.h"
#include "FFTtools.h"

using namespace std;


class TimeSequenceFilter{

double qualityParameter;
std::vector< std::vector < double > > hits;
double chNoiseEnergy[16];
double chNoiseEnergyRMS[16];
double chSignalEnergy[16];
double chSignalMaxEnergyPosition[16];
double ch_noise_std_dev[16];
double ch_noise_average[16];
std::vector<std::vector< double > > energyEnvelope;
std::vector<std::vector< double > > energyTime;
int totalHitCount;
int binnedHitCount;
double blockCount;


 double timeStep;
 int same_speed_t[800];
 int same_speed_v[800];
 int same_speed_d1[800];
 int same_speed_d2[800];
 int same_speed_h[800];


 std::vector<std::vector<std::vector< int > > > antuses;

public:

 TimeSequenceFilter();
 double getHorizontalHitDistance(int hit1, int hit2);
 double getQualityParameter(std::vector<TGraph*> gr_in, std::vector<std::vector<double> > ant_loc, int stationId, double * qualArray);
 double gen_trig_histo_single(std::vector<TGraph *> gr, TH2D *trig_pat, double *ch_std_dev, double *ch_average, double th_factor, int sum_time, int stationId);
 
 int gen_noise_std_dev(std::vector<TGraph *> gr, double *ch_noise_std_dev, double *ch_noise_average, int sum_time);
 
 double pattern_check(TH2D *trig_pat, double *pat_check_result, int cut_value, std::vector< std::vector< double > > antloc);
 
 int inter_string_check(TH2D *trig_pat, int *same_speed_t, int *same_speed_d1, int *same_speed_d2, int *same_speed_v, int *same_speed_h, int *c_s, int *c_d, std::vector< std::vector< double > > antloc);
 
 double getMaxEnergy(std::vector<TGraph *> gr);


 int GetAntUses(int i, int j, int k);
 //vector<TH1D*> GetQPHistograms();
 void GetQPHistograms(TH1D *hist0, TH1D *hist1, TH1D *hist2, TH1D *hist3, TH1D *hist4);
 /*
 double getPolarization(int stationId);
 TGraph * getEnergyGraph(int stationId, int channel);
 double getAverageNoiseEnergy(int stationId, int channel);
 double getAverageNoiseEnergyRMS(int stationId, int channel);
 double getChannelMaxEnergy(int stationId, int channel);
 double getSNR(int stationId);
 
 double getTotalHitCountValue();
 double getBinnedHitCountValue();
 */
};


#endif
