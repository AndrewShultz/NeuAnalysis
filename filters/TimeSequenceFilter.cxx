
#include "TimeSequenceFilter.h"

TimeSequenceFilter::TimeSequenceFilter() 
{
   //Default Constructor
   hits.clear();
   antuses.resize(16);
   for(int i=0; i<16; i++) {
     antuses[i].resize(5);
     for(int j=0; j<5; j++) {
       for(int k=0; k<16; k++) {
	 antuses[i][j].push_back(0);
       }
     }
   }

   for(int i=0; i<800; i++) {
     same_speed_t[i] = 0;
     same_speed_v[i] = 0;
     same_speed_d1[i] = 0;
     same_speed_d2[i] = 0;
     same_speed_h[i] = 0;
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
//
//
int TimeSequenceFilter::GetAntUses(int i, int j, int k) {
  return antuses[i][j][k];
}
//
//
//
//
//
//
//
//
//
//
void TimeSequenceFilter::GetQPHistograms(TH1D *hist0, TH1D *hist1, TH1D *hist2, TH1D *hist3, TH1D *hist4) {

  int colors[7] = {kPink-1,kBlue-8,kGreen+2,kOrange+7,kCyan+3,kBlue,kMagenta};
  /*vector<TH1D*> hist;
  for(int i=0; i<5; i++) {
    hist.push_back(new TH1D(Form("hist%i",i),"",800,-0.5,799.5));
    }*/
  printf("Time Step: %.6f \n",timeStep);
  for(int i=0; i<800; i++) hist0->Fill(timeStep/15.0*double(i-400),same_speed_t[i]);
  for(int i=0; i<800; i++) hist1->Fill(timeStep/15.0*double(i-400),same_speed_v[i]);
  for(int i=0; i<800; i++) hist2->Fill(timeStep/15.0*double(i-400),same_speed_d1[i]);
  for(int i=0; i<800; i++) hist3->Fill(timeStep/15.0*double(i-400),same_speed_d2[i]);
  for(int i=0; i<800; i++) hist4->Fill(timeStep/15.0*double(i-400),same_speed_h[i]);
			    
  //return hist;

}
//
//
//
//
//
//
//
//
//
//
double TimeSequenceFilter::getHorizontalHitDistance(int hit1, int hit2) {
  return sqrt( (hits[hit1][1]-hits[hit2][1])*(hits[hit1][1]-hits[hit2][1]) + (hits[hit1][2]-hits[hit2][2])*(hits[hit1][2]-hits[hit2][2]) );
}
//
//
//
//
//
//
//
//
//
//
int TimeSequenceFilter::inter_string_check(TH2D *trig_pat, int *same_speed_t, int *same_speed_d1, int *same_speed_d2, int *same_speed_v, int *same_speed_h, int *c_s, int *c_d, std::vector< std::vector< double > > antloc) {

  hits.clear();
  std::vector< double > singleHit;
  double string_pos[4];
  //Unfortunately the string numbering is different for ARA02 and ARA03. 
  //The following lines assine a new number to each string, which is consistent throughout the stations.
  //We need to be careful with this, since it is based on the geometry. Still it should be safe.
  for(int st=0;st<4;st++) {
    if(antloc[st][0]>0.0 && antloc[st][1]> 0.0 )string_pos[st] = 1;
    else if(antloc[st][0]>0.0 && antloc[st][1]<0.0 )string_pos[st] = 2;
    else if(antloc[st][0]<0.0 && antloc[st][1]>0.0 )string_pos[st] = 3;
    else if(antloc[st][0]<0.0 && antloc[st][1]<0.0 )string_pos[st] = 4;
  }
  
  //A few variables which will be used in the following
  int max_bin = trig_pat->GetNbinsX();				//The number of bins
  timeStep = trig_pat->GetXaxis()->GetBinWidth(1);		//the binsize
  double timeRange = trig_pat->GetXaxis()->GetXmax();		//the time range
  int count = 0;
  
  //	cout << "Bin number: " << max_bin << endl;  //DEBUG
  
  //Now the hits vector is filled:
  //All hits are written into a one dimensional vector with information about:
  //	1)hit time
  //	2)antenna x-coord
  //	3)antenna y-coord
  //	4)antenna z-coord
  //	5)string number (the new assigned one)
  //	6)antenna number
  for(int i=1;i<max_bin+1;i++) {
    for(int j=1;j<17;j++) {
      if(trig_pat->GetBinContent(i,j)>0) {
	singleHit.push_back( i*timeStep );  //
	singleHit.push_back( antloc[j-1][0] );
	singleHit.push_back( antloc[j-1][1] );
	singleHit.push_back( antloc[j-1][2] );
	singleHit.push_back( string_pos[(j-1)%4] );
	singleHit.push_back( (j-1)/4 );
	singleHit.push_back(j-1);
	hits.push_back(singleHit);
	singleHit.clear();
      }
    }
  }
  
  //This gives us the number of hits.
  count = hits.size();
  
  double distance = 0;
  
  //Now loop through all possible hit combinations and check their relation to each other.
  //If they fullfill the appolied conditions, they are selected to be histogrammed. The histograms
  //are separated for:
  //	hit pair antennas on the same string (vertical pairs)
  //	hit pair antennas at the same depth level (Note: it is not exactly the same depth, just the antenna number on the string)
  //		These are separated again in the four horizontal directions: two lateral, two diagonal
  //	In total we have 5 histograms
  for(int gg = 0; gg<count;gg++) {
    for(int kk=gg+1;kk<count;kk++) {
      //First check if pair has the same antenna number, means they are at the same depth level.
      if(hits[kk][5]==hits[gg][5]) {
	
	if(TimeSequenceFilter::getHorizontalHitDistance(kk,gg)!=0) {//This is to make sure that they are not the same antenna: Could be replaced by a string check
	  
	  //Now the string numbers are used to determine the horizontal direction:
	  //	lateral 1) string 1-2 and string 3-4
	  if(TMath::Abs(hits[kk][4]-hits[gg][4])==1 && hits[kk][4]+hits[gg][4]!=5) {
	    distance =(-hits[kk][4]+hits[gg][4])*TimeSequenceFilter::getHorizontalHitDistance(kk,gg);
	    same_speed_t[400 + TMath::Nint(15.0/timeStep*(hits[kk][0] - hits[gg][0])/distance) ]++;
	    antuses[hits[gg][6]][0][hits[kk][6]]++;
	  }
	  //	diagonal 1) string 2-3 
	  else if(hits[kk][4] + hits[gg][4] == 5 && TMath::Abs(hits[kk][4]-hits[gg][4])==1) {
	    distance =(-hits[kk][4]+hits[gg][4])*TimeSequenceFilter::getHorizontalHitDistance(kk,gg);
	    same_speed_d1[400 + TMath::Nint(15.0/timeStep*(hits[kk][0] - hits[gg][0])/distance) ]++;	
	    antuses[hits[gg][6]][1][hits[kk][6]]++;
	  }
	  //	lateral 2) string 2-4 and string 1-3 
	  else if( TMath::Abs(hits[kk][4]-hits[gg][4])==2) {
	    distance = (-hits[kk][4]+hits[gg][4])/2.0*TimeSequenceFilter::getHorizontalHitDistance(kk,gg);
	    same_speed_v[400 + TMath::Nint(15.0/timeStep*(hits[kk][0] - hits[gg][0])/distance) ]++;
	    antuses[hits[gg][6]][2][hits[kk][6]]++;
	  }
	  //	diagonal 2) string 1-4
	  else if(hits[kk][4] + hits[gg][4] == 5 && TMath::Abs(hits[kk][4]-hits[gg][4])==3) {
	    distance = (-hits[kk][4]+hits[gg][4])/3.0*TimeSequenceFilter::getHorizontalHitDistance(kk,gg);
	    same_speed_d2[400 + TMath::Nint(15.0/timeStep*(hits[kk][0] - hits[gg][0])/distance) ]++;
	    antuses[hits[gg][6]][3][hits[kk][6]]++;
	  }
	}
      }
      //And here we get to the vertical pairs:
      if(hits[kk][4]==hits[gg][4] && hits[kk][5]!=hits[gg][5] && (int)hits[kk][5]/2==(int)hits[gg][5]/2) {
	distance = hits[gg][3]-hits[kk][3];
	same_speed_h[400 + TMath::Nint(18.0/timeStep*(hits[kk][0] - hits[gg][0])/distance) ]++;
	antuses[hits[gg][6]][4][hits[kk][6]]++;
      }
    }
  }
  return count;
}
//
//
//
//
//
//
//
//
//
//
double TimeSequenceFilter::pattern_check(TH2D *trig_pat, double *pat_check_result, int cut_value, std::vector< std::vector< double > > antloc) {
  int binned_hit_count=0;

   for(int i=0; i<800; i++) {
     same_speed_t[i] = 0;
     same_speed_v[i] = 0;
     same_speed_d1[i] = 0;
     same_speed_d2[i] = 0;
     same_speed_h[i] = 0;
   }
  
  int c_s[800] = {0};
  int c_d[800] = {0};
  int speed_cut_value=0.0;
  
  
  int speed_t_max=0;
  int speed_v_max=0;
  int speed_h_max=0;
  int speed_d_max1=0;
  int speed_d_max2=0;
  int cs_max=0;
  int cd_max=0;
  
  int temp_speed_t_max=0;
  int temp_speed_v_max=0;
  int temp_speed_h_max=0;
  int temp_speed_d1_max=0;
  int temp_speed_d2_max=0;
  int temp_cs_max=0;
  int temp_cd_max=0;
  
  int max_speed_t_pos=-400;
  int max_speed_v_pos=-400;
  int max_speed_h_pos=-400;
  int max_speed_d1_pos=-400;
  int max_speed_d2_pos=-400;
  int max_cs_pos=-400;
  int max_cd_pos=-400;
  
  binned_hit_count = inter_string_check(trig_pat, &same_speed_t[0], &same_speed_d1[0], &same_speed_d2[0], &same_speed_v[0], &same_speed_h[0],&c_s[0], &c_d[0], antloc);
  
  binnedHitCount = binned_hit_count;
  
  speed_cut_value =same_speed_d1[416]+c_s[16] +same_speed_d1[415]+c_s[15] +  same_speed_d1[417]+c_s[17] ;

  int qpBinWidth = 1;
  for(int tt=qpBinWidth; tt<(800-qpBinWidth); tt++) {
    temp_speed_t_max = same_speed_t[tt];
    temp_speed_v_max = same_speed_v[tt]+same_speed_v[tt-qpBinWidth]+same_speed_v[tt+qpBinWidth];
    temp_speed_h_max = same_speed_h[tt]+same_speed_h[tt-qpBinWidth]+same_speed_h[tt+qpBinWidth];
    temp_speed_d1_max = same_speed_d1[tt]+same_speed_d1[tt-qpBinWidth]+same_speed_d1[tt+qpBinWidth];
    temp_speed_d2_max = same_speed_d2[tt]+same_speed_d2[tt-qpBinWidth]+same_speed_d2[tt+qpBinWidth];
    temp_cs_max = c_s[tt];
    temp_cd_max = c_d[tt];
    
    for(int ibin=1; ibin<=qpBinWidth; ibin++) {
      temp_speed_t_max += same_speed_t[tt-qpBinWidth]+same_speed_t[tt+qpBinWidth];
      temp_speed_v_max += same_speed_v[tt-qpBinWidth]+same_speed_v[tt+qpBinWidth];
      temp_speed_h_max += same_speed_h[tt-qpBinWidth]+same_speed_h[tt+qpBinWidth];
      temp_speed_d1_max += same_speed_d1[tt-qpBinWidth]+same_speed_d1[tt+qpBinWidth];
      temp_speed_d2_max += same_speed_d2[tt-qpBinWidth]+same_speed_d2[tt+qpBinWidth];
      temp_cs_max += c_s[tt-qpBinWidth] + c_s[tt+qpBinWidth];
      temp_cd_max += c_d[tt-qpBinWidth] + c_d[tt+qpBinWidth];
    }
    
    if(temp_speed_t_max>speed_t_max ) {speed_t_max = temp_speed_t_max; max_speed_t_pos = tt - 400;}
    if(temp_speed_v_max>speed_v_max) {speed_v_max = temp_speed_v_max; max_speed_v_pos = tt - 400;}
    if(temp_speed_h_max>speed_h_max) {speed_h_max = temp_speed_h_max; max_speed_h_pos = tt - 400;}
    if(temp_speed_d1_max>speed_d_max1) {speed_d_max1 = temp_speed_d1_max; max_speed_d1_pos = tt -400;}
    if(temp_speed_d2_max>speed_d_max2) {speed_d_max2 = temp_speed_d2_max; max_speed_d2_pos = tt - 400;}
    if(temp_cs_max>cs_max){ cs_max= temp_cs_max; max_cs_pos = tt;}
    if(temp_cd_max>cd_max){ cd_max= temp_cd_max; max_cd_pos = tt;}
  }
  double pos_side = TMath::Sqrt( max_speed_t_pos*max_speed_t_pos/15.0/15.0 + max_speed_v_pos*max_speed_v_pos/15.0/15.0 + max_speed_h_pos*max_speed_h_pos/18.0/18.0);
  double pos_diagonal = TMath::Sqrt(max_speed_d1_pos*max_speed_d1_pos/15.0/15.0 + max_speed_d2_pos*max_speed_d2_pos/15.0/15.0 + max_speed_h_pos*max_speed_h_pos/18./18.0);
  pat_check_result[0] = speed_cut_value;		// Time diff complete
  pat_check_result[1] = speed_t_max;
  pat_check_result[2] = speed_v_max;
  pat_check_result[3] = speed_h_max;
  pat_check_result[4] = speed_d_max1;
  pat_check_result[5] = speed_d_max2;
  pat_check_result[6] = pos_side; 
  pat_check_result[7] = pos_diagonal;
  pat_check_result[8] = max_speed_t_pos;
  pat_check_result[9] = max_speed_v_pos;
  pat_check_result[10] = max_speed_h_pos;
  pat_check_result[11] = max_speed_d1_pos;
  pat_check_result[12] = max_speed_d2_pos; 
  pat_check_result[13] = binned_hit_count;
  
  return pos_side;

}
//
//
//
//
//
//
//
//
//
//
///int gen_noise_std_dev(UsefulAtriStationEvent *theEvent, double *ch_noise_std_dev, double *ch_noise_average, int sum_time)
int TimeSequenceFilter::gen_noise_std_dev(std::vector<TGraph *> gr, double *ch_noise_std_dev, double *ch_noise_average, int sum_time) {
  Double_t *times = 0;
  Double_t *volts = 0;
  double en_average=0;
  double en_var=0;
  double en_std_dev=0;
  double en_sum1=0;
  double en_sum2=0;
  double en_sum3=0;
  double en_sum4=0;
  double en_sq_sum1=0;
  double en_sq_sum2=0;
  double en_sq_sum3=0;
  double en_sq_sum4=0;
  double min_sum1 =0;
  double min_sum2 =0;
  double min_sq_sum1 =0;
  double min_sq_sum2 =0;
  int en_bins = 0;
  int bin;
  for(int rfchan=0;rfchan<16;rfchan++) {
    times = gr[rfchan]->GetX();
    volts = gr[rfchan]->GetY();
    bin =gr[rfchan]->GetN() - 1; 
    double getx[bin];
    double gety[bin];
    double prop_energy[bin];
    en_average=0;
    en_var=0;
    en_std_dev=0;
    min_sum1 =0;
    min_sum2 =0;
    min_sq_sum1 =0;
    min_sq_sum2 =0;
    en_sum1=0;
    en_sq_sum1=0;
    en_sum2=0;
    en_sq_sum2=0;
    en_sum3=0;
    en_sq_sum3=0;
    en_sum4=0;
    en_sq_sum4=0;
    en_bins = 0;
    double en_time[bin];
    prop_energy[en_bins]=0;
    for (int l=0; l<bin; l++) {
      getx[l] = times[l]; 
      gety[l] = volts[l];
      if(l>sum_time) {
	for(int ie =0; ie<sum_time;ie++) {
	  prop_energy[en_bins]+=gety[l-sum_time+ie]*gety[l-sum_time+ie];
	}
      }
      prop_energy[en_bins] = TMath::Sqrt(prop_energy[en_bins]/sum_time);
      en_time[en_bins] = getx[l];
      
      if(l<bin/4)en_sum1+=prop_energy[en_bins];
      if(l>bin/4-1 && l<bin/2)en_sum2+=prop_energy[en_bins];
      if(l>bin/2-1 && l<bin/4*3.0)en_sum3+=prop_energy[en_bins];
      if(l>bin/4*3.0-1)en_sum4+=prop_energy[en_bins];
      
      if(l<bin/4)en_sq_sum1+=prop_energy[en_bins]*prop_energy[en_bins];
      if(l>bin/4-1 && l<bin/2)en_sq_sum2+=prop_energy[en_bins]*prop_energy[en_bins];
      if(l>bin/2-1 && l<bin/4*3.0)en_sq_sum3+=prop_energy[en_bins]*prop_energy[en_bins];
      if(l>bin/4*3.0-1)en_sq_sum4+=prop_energy[en_bins]*prop_energy[en_bins];
      en_bins++;
      prop_energy[en_bins] = 0;
    }

    en_sum1 = en_sum1/(en_bins/4.0 - sum_time);
    en_sum2 = en_sum2/(en_bins/4.0 );
    en_sum3 = en_sum3/(en_bins/4.0 );
    en_sum4 = en_sum4/(en_bins/4.0 );
    en_sq_sum1 = en_sq_sum1/(en_bins/4.0 - sum_time);
    en_sq_sum2 = en_sq_sum2/(en_bins/4.0 );
    en_sq_sum3 = en_sq_sum3/(en_bins/4.0 );
    en_sq_sum4 = en_sq_sum4/(en_bins/4.0 );
    
    min_sum1 = en_sum1;
    if(en_sum2<min_sum1) {min_sum1 = en_sum2; min_sum2 = en_sum1;}
    else min_sum2 = en_sum2;
    if(en_sum3<min_sum1) {min_sum2 = min_sum1; min_sum1 = en_sum3;}
    else if(en_sum3<min_sum2) {min_sum2 = en_sum3;}
    if(en_sum4<min_sum1) {min_sum2 = min_sum1; min_sum1 = en_sum4;}
    else if(en_sum4<min_sum2) {min_sum2 = en_sum4;}
    
    min_sq_sum1 = en_sq_sum1;
    if(en_sq_sum2<min_sq_sum1) {min_sq_sum1 = en_sq_sum2; min_sq_sum2 = en_sq_sum1;}
    else min_sq_sum2 = en_sq_sum2;
    if(en_sq_sum3<min_sq_sum1) {min_sq_sum2 = min_sq_sum1; min_sq_sum1 = en_sq_sum3;}
    else if(en_sq_sum3<min_sq_sum2) {min_sq_sum2 = en_sq_sum3;}
    if(en_sq_sum4<min_sq_sum1) {min_sq_sum2 = min_sq_sum1; min_sq_sum1 = en_sq_sum4;}
    else if(en_sq_sum4<min_sq_sum2) {min_sq_sum2 = en_sq_sum4;}
    
    en_average = (double)(min_sum1 + min_sum2)/2.0;
    en_var = 1.0*(min_sq_sum1 + min_sq_sum2)/2.0 - en_average*en_average;
    en_std_dev = TMath::Sqrt(en_var);
    ch_noise_std_dev[rfchan] = en_std_dev;
    ch_noise_average[rfchan] = en_average;
    chNoiseEnergyRMS[rfchan] = ch_noise_std_dev[rfchan];  //FIXME: modified
  }
  return en_average;
} 
//
//
//
//
//
//
//
//
//
//
double TimeSequenceFilter::gen_trig_histo_single(std::vector<TGraph *> gr, TH2D *trig_pat, double *ch_std_dev, double *ch_average, double th_factor, int sum_time, int stationId) {
  Double_t *times = 0;
  Double_t *volts = 0;
  int bin = 0;
  int hit_count =0;
  int rfchan =0;
  int en_bins = 0;
  for(int rfchan=0;rfchan<16;rfchan++) {//loop channels
    bool hitb = false; // My Add
    times = gr[rfchan]->GetX();
    volts = gr[rfchan]->GetY();
    bin = gr[rfchan]->GetN()-1;
    double getx[bin];
    double gety[bin];
    double prop_energy[bin];
    double en_time[bin];
    en_bins = 0;
    prop_energy[en_bins] = 0;
    for (int l=0; l<bin; l++) {
	getx[l] = times[l];
	gety[l] = volts[l];
	if(l>sum_time) {
	  for(int ie =0; ie<sum_time;ie++) {
	    prop_energy[en_bins]+=gety[l-sum_time+ie]*gety[l-sum_time+ie];
	  }
	}
	prop_energy[en_bins] = TMath::Sqrt(prop_energy[en_bins]/sum_time); //division is wrong, doesn't matter - changed to have the energy outstanding from the waveform in the plot
	en_time[en_bins] = getx[l];
	en_bins++;
	prop_energy[en_bins] = 0;
    }	
    for(int tr=0;tr<en_bins;tr++) {//Fill the trigger pattern histogram
      if(prop_energy[tr]>th_factor*ch_std_dev[rfchan]+ch_average[rfchan] && !(stationId==2 && (rfchan/4==3 && rfchan%4==3))) {
	trig_pat->Fill(en_time[tr],rfchan);
	if(!hitb) { hit_count++; hitb = true; }
      }
    }
  }//end loop channels
  return hit_count;
}
//
//
//
//
//
//
//
//
//
//
double TimeSequenceFilter::getQualityParameter(std::vector<TGraph*> gr_in, std::vector<std::vector<double> > ant_loc, int stationId, double *qualArray) {

  double quality = 0;
  double th_factor = 4.0;// atof(argv[2])/10.0;
  int sum_int_time = 10;//atoi(argv[3]);
  int pattern_bin_no = 500/(sum_int_time/2.0);
  std::vector<TGraph *> gr;
  
  for(int u=0;u<16;u++) {
    ch_noise_std_dev[u] = 0;
    ch_noise_average[u] = 0;
  }
  
  int average_count = 0;
  
  
  for(int a=0;a<16;a++) {
    gr.push_back( FFTtools::getInterpolatedGraph(gr_in[a],0.5));
  }
  double timeLength = 0;
  timeLength = gr[0]->GetN()/2.0;			//this is the duration of the waveform
  blockCount = timeLength/20.0*4.0;	//actually the count of blocks. Factor 4 is historical. Just a scaling...
  
  double pat_check_results[20];
  int cut_value = 0;
  int hist_hit_count =0;
  double ch_energy[16] = {0};
  
  TH2D * trig_pat = new TH2D("trig_pat", "The event trigger pattern", pattern_bin_no,0,500, 16,0,16);
  average_count = gen_noise_std_dev(gr, &ch_noise_std_dev[0], &ch_noise_average[0],sum_int_time);
  getMaxEnergy(gr);
  hist_hit_count = gen_trig_histo_single(gr, trig_pat, &ch_noise_std_dev[0], &ch_noise_average[0], th_factor,sum_int_time, stationId);
  pattern_check(trig_pat, &pat_check_results[0], cut_value, ant_loc);
  double side_speed = 1.0/(pat_check_results[6]*4.0)*10.0;
  double diagonal_speed = 1.0/(pat_check_results[7]*4.0)*10.0;
  double theta = TMath::ACos( pat_check_results[10]/18.0*(sum_int_time/2.0)*0.171 );
  totalHitCount = hist_hit_count;
  quality = (pat_check_results[1] + pat_check_results[2] + pat_check_results[3])*1.0/blockCount;
		
  qualArray[0] = (pat_check_results[1] + pat_check_results[2] + pat_check_results[3])*1.0/blockCount;
  qualArray[1] = (side_speed)*(1.0);
  qualArray[2] = (diagonal_speed)*(1.0);
  qualArray[3] = theta;
  qualArray[4] = totalHitCount;
  

  qualityParameter = qualArray[0];
  
  delete trig_pat;
  for(int i=0;i<16;i++){
    delete gr[i];
  }
  
  return quality;

}
//
//
//
//
//
//
//
//
//
//
double TimeSequenceFilter::getMaxEnergy(std::vector<TGraph *> gr) {
  int bin = 0;
  int hit_count =0;
  int rfchan =0;
  int en_bins = 0;
  int sum_time = 80;
  double tempEnergy = 0;
  double tempPosition = 0;
  std::vector<double> chEnvelope;
  std::vector<double> chTime;
  for(int rfchan=0;rfchan<16;rfchan++) {//loop channels
    tempEnergy = 0;
    bin = gr[rfchan]->GetN();
    double getx[bin];
    double gety[bin];
    double prop_energy[bin];
    double en_time[bin];
    en_bins = 0;
    prop_energy[en_bins] = 0;
    for (int l=0; l<bin; l++) {
      getx[l] = gr[rfchan]->GetX()[l];
      gety[l] = gr[rfchan]->GetY()[l];
      if(l>sum_time){
	for(int ie =0; ie<sum_time;ie++){
	  prop_energy[en_bins]+=gety[l-sum_time+ie]*gety[l-sum_time+ie];
	}
      }
      prop_energy[en_bins] = TMath::Sqrt(prop_energy[en_bins]/sum_time); //division is wrong, doesn't matter - changed to have the energy outstanding from the waveform in the plot
      en_time[en_bins] = getx[l];
      chEnvelope.push_back(prop_energy[en_bins]);
      chTime.push_back(getx[l-sum_time/2]);
      en_bins++;
      prop_energy[en_bins] = 0;
    }
    chEnvelope.clear();	
    chTime.clear();	
    for(int tr=0;tr<en_bins;tr++) {//Fill the trigger pattern histogram
      if(prop_energy[tr]>tempEnergy && en_time[tr]<310.0 ) {
	tempEnergy = prop_energy[tr];
	tempPosition = en_time[tr];
      }
    }
    
    if(tempPosition - en_time[0] > sum_time + 100 ) chNoiseEnergy[rfchan] = prop_energy[sum_time+10];
    else chNoiseEnergy[rfchan] = prop_energy[en_bins-20]; //20ns before the end -->this is the dangerous part: Signal could be very long
    
    chSignalEnergy[rfchan] = tempEnergy;
    chSignalMaxEnergyPosition[rfchan] = tempPosition;
  }//end loop channels

  return hit_count;

}
