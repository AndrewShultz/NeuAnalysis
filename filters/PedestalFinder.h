
#ifndef PEDESTALFINDER_H
#define PEDESTALFINDER_H

#include <string>
#include <fstream>

using namespace std;

string findPedestalFile(int stationID, int runNum) {

  string line;
  ifstream ifs(Form("/data/user/aschultz/PedestalHandling/goodPedestalsARA0%i.txt",stationID));

  string pedfile = "";

  while(getline(ifs,line)) {

    std::string::size_type pos1 = line.find("run");
    std::string::size_type pos2 = line.find(".dat");

    if( pos1==std::string::npos || pos2==std::string::npos) continue;

    string sustr = line.substr(pos1+3,6);
    int intRunNum = atoi(sustr.c_str());

    if(intRunNum>runNum) break;

    pedfile = line;

  }

  return pedfile;

}

#endif
