#ifndef LOADBAR_H
#define LOADBAR_H

//#include <sys/stat.h>
#include <sys/ioctl.h>
//#include <unistd.h>
//#include <string>
#include <iostream>
//#include <fstream>
//#include <stdlib.h>
//#include <stdio.h>

using namespace std;

struct winsize size;
bool initalized = false;

void initLoadBar();

void loadBar(double current, double end,string outputString="");
void loadBar(int current, int end,string outputString="");
void loadBar(double percent,string outputString="");
void loadBarEnd();

void initLoadBar() {
  ioctl(STDOUT_FILENO,TIOCGWINSZ,&size);
  //printf("Columns = %i \n",size.ws_col);
  printf("\n\n\n");
  initalized = true;
}

void loadBar(double percent,string outputString) {

  if(!initalized) initLoadBar();

  if(isatty(1)) {
    printf("\033[1A\033[2K");
    if(outputString!="") std::cout << outputString << std::endl;
    if(percent!=100.0) {
      //printf("\033[2A\033[K",size.ws_row);
      int greenspace = (size.ws_col*0.01*percent);
      int redspace = size.ws_col-greenspace;
      printf("\033[42m");
      for(int i=0; i<greenspace; i++) printf(" ");
      printf("\033[41m");
      for(int i=0; i<redspace; i++) printf(" ");
      printf("\033[0m\n");
    }
  }

}

void loadBar(double current, double end,string outputString) {
  loadBar(current/end*100.0,outputString);
}

void loadBar(int current, int end,string outputString) {
  loadBar(double(current)/double(end)*100.0,outputString);
}

void loadBarEnd() {
  loadBar(1,1);
}

#endif
