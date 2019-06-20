{
  cout << "Executing rootlogon.C" << endl;

  // Set up includes to AraEvent
  gSystem->AddIncludePath("-I$ARA_UTIL_INSTALL_DIR/include/");
  // Load the AraEvent library
  gSystem->Load("$ARA_UTIL_INSTALL_DIR/lib/libAraEvent.so");

  gSystem->AddIncludePath("-I/home/aschultz/Analysis/filters/");
  gSystem->Load("/home/aschultz/Analysis/filters/libASFilters.so");

}
