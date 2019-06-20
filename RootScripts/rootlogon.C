{
  cout << "Executing rootlogon.C" << endl;

  // Set up includes to AraEvent
  gSystem->AddIncludePath("-I$ARA_UTIL_INSTALL_DIR/include/");
  // Load the AraEvent library
  gSystem->Load("$ARA_UTIL_INSTALL_DIR/lib/libAraEvent.so");
  // Add filter includes and load ASFilters shared library
  gSystem->AddIncludePath("-I$AS_FILTERS/filters/");
  gSystem->Load("$AS_FILTERS/filters/libASFilters.so");

}
