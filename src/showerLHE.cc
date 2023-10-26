// main11.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; LHE file;

// This is a simple test program.
// It illustrates how Les Houches Event File input can be used in Pythia8.
// It uses the ttsample.lhe(.gz) input file, the latter only with 100 events.

#include "chrono"
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"
#include "Util.h"
#include "Pythia8Plugins/ResonanceDecayFilterHook.h"


using namespace std;
using namespace Pythia8;
int main( int argc, char* argv[]  ) {

  // You can always read an plain LHE file,
  // but if you ran "./configure --with-gzip" before "make"
  // then you can also read a gzipped LHE file.
  bool useGzip = false;

  // Generator. We here stick with default values, but changes
  // could be inserted with readString or readFile.
  Pythia pythia;
  
  // Setup UserHooks
  auto myUserHooks = make_shared<ResonanceDecayFilterHook>(pythia.settings);
  pythia.setUserHooksPtr( myUserHooks);

  string Pythia8Data(getenv("PYTHIA8DATA"));
  std::cout<<"Setting default PYTHIA8DATA as "<<Pythia8Data<<"\n";
  string PythiaCfgFile("data/pythia.config");
  string inputFname("ipt.root");
  string outputFname("output.root");
  string ChannelDesc("Pythia GEN");
  Long64_t maxEvents(500);
  Int_t reportEvery(100);
  TString tag;
  Int_t cmdIDX=1;

  // pythia.readString("Print:verbosity = 3");
  // Initialize Les Houches Event File run. List initialization information.
  pythia.readString("Beams:frameType = 4");
  //if (useGzip) pythia.readString("Beams:LHEF = ttbar.lhe.gz");
  //else         pythia.readString("Beams:LHEF = ttbar.lhe");
   
   /***********************************************************************************************************/
    std::cout<<"Reading the config file  : "<<argv[cmdIDX]<<"\n";
    fstream cfgFile(argv[cmdIDX],ios::in); cmdIDX+=1 ;
	string line;
	bool cfgModeFlag=false;

    Double_t aDouble;
    cfgFile.clear();
    cfgFile.seekg(0,ios::beg);
    cfgModeFlag=false;
    std::istringstream strStream;
    std::string field;
    Int_t tmpI;

	while(std::getline(cfgFile,line))
	{
	   if(line=="#PARAMS_BEG") {cfgModeFlag=true;continue;}
	   if(line=="#PARAMS_END") {cfgModeFlag=false;continue;}
	   if(not cfgModeFlag) continue;
	   if(line=="") continue;
	   if(line=="#END") break;
       strStream.clear();
       strStream.str(line);
       while (getline(strStream, field,'='))
       {
           getIntFromTag("ReportEvery"        ,strStream , field , reportEvery);
           getIntFromTag("MaxEvents"          ,strStream , field , maxEvents);
           getStringFromTag("PythiaCfgFile"        ,strStream , field , PythiaCfgFile);
           getStringFromTag("Pythia8Data"        ,strStream , field , Pythia8Data);
           getStringFromTag("ChannelDesc"        ,strStream , field ,ChannelDesc);
           getStringFromTag("OutputFname"        ,strStream , field ,outputFname);
           getStringFromTag("InputLHEFile"    ,strStream , field , inputFname);
       }
    }
    cfgFile.close();

   /***********************************************************************************************************/
   if ( argc > ( cmdIDX ) ) {    
        PythiaCfgFile=argv[cmdIDX];   
        std::cout<<"[ "<<cmdIDX<<" ] "<<"** Setting input pythia config file name as : "<<PythiaCfgFile<<"\n";
   } ; cmdIDX+=1;
   
   if ( argc > ( cmdIDX ) ) {    
        inputFname=argv[cmdIDX];   
        std::cout<<"[ "<<cmdIDX<<" ] "<<"** Setting input file name as : "<<inputFname<<"\n";
   } ; cmdIDX+=1;
  
  if (argc > (cmdIDX) ) {    
        outputFname=argv[cmdIDX]; 
        std::cout<<"[ "<<cmdIDX<<" ] "<<"** Setting output file name as : "<<outputFname<<"\n";
   } ; cmdIDX+=1;
    
   std::cout<<" ReportEvery  : "<<reportEvery<<"\n";
   std::cout<<" MaxEvents    : "<<maxEvents<<"\n";
   std::cout<<" Pythia8Data  : "<<Pythia8Data<<"\n";
   std::cout<<" PythiaCfgFile  : "<<PythiaCfgFile<<"\n";
   std::cout<<" ChannelDesc  : "<<ChannelDesc<<"\n";
   std::cout<<" InputLHEFile  : "<<inputFname<<"\n";
   
   pythia.readFile(PythiaCfgFile);
   pythia.readString("Beams:LHEF = "+inputFname);
   pythia.init();

  // Book HEPMC Output
  Pythia8::Pythia8ToHepMC topHepMC(outputFname);
  
  // Book histogram.
  Hist nCharged("charged particle multiplicity",100,-0.5,399.5);

  // Allow for possibility of a few faulty events.
  int nAbort = 10;
  int iAbort = 0;

  // Begin event loop; generate until none left in input file.
  auto t_start = std::chrono::high_resolution_clock::now();
  auto t_end = std::chrono::high_resolution_clock::now();

  for (int iEvent = 0; ; ++iEvent) {
    if ( (maxEvents > 0 ) and ( iEvent >= maxEvents ) ) break;
     
    if(iEvent%reportEvery == 0 )
       {
          t_end = std::chrono::high_resolution_clock::now();
          std::cout<<"Processing Entry in event loop : "<<iEvent<<" / "<<maxEvents<<"  [ "<<100.0*iEvent/maxEvents<<"  % ]  "
             << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
             <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents - iEvent)/(1e-9 + iEvent)* 0.001
             <<std::endl;
       }



    // Generate events, and check whether generation failed.
    if (!pythia.next()) {
      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;
      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }

    // Sum up final charged multiplicity and fill in histogram.
        int nChg = 0;
        for (int i = 0; i < pythia.event.size(); ++i)
            if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
                ++nChg;
        nCharged.fill(nChg);
    

    // Regitering the event to HEPMC output
    topHepMC.writeNextEvent( pythia );

  // End of event loop.
  }

  // Give statistics. Print histogram.
  pythia.stat();
  cout << nCharged;

  // Done.
  return 0;
}
