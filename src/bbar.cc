
#include "chrono"
#include "TreeManager.h"
#include "Util.h"
#include "TTree.h"
#include "TFile.h"

#include "Pythia8/Pythia.h"


using namespace Pythia8;
int main(int argc, char* argv[]) {

  // Check arguments.
  if (argc < 2) {
    cerr << " Unexpected number of command-line arguments. \n You are"
         << " expected to provide the arguments \n"
         << " 1. CfgFile \n"
         << " Program stopped. " << endl;
    return 1;
  }
  
  std::map<TString, Float_t > storageFloat;

  string Pythia8Data(getenv("PYTHIA8DATA"));
  std::cout<<"Setting default PYTHIA8DATA as "<<Pythia8Data<<"\n";
  string PythiaCfgFile("data/pythia.config");
  string outputFname("output.root");
  string ChannelDesc("Pythia GEN");
  Long64_t maxEvents(500);
  Int_t reportEvery(100);
  TString tag;

 /***********************************************************************************************************/
    std::cout<<"Reading the config file  : "<<argv[1]<<"\n";
    fstream cfgFile(argv[1],ios::in);
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
           getStringFromTag("PythiaCfgFile"   ,strStream , field , PythiaCfgFile);
           getStringFromTag("Pythia8Data"     ,strStream , field , Pythia8Data);
           getStringFromTag("ChannelDesc"     ,strStream , field , ChannelDesc);
           getStringFromTag("OutputFname"     ,strStream , field , outputFname);
       }
    }
    cfgFile.close();

   /***********************************************************************************************************/


  // Intialize Pythia.
  Pythia pythia;
  pythia.readFile(PythiaCfgFile);
  pythia.init();

  setenv("PYTHIA8DATA", Pythia8Data.c_str(), 1);

  // Set up the ROOT TFile and TTree.
  TFile *file = TFile::Open(outputFname.c_str(),"recreate");
  Event *eventP = &pythia.event;
  
  TreeManager treeManager;
  TTree *T = new TTree("events",ChannelDesc.c_str());
  
  // Setup Particles To Store
  treeManager.AddParticleBranch(T,"g1");
  treeManager.AddParticleBranch(T,"g2");
  treeManager.AddParticleBranch(T,"b1");
  treeManager.AddParticleBranch(T,"b2");
  treeManager.AddParticleBranch(T,"q1");
  treeManager.AddParticleBranch(T,"q2");

    tag = "bbar_mass"    ;  storageFloat[tag] = 0.0; T->Branch(tag,&storageFloat[tag]) ; 

  // The event loop.
  Hist mass("pT(bbar) [GeV]", 100, 0., 2.);
  Float_t dTemp;
  auto t_start = std::chrono::high_resolution_clock::now();
  auto t_end = std::chrono::high_resolution_clock::now();
  for (int iEvent = 0; iEvent < maxEvents; ++iEvent) {
 
       if(iEvent%reportEvery == 0 )
       {
          t_end = std::chrono::high_resolution_clock::now();
          std::cout<<"Processing Entry in event loop : "<<iEvent<<" / "<<maxEvents<<"  [ "<<100.0*iEvent/maxEvents<<"  % ]  "
             << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
             <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents - iEvent)/(1e-9 + iEvent)* 0.001
             <<std::endl;
       }

    // Generate the event.
    if (!pythia.next()) continue;
    // Analyze the event.
    Event &event = pythia.event;
    std::cout<<"Event : "<<iEvent<<" : \n";
    for (int iPrt = 0; iPrt < (int)event.size(); ++iPrt) {
      
      if (iPrt > 7) break;
      std::cout<<"\t i = "<<iPrt<<" | "<<event[iPrt].id()<<" | "<<" parents : "<<event[iPrt].mother1()<<" , "<<event[iPrt].mother2()<<" | daughters : "<<event[iPrt].daughter1()<<" , "<<event[iPrt].daughter2()<<" pT : "<<event[iPrt].pT()<<"\n";
    
    }
    std::cout<<"\t g1 : "<<event[3].id()<<" pT : "<<event[3].pT()<<"\n";
    std::cout<<"\t g2 : "<<event[4].id()<<" pT : "<<event[4].pT()<<"\n";
    std::cout<<"\t g1 + g2  pT : "<<(event[3].p() + event[4].p() ).pT() <<"\n";
    mass.fill((event[3].p() + event[4].p()).pT());
    T->Fill();
  }

  // Print the statistics and histogram.
  pythia.stat();
  mass /= mass.getEntries();
  cout << mass;
  
  T->Write();
  file->Write();
  file->Close();
  return 0;
}
