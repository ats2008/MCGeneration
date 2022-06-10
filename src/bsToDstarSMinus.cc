
#include "chrono"
#include "TreeManager.h"
#include "Util.h"
#include "TTree.h"
#include "TFile.h"

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/EvtGen.h"

#define HEPMC2

#ifndef HEPMC2
#include "Pythia8Plugins/HepMC3.h"
#else
#include "Pythia8Plugins/HepMC2.h"
#endif



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
  string EvtGenDecayFile("data/DECAY.DEC");
  string EvtGenParicletDataFile("data/evt.pdl");
  string ParticleDacayFile("data/customDecays/B02DxMuChain.DEC");
  string outputFname("output.root");
  string ChannelDesc("Pythia GEN");
  string hepMCOutputName("outputHepMC.lhe");
  vector<string> pythiaCustomization;
  Bool_t SaveLHE(false);
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
           getStringFromTag("PythiaCfgFile"        ,strStream , field , PythiaCfgFile);
           getStringFromTag("Pythia8Data"        ,strStream , field , Pythia8Data);
           getStringFromTag("EvtGenDecayFile"        ,strStream , field ,EvtGenDecayFile);
           getStringFromTag("EvtGenParicletDataFile"        ,strStream , field ,EvtGenParicletDataFile);
           getStringFromTag("ParticleDacayFile"        ,strStream , field , ParticleDacayFile);
           getStringFromTag("ChannelDesc"        ,strStream , field ,ChannelDesc);
           getStringFromTag("OutputFname"        ,strStream , field ,outputFname);
           getBoolFromTag("SaveLHE"            ,strStream , field ,SaveLHE);
           getStringFromTag("OutputHEPFname"     ,strStream , field ,hepMCOutputName);
       }
    }
    getVetorFillledFromConfigFile(cfgFile , pythiaCustomization, "#PYTHIA_STR_BEG","#PYTHIA_STR_END",true);
    cfgFile.close();

   /***********************************************************************************************************/


  // Intialize Pythia.
  Pythia pythia;
  pythia.readFile(PythiaCfgFile);
  for(auto customization : pythiaCustomization)
  {
    pythia.readString(customization);
  }
  pythia.init();
  
  Pythia8::Pythia8ToHepMC saveToHepMC;
  if(SaveLHE) saveToHepMC.setNewFile(hepMCOutputName);

  // Initialize EvtGen.
  EvtGenDecays *evtgen = 0;
  setenv("PYTHIA8DATA", Pythia8Data.c_str(), 1);
  evtgen = new EvtGenDecays(&pythia, EvtGenDecayFile,EvtGenParicletDataFile);
  evtgen->readDecayFile(ParticleDacayFile);

  // Set up the ROOT TFile and TTree.
  TFile *file = TFile::Open(outputFname.c_str(),"recreate");
  Event *eventP = &pythia.event;
  
  TreeManager treeManager;
  TTree *T = new TTree("events",ChannelDesc.c_str());
  
  // Setup Particles To Store
  treeManager.AddParticleBranch(T,"Bs");
  treeManager.AddParticleBranch(T,"Ds_Star1");
  treeManager.AddParticleBranch(T,"Mu1");
  treeManager.AddParticleBranch(T,"Nu1");
  treeManager.AddParticleBranch(T,"Gamma");
  treeManager.AddParticleBranch(T,"Ds");
  treeManager.AddParticleBranch(T,"Mu2");
  treeManager.AddParticleBranch(T,"Nu2");

    tag = "Mu1Mu2_deltaR"  ;  storageFloat[tag] = 0.0; T->Branch(tag,&storageFloat[tag]) ; 
    tag = "Mu1Mu2_deltaEta";  storageFloat[tag] = 0.0; T->Branch(tag,&storageFloat[tag]) ; 
    tag = "Mu1Mu2_deltaPhi";  storageFloat[tag] = 0.0; T->Branch(tag,&storageFloat[tag]) ; 
    tag = "Mu1Mu2_mass"    ;  storageFloat[tag] = 0.0; T->Branch(tag,&storageFloat[tag]) ; 
    tag = "Mu1Gamma_mass"    ;  storageFloat[tag] = 0.0; T->Branch(tag,&storageFloat[tag]) ; 
    tag = "Mu2Gamma_mass"    ;  storageFloat[tag] = 0.0; T->Branch(tag,&storageFloat[tag]) ; 
    tag = "Mu1Mu2Gamma_mass" ;  storageFloat[tag] = 0.0; T->Branch(tag,&storageFloat[tag]) ; 

  // The event loop.
  Hist mass("m(e, pi0) [GeV]", 100, 0., 6.);
   Float_t dTemp;
  auto t_start = std::chrono::high_resolution_clock::now();
  auto t_end = std::chrono::high_resolution_clock::now();
  Long64_t nEvents(0);
  Long64_t nParent(0);
  Long64_t nDecays(0);
  Bool_t hasDecay(false);
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
    nEvents++;
    // Perform the decays with EvtGen.
    if (evtgen) evtgen->decay();
    // Analyze the event.
    Event &event = pythia.event;

    hasDecay=false;
    for (int iPrt = 0; iPrt < (int)event.size(); ++iPrt) {
      
      if (event[iPrt].idAbs() != 531) continue;
      nParent++;
      int iBs(event[iPrt].iBotCopyId());
      int iDsStar1(-1);
      int iDs(-1);
      int iMu1(-1);
      int iNuMu1(-1);
      int iMu2(-1);
      int iNuMu2(-1);
      int iG1(-1);

      //std::cout<<"Bs Meson Found : "<<event[iPrt].idAbs()<<"\n";
      for (int iDtr = event[iBs].daughter1(); iDtr <= event[iBs].daughter2(); ++ iDtr) {
        if (event[iDtr].idAbs() ==  13 )  iMu1    = event[iDtr].iBotCopyId();
        if (event[iDtr].idAbs() ==  14 )  iNuMu1  = event[iDtr].iBotCopyId();
        if (event[iDtr].idAbs() ==  433)  iDsStar1= event[iDtr].iBotCopyId();
      }
      //std::cout<<"\t\t\t iMu1 , iNu1 , iDs*1 "<<iMu1<<" , "<<iNuMu1<<" , "<<iDsStar1<<"\n"; 
      
      if(iBs ==-1 ) continue;
      if( iDsStar1  ==-1 ) continue;
      if( iMu1  ==-1 ) continue;
      if( iNuMu1  ==-1 ) continue;
        
      if(event[iMu1].pT() < 3.5 ) continue;
      if(event[iMu1].eta() < -2.5 ) continue;
      if(event[iMu1].eta() >  2.5 ) continue;
      
      for (int iDtr = event[iDsStar1].daughter1(); iDtr <= event[iDsStar1].daughter2(); ++ iDtr) {
        if (event[iDtr].idAbs() ==  431)  iDs   = event[iDtr].iBotCopyId();
        if (event[iDtr].idAbs() ==  22 )  iG1  = event[iDtr].iBotCopyId();
      }
      //std::cout<<"\t\t\t iDs , iG1 "<<iDs<<" , "<<iG1<<"\n"; 
      if( iDs  ==-1 ) continue;
      if( iG1  ==-1 ) continue;
      if(event[iG1].pT() < 3.5 ) continue;
      if(event[iG1].eta() < -1e4 ) continue;
      if(event[iG1].eta() > 3.0 ) continue;


      for (int iDtr = event[iDs].daughter1(); iDtr <= event[iDs].daughter2(); ++ iDtr) {
        if (event[iDtr].idAbs() ==  13 )  iMu2    = event[iDtr].iBotCopyId();
        if (event[iDtr].idAbs() ==  14 )  iNuMu2  = event[iDtr].iBotCopyId();
      }

      //std::cout<<"\t\t\t iMu2 , iNu2 "<<iMu2<<" , "<<iNuMu2<<"\n"; 
      if( iMu2  ==-1   ) continue;
      if( iNuMu2  ==-1 ) continue;
      if(event[iMu2].pT() < 3.5 ) continue;
      if(event[iMu2].eta() < -2.5 ) continue;
      if(event[iMu2].eta() >  2.5 ) continue;
      //std::cout<<"\t  mu2 pt"<<event[iMu2].pT() <<" /  3.5 "<<"\n";
      

      if(event[iMu1].pT() < event[iMu2].pT())
      {
        auto tmp=iMu1; iMu1=iMu2 ; iMu2=tmp;
      }
      
      nDecays++;
      storageFloat["Mu1Mu2_deltaR"  ] = getDR(event[iMu1].eta(),event[iMu1].phi(),event[iMu2].eta(),event[iMu2].phi());
      storageFloat["Mu1Mu2_deltaEta"] = getDETA(event[iMu1].eta(),event[iMu2].eta());
      storageFloat["Mu1Mu2_deltaPhi"] = getDPHI(event[iMu1].phi(),event[iMu2].phi());
      storageFloat["Mu1Mu2_mass"    ] = (event[iMu1].p() + event[iMu2].p()).mCalc() ;
      storageFloat["Mu1Gamma_mass"    ] = (event[iMu1].p() + event[iG1].p()).mCalc() ;
      storageFloat["Mu2Gamma_mass"    ] = (event[iMu2].p() + event[iG1].p()).mCalc() ;
      storageFloat["Mu1Mu2Gamma_mass" ] = (event[iMu1].p() +  event[iMu2].p() + event[iG1].p()).mCalc() ;

      mass.fill((event[iMu1].p() + event[iMu2].p()).mCalc());
  
      treeManager.FillParticleBranches("Bs",event,iBs);
      treeManager.FillParticleBranches("Ds_Star1",event,iDsStar1);
      treeManager.FillParticleBranches("Gamma",event,iG1);
      treeManager.FillParticleBranches("Ds",event,iDs);
      treeManager.FillParticleBranches("Mu1",event,iMu1);
      treeManager.FillParticleBranches("Nu1",event,iMu2);
      treeManager.FillParticleBranches("Mu2",event,iMu1);
      treeManager.FillParticleBranches("Nu2",event,iMu2);
      
      hasDecay=true;
    }
    if(hasDecay) T->Fill();
    if(SaveLHE) saveToHepMC.writeNextEvent( pythia );
  }

  // Print the statistics and histogram.
  pythia.stat();
  mass /= mass.getEntries();
  cout << mass;
  
  std::cout<<"         ==========     Summary    =========== "<<"\n"
           <<"            N Events Made : "<<nEvents          <<"\n"
           <<"            N Parent Made : "<<nParent          <<"\n"
           <<"      N Events With Decay : "<<nDecays          <<"\n"
           <<"         ===================================== "<<"\n";

  T->Write();
  file->Write();
  file->Close();
  if (evtgen) delete evtgen;
  return 0;
}
