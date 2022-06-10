#ifndef __TREEMANAGER_H_
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/EvtGen.h"


using namespace std;
using namespace Pythia8;

class TreeManager 
{

 public :
    std::map<TString, Float_t > storageFloat;
    
    void AddParticleBranch(TTree* aTree,TString tag);
    void ClearAllVar();
    void FillParticleBranches(TString tag,Event &event,Int_t index);

};

#include "../src/TreeManager.cc"

#endif
