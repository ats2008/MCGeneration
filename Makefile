CXX 		 =g++
LCC_RELEASE  =LCG_104
PLATFORM     =x86_64-centos7-gcc11-opt


ROOT_LIBS	 =`root-config --libs` 
PYTHIA_LIBS  =`pythia8-config --cxxflags --ldflags --libs`
LIBS_EXTRA_BASE     = -L/cvmfs/sft.cern.ch/lcg/views/$(LCC_RELEASE)/$(PLATFORM)/lib  -L/cvmfs/sft.cern.ch/lcg/views/$(LCC_RELEASE)/$(PLATFORM)/lib64  
LIBS_EXTRA_EVTGEN   = -lEvtGen -lEvtGenExternal  -lgfortran -lPhotospp 
LIBS_EXTRA_HEPMC    =  -lHepMC3 

LIBS_EXTRA = ${LIBS_EXTRA_BASE} ${LIBS_EXTRA_EVTGEN} ${LIBS_EXTRA_HEPMC}

CXXFLAGS 	=-g3 -std=c++17  $(ROOT_LIBS) $(PYTHIA_LIBS) $(LIBS_EXTRA_BASE) 
INCPATH     = include -I/cvmfs/sft.cern.ch/lcg/views/$(LCC_RELEASE)/$(PLATFORM)/include

TARGET		= evtGenExample pythiaExample
all:	$(TARGET) 

showerLHE : src/showerLHE.cc  
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS) ${LIBS_EXTRA_HEPMC}

showerLHE_lepVeto : src/showerLHE_lepVeto.cc
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS) ${LIBS_EXTRA_HEPMC}

evtGenExample : src/example1.cc  
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

pythiaExample : src/main48.cc
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

bbar : src/bbar.cc
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

bbarHEPMC : src/bbarHEPMC.cc
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

Jpsi2MuMuGamma : src/jpsiToMuMuGamma.cc
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o bin/$@.exe $<  $(LIBS)

bsToMuMuGamma : src/bsToMuMuGamma.cc
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o bin/$@.exe $<  $(LIBS)

b0ToDstarSMinus : src/b0ToDstarSMinus.cc
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o bin/$@.exe $<  $(LIBS)

b0ToDstarSMinusFullAcceptance : src/b0ToDstarSMinusFullAcceptance.cc
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o bin/$@.exe $<  $(LIBS)

bsToDstarSMinus : src/bsToDstarSMinus.cc
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o bin/$@.exe $<  $(LIBS)

bsToDstarSMinusFullAcceptance : src/bsToDstarSMinusFullAcceptance.cc
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o bin/$@.exe $<  $(LIBS)

clean :
	@rm *.exe

