MC Generation using standalone and cross compiled generator binaries from **LCG releases from cvmfs**

Soucrce the LCG Release
```
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_101 x86_64-centos7-gcc11-opt 
```
For Photos++
```
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_101_ATLAS_13  x86_64-centos7-gcc11-opt
```

## Event Generation using MG5 amc@NLO
- Workflow to make events from MG5 --> Shower in Pythia --> Do Delphes Sim
  - Prerequisits
    - Make the showering executable. All the pythia related stuff has to be configured here 
        ```
        # Make the showering code's exe
        make showerLHE
        # Usage
        ./showerLHE.exe configs/lhe.cfg <PYTHIA_CFG> <INPUT_LHE> <OUTPUT_HEPMC3>
        ```
    - Setup and validate a Madgraph intallation
        
  - Define proc card ( sample proc card : `data/proc_cards/vHH.dat` )
  - Define Forced decay file , to be used during the showering step. All the pythia related configurations can go here ( sample proc card : `data/vvHTo2G_lheDecay.configs` )
  - Have a valid delphis card ( sample cms card here : `data/delphes_card_CMS.tcl` )
  - Generate Condor jobs for event generation
  ```
  python3 makeMagraphEvents.py  -h
  python3 makeMagraphEvents.py  -v v0 -n 100 -e 4000 -s
  ```

