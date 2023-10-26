#!/bin/bash
source /grid_mnt/t3home/athachay/.bashrc
source /cvmfs/cms.cern.ch/cmsset_default.sh 
set -x
export HOME=@@HOME
export X509_USER_PROXY=@@proxy_path
cd @@DIRNAME
set +x
set -x
TMPDIR=`mktemp -d`
cd $TMPDIR
cp @@DIRNAME/@@PROCFILENAME .
cp @@DIRNAME/@@LAUNCHFILENAME .
mv @@RUNSCRIPT @@RUNSCRIPT.busy 

@@MG5_EXECUTABLE @@PROCFILENAME

cat @@PTAG/Cards/run_card.dat | grep -E "iseed|nevents"
sed -i 's/.*=.*nevents.*!/ @@MAXEVENTS = nevents !/g' @@PTAG/Cards/run_card.dat
sed -i 's/0.*=.*iseed.*!/ @@IDX = iseed       !/g' @@PTAG/Cards/run_card.dat
cat @@PTAG/Cards/run_card.dat | grep -E "iseed|nevents"
@@MG5_EXECUTABLE @@LAUNCHFILENAME

set +x
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_104  x86_64-centos7-gcc11-opt
set -x 

cp @@PTAG/Events/run_@@IDX/events.lhe.gz .
gunzip events.lhe.gz 
/grid_mnt/t3storage3/athachay/trippleHiggs/delphisStudies/MCGeneration/showerLHE.exe /grid_mnt/t3storage3/athachay/trippleHiggs/delphisStudies/MCGeneration/configs/lhe.cfg @@DECAY_FILE events.lhe events.hepmc3
/grid_mnt/t3storage3/athachay/trippleHiggs/delphisStudies/Delphes-3.5.0/DelphesHepMC3  @@DEPHISCARD events.root events.hepmc3
if [ $? -eq 0 ]; then 
#    cp @@PTAG/Events/run_@@IDX/* @@DESTINATION
    cp *.lhe *.root @@PTAG/Events/run_@@IDX/run*.txt @@DESTINATION
    if [ $? -eq 0 ] ; then
        mv @@RUNSCRIPT.busy @@RUNSCRIPT.sucess 
        echo OK
    else
        mv @@RUNSCRIPT.busy @@RUNSCRIPT
        echo FAIL AT COPY
    fi
else
    mv @@RUNSCRIPT.busy @@RUNSCRIPT 
    echo FAIL AT EXE
fi
