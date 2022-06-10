#!/usr/bin/env python3
import os
import sys
import time
rSeed=int(time.time())%25000

version='v1'

"""
    Usage
    ./makeCondorJobForAnalysis.py <EXECUTABLE> <InputFileListFname> <CFG_TEMPLATE> <analysisOption> <destination> <NJOBS> <FILES_PER_JOB> <jobPrefix>"

"""

def loadConfigTemplate(fname):
    f=open(fname,'r')
    txt=f.readlines()
    f.close()
    cfgString="".join(txt)
    return cfgString    


executable='analysisNtupleMaker.exe'
analysisOption='1'


NJOBS=20000
NEVENTS_PER_JOB = -1
ZERO_OFFSET=0
FILES_PER_JOB=40
rSeed=-1

pwd=os.environ['PWD']
proxy_path=os.environ['X509_USER_PROXY']
HOME=os.environ['HOME']
xrdRedirector="root://cms-xrd-global.cern.ch/"

FileSource ="bmm5FileList.txt"
destination='/grid_mnt/t3storage3/athachay/bs2mumug/run2studies/CMSSW_10_6_19_patch2/src/BsMMGAnalysis/MergeWithBMMNtuples/RunLumiEventFileMaker/runLumiList/'
cfgTemplateFname="analysis2018Temp.cfg"
tag=""
argC=1

if len(sys.argv) > 2:
    executableSrc=sys.argv[argC]  
    executable=executableSrc.split('/')[-1]
else:
    print("Usage\n\t./makeCondorJobForAnalysis.py <EXECUTABLE> <CFG_TEMPLATE>  <destination> <NJOBS> <EVENTS_PER_JOB> <jobPrefix>")
    sys.exit(1)
argC+=1
if len(sys.argv) > argC:
    cfgTemplateFname=sys.argv[argC]  
argC+=1
if len(sys.argv) > argC:
    destination=sys.argv[argC]  
argC+=1
if len(sys.argv) > argC:
    NJOBS=int(sys.argv[argC])  
argC+=1
if len(sys.argv) > argC:
    NEVENTS_PER_JOB=int(sys.argv[argC])  
argC+=1
if len(sys.argv) > argC:
    tag=sys.argv[argC]  
argC+=1
if len(sys.argv) > argC:
    ZERO_OFFSET=int(sys.argv[argC])  
argC+=1
if len(sys.argv) > argC:
    rSeed=int(sys.argv[argC])  

if rSeed<0:
    rSeed=int(time.time())%25000
if(not os.path.exists(destination)):
    os.system("mkdir -p "+destination)
destination=os.path.abspath(destination)

head='Condor/Jobs'+tag
if not os.path.exists(head ):
    os.system('mkdir -p '+head)

print("Executable ",executable)
print("Source file list ",FileSource)
print("CFG template  file ",cfgTemplateFname)
print("analysis Type ",analysisOption)
print("destination : ",destination)
print("NJOBS : ",NJOBS)
print("FILES_PER_JOB : ",FILES_PER_JOB)
print("NEVENTS_PER_JOB : ",NEVENTS_PER_JOB)
print("tag : ",tag)
print("ZERO_OFFSET : ",ZERO_OFFSET)

configurationTxt= loadConfigTemplate(cfgTemplateFname)
configurationTxt= configurationTxt.replace("@@PWD",pwd)
configurationTxt= configurationTxt.replace("@@VERSION",version)
configurationTxt= configurationTxt.replace("@@TAG",tag)

condorScriptString="\
executable = $(filename)\n\
output = $Fp(filename)run.$(Cluster).stdout\n\
error = $Fp(filename)run.$(Cluster).stderr\n\
log = $Fp(filename)run.$(Cluster).log\n\
+JobFlavour = \"espresso\"\n\
"

condorScriptName=head+'/job'+tag+'.sub'
condorScript=open(condorScriptName,'w')
condorScript.write(condorScriptString)


runScriptTxt="\
#!/bin/bash\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh \n\
export HOME="+HOME+"\n\
export X509_USER_PROXY="+proxy_path+"\n\
cd @@DIRNAME \n\
#source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_101 x86_64-centos7-gcc11-opt \n\
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_101_ATLAS_13 x86_64-centos7-gcc11-opt\n\
set -x\n\
TMPDIR=`mktemp -d`\n\
cd $TMPDIR\n\
cp  "+pwd+"/"+executableSrc+" .\n\
cp @@DIRNAME/@@CFGFILENAME .\n\
mv @@RUNSCRIPT @@RUNSCRIPT.busy \n\
./"+executable+" @@CFGFILENAME "+analysisOption+"\n\
if [ $? -eq 0 ]; then \n\
    mv *.root "+destination+"\n\
    if [ $? -eq 0 ] ; then \n\
        mv @@CFGFILENAME " + destination + "\n\
        mv @@RUNSCRIPT.busy @@RUNSCRIPT.sucess \n\
        echo OK\n\
    else\n\
        mv @@RUNSCRIPT.busy @@RUNSCRIPT \n\
        echo FAIL\n\
    fi\n\
else\n\
    mv @@RUNSCRIPT.busy @@RUNSCRIPT \n\
    echo FAIL\n\
fi\n\
"

print("Making ",NJOBS," Jobs ")

njobs=0
for ii in range(NJOBS):
    i=ii+ZERO_OFFSET
    
    dirName= pwd+'/'+head+'/Job_'+str(i)
    
    if(ii%10==0) : print("\nJob Made : ",end = " " )
    print(ii,end =" ")

    if os.path.exists(dirName):
        k=True
    else:
        os.system('mkdir '+dirName)
    cfgFileName='mcCfg_'+str(i)+'.cfg'
    cfgFile=open(dirName+'/'+cfgFileName,'w')
    tmp=""
    tmp=configurationTxt.replace("@@PWD",pwd)
    tmp=tmp.replace("@@TAG",tag)
    tmp=tmp.replace("@@IDX",str(i))
    tmp=tmp.replace("@@RSEED",str(i+rSeed))
    tmp=tmp.replace("@@MAXEVENTS",str(NEVENTS_PER_JOB))
    cfgFile.write(tmp)
    cfgFile.close()   
    
    runScriptName=dirName+'/'+tag+'run'+str(i)+'.sh'
    if os.path.exists(runScriptName+'.sucess'):
       os.system('rm '+runScriptName+'.sucess')
    runScript=open(runScriptName,'w')
    tmp=runScriptTxt.replace("@@DIRNAME",dirName)
    tmp=tmp.replace("@@IDX",str(i))
    tmp=tmp.replace("@@CFGFILENAME",cfgFileName)
    tmp=tmp.replace("@@RUNSCRIPT",runScriptName)
    runScript.write(tmp)
    runScript.close()
    os.system('chmod +x '+runScriptName)
    condorScript.write("queue filename matching ("+runScriptName+")\n")
    njobs+=1
print()
print(" Number of jobs made : ", njobs)
print(" Condor  submit file name : ", condorScriptName)
condorScript.close()

