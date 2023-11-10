#!/usr/bin/env python 
from __future__ import print_function
import os,glob
import copy,json
import argparse

condorScriptString="\
executable = $(filename)\n\
output = $Fp(filename)run.$(Cluster).stdout\n\
error = $Fp(filename)run.$(Cluster).stderr\n\
log = $Fp(filename)run.$(Cluster).log\n\
+JobFlavour = \"longlunch\"\n\
"
pwd=os.environ['PWD']
proxy_path=os.environ['X509_USER_PROXY']
HOME=os.environ['HOME']
NJOBS=5
NEVENTS_PER_JOB = 1000
ZERO_OFFSET=1
FILES_PER_JOB=1
maxMeterialize=100
offsetStep=0
RESULT_BASE='/grid_mnt/t3storage3/athachay/trippleHiggs/delphisStudies/MCGeneration/results'
JOB_TYPE='mg5Generation'

parser = argparse.ArgumentParser()
parser.add_argument('-s',"--submit", help="Submit file to condor pool", action='store_true' )
parser.add_argument('-r',"--resubmit", help="Re-Submit file to condor pool", action='store_true' )
parser.add_argument('-t',"--test", help="Test Job", action='store_true' )
parser.add_argument('-n',"--njobs", help="Number of jobs to make",default=-6000,type=int)
parser.add_argument('-e',"--nevts", help="Number of events per job",default=-1,type=int)
parser.add_argument('-c',"--config", help="Configuration file")
parser.add_argument(     "--cats", help="JobTags to process",default='')
parser.add_argument(     "--offset", help="offset to the jobs / random seeds  ",default=1, type=int)
parser.add_argument('-v',"--version", help="Vesion of the specific work",default='TEST')

args = parser.parse_args()

version=args.version

njobs=int(args.njobs)
maxevents=int(args.nevts)
max_meterialize=250
jobsToProcess=None
submit2Condor=args.submit
resubmit2Condor=args.resubmit
isTest=args.test
job_hash=f'mg5Gen_{args.version}'
ZERO_OFFSET=args.offset

print(" submit jobs ",submit2Condor)
print(" resubmit jobs ",resubmit2Condor)
print(" isTest ",isTest)
print(" njobs ",njobs)
print(" maxEvt ",maxevents)
print(" rand seed offset  ",ZERO_OFFSET)

if submit2Condor or resubmit2Condor:
    choice=input("Do you really want to submit the jobs to condor pool ? ")
    if 'y' not in choice.lower():
        print("Exiting ! ")
        exit(0)

jobDict={}
jobConfigJson= args.config
with open(jobConfigJson) as f:
    print("Loading Job Config ",jobConfigJson)
    jobDict=json.load(f)
parameters={}    
if "PARAMETERS" in jobDict:
    parameters.update(jobDict.pop("PARAMETERS"))
executable=parameters["MG5_EXE"]
if isTest :
    args.njobs=2
    args.nevts=1000

allCondorSubFiles=[]
if jobsToProcess==None:
    jobsToProcess=list( jobDict.keys() )
if args.cats:
    jobsToProcess=args.cats.split(",")

print()
print(" Processing Job categories : ",jobsToProcess)
print()
i=0

for jobTag in jobsToProcess:
    jobTag=jobTag
    if jobTag not in jobsToProcess:
        print(f"Skipping {jobTag}")

    runScriptTemplate = jobDict[jobTag]['script'] 
    runScriptTxt=[]
    with open(runScriptTemplate,'r') as f:
        runScriptTxt=f.readlines()
    runScriptTxt=''.join(runScriptTxt)
    
    proc_card_pth=jobDict[jobTag]['proc_card']
    proc_card=[]
    with open(proc_card_pth,'r') as f:
        proc_card=f.readlines()
    proc_tag="tproc"
    for l in proc_card:
        if 'output' in l:
            proc_tag=l.strip().split(" ")[-1].strip()
    print(f"proc tag set to {proc_tag}")
            
    proc_card=''.join(proc_card)
    
    launch_card_pth=jobDict[jobTag]['launch_card']
    launch_card=[]
    with open(launch_card_pth,'r') as f:
        launch_card=f.readlines()
    launch_card=''.join(launch_card)
    
    #if proc not in procToProcess:
    #    continue
    head=pwd+f'/Condor/{JOB_TYPE}/{job_hash}/{jobTag}/'
    njobs_toMake=NJOBS 
    if 'njobs' in jobDict[jobTag]:  njobs_toMake=jobDict[jobTag]['njobs']
    if args.njobs > 0: njobs_toMake = args.njobs
    maxevents=NEVENTS_PER_JOB
    if 'maxevents' in jobDict[jobTag]:maxevents=jobDict[jobTag]['maxevents']
    if args.nevts > 0 : maxevents=args.nevts
        
    print(f"Making {njobs_toMake} Jobs with {maxevents} each , for {job_hash}/{jobTag}")
    njobs_made=0
    for ii in range(njobs_toMake):
        if(ii%10==0) : print("\nJob Made : ",end = " " )
        print(ii,end =" ")
        i = ii +ZERO_OFFSET

        njobs_made+=1
        dirName =f'{head}/Job_{i}/'
        if not os.path.exists(dirName):
            os.system('mkdir -p '+dirName)
        destination=f'{RESULT_BASE}/{JOB_TYPE}/{job_hash}/{jobTag}/out_{i}/'
        if not os.path.exists(destination):
            os.system('mkdir -p '+destination)
        
        ## Setting PROC CARD 
        procFileName='proc_'+str(i)+'.dat'
        with open(dirName+'/'+procFileName,'w') as f:
            f.write(proc_card)
        
        ## Setting LAUNCH CARD 
        tmp=launch_card.replace("@@IDX",f"{i}")
        tmp=tmp.replace("@@TAG",proc_tag)
        launchFileName='launch_'+str(i)+'.dat'
        with open(dirName+'/'+launchFileName,'w') as f:
            f.write(tmp)

        ## RUNNING WRAPPER SCRIPT
        runScriptName=dirName+f'/{job_hash}_{proc_tag}_{i}_run.sh'
        if os.path.exists(runScriptName+'.sucess'):
           os.system('rm '+runScriptName+'.sucess')
        runScript=open(runScriptName,'w')
        tmp=runScriptTxt.replace("@@DIRNAME",dirName)
        tmp=tmp.replace("@@IDX",str(i))
        tmp=tmp.replace("@@PTAG",proc_tag)
        tmp=tmp.replace("@@MAXEVENTS"  ,str(maxevents))
        tmp=tmp.replace("@@PROCFILENAME"  ,procFileName)
        tmp=tmp.replace("@@LAUNCHFILENAME",launchFileName)
        tmp=tmp.replace("@@RUNSCRIPT",runScriptName)
        tmp=tmp.replace("@@MG5_EXECUTABLE",executable)
        tmp=tmp.replace("@@proxy_path",proxy_path)
        tmp=tmp.replace("@@HOME",HOME)
        tmp=tmp.replace("@@DESTINATION",destination)
        if 'decayFile' in jobDict[jobTag]:
            tmp=tmp.replace("@@DECAY_FILE",jobDict[jobTag]['decayFile'])
        if 'delphesCard' in jobDict[jobTag]:
            tmp=tmp.replace("@@DEPHISCARD",jobDict[jobTag]['delphesCard'])
            
        runScript.write(tmp)
        runScript.close()
        os.system('chmod +x '+runScriptName)

    print()
    condorScriptName=head+f'/job_{job_hash}_{jobTag}.sub'
    with open(condorScriptName,'w') as condorScript:
        condorScript.write(condorScriptString)
        condorScript.write("queue filename matching ("+head+"/*/*.sh)\n")
        rpth=os.path.relpath(condorScriptName)
        rFolder ='/'.join(rpth.split("/")[:-1])
        print(f"Condor {njobs_made} Jobs made !\n\t Job Directory : {rFolder}/ \n\t submit file  : {rpth}")
    allCondorSubFiles.append(condorScriptName)

    if resubmit2Condor:
        continue
print("")
print("")
print("All condor submit files to be submitted ")
for fle in allCondorSubFiles:
    print('condor_submit '+fle)
    if submit2Condor or resubmit2Condor:
        os.system('condor_submit '+fle)
print("")
print("")
