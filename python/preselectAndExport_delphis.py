import awkward as ak
import numpy as np
import matplotlib.colors as cols
import vector
import matplotlib.pyplot as plt
import mplhep as hep
import json,glob,argparse,os,time
import ROOT;
vector.register_awkward()

import uproot as urt


allKeys=['Jet','Photon','Electron','Photon','Muon','MissingET','ScalarHT','GenJet','Particle']
keysToStore=["PID","Status","pt","eta","phi","mass"]
blackListed=['fUniqueID' , 'fBits' , '[',']','Constituents','Particles']
particlesToStore={ 
                   'higgs'   : 25,
                   'wNonTop'  : 24,
                   'zBoson'     : 23,
                   'wBoson'  : 24,
                   'topQ'    : 6
                }

def updtaeteTopDaughters(parents,daughter1,daughter2,data):
    print(parents)
    dauArrys=[]
    for EID in range(len(parents)):
        dauArr=[]
        if EID%100==0:
            print(f"\r Updating Top {EID} / {len(parents) }   ",end="")
        topF=0
        for TID in range(len(parents[EID])):
            d1 = daughter1[EID][TID].PID#data.Particle[EID][ parents[EID][TID].D1 ].PID
            d2 = daughter2[EID][TID].PID#data.Particle[EID][ parents[EID][TID].D2 ].PID
            dau3=None
            if abs(d1) == 24:
                mother = daughter1[EID][TID]    #data.Particle[EID][ parents[EID][TID].D1 ]
                dau3   = daughter2[EID][TID]    #data.Particle[EID][ parents[EID][TID].D2 ] 
            elif abs(d2) == 24:
                mother = daughter2[EID][TID]    #data.Particle[EID][ parents[EID][TID].D2 ]
                dau3   = daughter1[EID][TID]    #data.Particle[EID][ parents[EID][TID].D1 ] 
            else:
                continue
            topF+=1
            sucess=False
            for  i in range(20):
                d1=mother.D1
                d2=mother.D2
                dau1=data.Particle[EID][d1]
                dau2=data.Particle[EID][d2]
                if (dau1.PID == mother.PID) :   # RADIATION SKIP
                    mother=data.Particle[EID][d1]
                elif (dau2.PID == mother.PID) :   # RADIATION SKIP
                    mother=data.Particle[EID][d2]
                else:
                    sucess=True
                    break
            if not sucess:
                print("W DECAY NOT FOUND")
            else:
                dauArr.append(ak.zip({
                                      "dau1" : data.Particle[EID][d1] ,
                                      "dau2" : data.Particle[EID][d2] ,
                                      "dau3" : dau3
                                     } ))
        if topF<2:
            print("Spurious Event!! " )
            for TID in range(len(parents[EID])):
                d1 = daughter1[EID][TID]
                d2 = daughter2[EID][TID]
                print(d1,d2)
        dauArrys.append(ak.Array(dauArr))
    print()
    return ak.Array(dauArrys)
#     break


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v',"--version", help="Base trigger set",default='')
    parser.add_argument(     "--cfg"    , help="Input configurateion")
    parser.add_argument(     "--flist"    , help="Force an input fileList",default=None)
    parser.add_argument(     "--procs"    , help="Procs to process", default=None)
    
    args = parser.parse_args()
    
    config={
        'fileListIn'  : None   ,
        'destination' : './' ,
        'saveOutPut' : False
    }
    
    with open(args.cfg) as f:
        config.update(json.load(f))
    if args.flist:
        config["fileList"] = args.flist
        print("Input file list set to ",config["fileListIn"])
    
    fileListIn= config["fileList"] #'../workarea/delphes_CMS/v0/flist.json'
    destination=config['destination']+'/'+args.version+'/'
    print(f"Destination Set to {destination}")
    fileDict={}
    print(f"Opening file List from {fileListIn}")
    with open(fileListIn) as f:
        fileDict=json.load(f)

    procsToProcess=list(fileDict.keys())
    if args.procs:
        procsToProcess=args.procs.split(",")
    for proc in fileDict:
        if proc not in procsToProcess:
            continue
        prefix=f"{destination}/{proc}/"
        if not os.path.exists(prefix):
            print(f"** maing directory : {prefix}")
            os.system(f'mkdir -p {prefix}')
        nFileMax=len(fileDict[proc])
        print(f"Processing {proc} [ {nFileMax} Files ]")
        fIdx=0
        for fileName in fileDict[proc]:
            fIdx+=1
            print(f"\t Processing file [ {fIdx} / {nFileMax} ] : {fileName}")
            with urt.open(fileName)  as f:
                dataStoreRAW={}
                for object in allKeys:
                    print(f"Adding the object {object}")
                    tree=f['Delphes'][object]
                    branches=[]
                    for ky in tree.keys():
                        skip=False
                        for chk in blackListed:
                            if chk in ky:
                                skip=True
                                break
                        if skip:
                            continue
                        branches.append(ky)
                    data=tree.arrays(branches)
                    
                    base_key=object
                    br_list=[]
                    for ky in data.fields:
                        br=ky.split('.')[-1]
                        if br=='PT' : br='pt'
                        if br=='MET' : br='pt'
                        if br=='HT' : br='pt'
                        if br=='Eta': br='eta'
                        if br=='Mass': br='mass'
                        if br=='Phi': br='phi'
                        data[br]=data[ky]
                        br_list.append(br)
                    if 'mass' not in data.fields:
                        data['mass']=ak.zeros_like(data.pt)
                        br_list.append('mass')
                    data=data[br_list]
                    data= ak.zip({ kyz : ak.flatten(data[kyz],axis=0) for kyz in data.fields} ,with_name='Momentum4D')
                    dataStoreRAW[object]=data
                    # break
    
                
                data = ak.Array(dataStoreRAW)
                print(f"\t Loaded {len(data)} Events ")
                diphotons = ak.combinations(
                    data.Photon,
                    2,
                    fields = ["LeadPhoton", "SubleadPhoton"]
                )
                
                diphotons['p4'] = diphotons.LeadPhoton +diphotons.SubleadPhoton
                sortedIdx=ak.argsort(diphotons.p4.pt,ascending=False)[:,:1]
                diphotons=diphotons[sortedIdx]
                
                for ky in ['pt','eta','phi','mass']:
                    diphotons[ky]=getattr(diphotons["p4"],ky)
                
                data['diphotons']=diphotons[['pt','eta','phi','mass',"LeadPhoton", "SubleadPhoton"]]
    
                ##  Selections
                
                # Lepton Veto
                
                electrons_sel = ( data.Electron.IsolationVar  < 0.1 ) & ( data.Electron.pt  > 25.0 )
                electrons = data.Electron[electrons_sel]
                mask_electron = ak.num(electrons) < 1
                print(f"Electron veto : P -> {np.sum(mask_electron)} , T -> {len(mask_electron)} , eff : {np.sum(mask_electron)/len(mask_electron):.3f}")
                data  = data[ mask_electron ]

                muons_sel  = ( data.Muon.IsolationVar  < 0.1 ) & ( data.Muon.pt  > 25.0 )
                muons = data.Muon[muons_sel]
                mask_muons = ak.num(muons) < 1
                print(f"Muon veto : P -> {np.sum(mask_muons)} , T -> {len(mask_muons)} , eff : {np.sum(mask_muons)/len(mask_muons):.3f}")
                data  = data[ mask_muons ]

                # Diphotons
                NEVTS = len(data)
                mgg_selection = (data.diphotons.mass > 80.0) &  (data.diphotons.mass < 180.0)
                data["diphotons" ]= data.diphotons[mgg_selection]
                mask_mgg =  ak.num( data.diphotons ) > 0
                print(f"Mass Selection : P -> {np.sum(mask_mgg)} , T -> {len(mask_mgg)} , eff : {np.sum(mask_mgg)/len(mask_mgg):.3f}")
                data = data[mask_mgg]
                data["diphotons" ]  = ak.flatten(data.diphotons)

                # Jets
                Jets = data.Jet
                drsL  =Jets.deltaR(data.diphotons.LeadPhoton)
                drsSL =Jets.deltaR(data.diphotons.SubleadPhoton)
                
                jetMask = Jets.pt > 20
                jetMask = jetMask & ( np.abs(Jets.eta) < 2.5 )
                jetMask = jetMask & ( drsL > 0.4 ) & ( drsSL > 0.4 )
                data["Jet"] = data.Jet[jetMask]
                
                mask_jets = ak.num(Jets) > 3
                print(f"Jet Selection : P -> {np.sum(mask_jets)} , T -> {len(mask_jets)} , eff : {np.sum(mask_jets)/len(mask_jets):.3f}")
                data = data[mask_jets]
    
                # Addding the GEN Info

                genDictStore={}
                print("loading daughter info ! ..")
                allParticlesD1 = ak.to_numpy(ak.fill_none(ak.pad_none(data.Particle.D1,np.max(ak.num(data.Particle.D1))),-1))
                allParticlesD2 = ak.to_numpy(ak.fill_none(ak.pad_none(data.Particle.D2,np.max(ak.num(data.Particle.D2))),-1))
                for particleName in particlesToStore:
                    genDict={}
                    pid = particlesToStore[particleName]
                    print(f"Making the {particleName} entries with PID {pid}")
                    
                    mask      = np.abs(data.Particle.PID) == pid
                    parents   = data.Particle[ mask ]
                    grandParent = data.Particle[parents.M1]
                    cleanMask = np.logical_not((grandParent.D1 == grandParent.D2) &  (grandParent.D1 == grandParent.D2) )
                    parents   = parents[ cleanMask ] 
                    daughter1 = data.Particle[ parents.D1 ]
                    daughter2 = data.Particle[ parents.D2 ]
                    
                    if True or (particleName=='wNonTop'):
                        ## Clean the decay chain
                        d1=ak.to_numpy(ak.fill_none(ak.pad_none(parents.D1,np.max(ak.num(parents.D1))),-1))
                        d2=ak.to_numpy(ak.fill_none(ak.pad_none(parents.D2,np.max(ak.num(parents.D2))),-1))
                        isDummy = np.logical_and( d1==d2 , d1 >=1)
                        resolved=np.sum(isDummy)
                        lIdx=0
                        while (resolved > 0) and (lIdx < 120 ):
                            eid,tid=np.where(isDummy)
                            d1Id = d1[eid,tid]
                            d2Id = d2[eid,tid]
                            d1[eid,tid] = allParticlesD1[eid,d1Id]
                            d2[eid,tid] = allParticlesD2[eid,d2Id]
                            isDummy = np.logical_and( d1==d2 , d1 >=1)
                            lIdx+=1
                            resolved=np.sum(isDummy)
                            print("\r --> ",lIdx," to be resolved : ",resolved,end="")
                        print()
                        if resolved >0 :
                            print("\t\tDECAY CHAINS NOT YET RESOLVED ! to be resolved : ",resolved)
                        nCounts=np.sum( d1 >= 0 ,axis=1)
                        d1=ak.unflatten(d1[ d1 >=0 ],nCounts)
                        d2=ak.unflatten(d2[ d2 >=0 ],nCounts)
                        
                        daughter1 = data.Particle[ d1 ]
                        daughter2 = data.Particle[ d2 ]
                    
                    grandParent = data.Particle[parents.M1]
                    cleanMask = np.logical_not((grandParent.D1 == grandParent.D2) &  (grandParent.D1 == grandParent.D2) )
                    if particleName=='topQ':
                        topDecaySelector = ( np.abs(daughter1.PID) < 6 ) | (np.abs(daughter2.PID) < 6)
                        cleanMask = cleanMask & topDecaySelector
                    if particleName=="wNonTop":
                        cleanMask = cleanMask & ( np.abs(data.Particle[ parents.M1 ].PID ) != 6 ) 
                        cleanMask = cleanMask & ( np.abs(data.Particle[ parents.M1 ].PID ) != 24 ) # cleaning the spurious w to w
                    if particleName=='higgs':
                        pass
                    if particleName=="wBoson":
                        cleanMask = cleanMask & ( np.abs(data.Particle[ parents.M1 ].PID ) != 24 ) # cleaning the spurious w to w
                        pass
                    
                    parents = parents[cleanMask]
                    daughters = {
                        "dau1" : daughter1[cleanMask],
                        "dau2" : daughter2[cleanMask], 
                    }
                    if particleName=='topQ':
                        daughters["dau3"] = daughters["dau2"] 
                        if np.sum(ak.num( parents ) ) > 0:
                            dauArrays=updtaeteTopDaughters(parents,daughters["dau1"],daughters["dau2"],data)
                            daughters = {
                                "dau1" : dauArrays.dau1,
                                "dau2" : dauArrays.dau2,
                                "dau3" : dauArrays.dau3
                            }
                #          for evtIn range(len(parents)):
                
                    genDict[particleName] = ak.zip({
                        ky : parents[ky] for ky in keysToStore
                    }, with_name = 'Momentum4D')
                    
                    for dky in daughters:
                        genDict[particleName][dky]=ak.zip({
                            ky : daughters[dky][ky] for ky in keysToStore
                        }, with_name = 'Momentum4D')
                            
                    genDictStore[particleName] = genDict[particleName]
                for ky in genDictStore:
                    data[ky]=genDictStore[ky]
                    print(f"Multiplicity of {ky} {np.unique(ak.num(data[ky]),return_counts=True)}")
                    #print(f"Average number of {ky} {np.average(ak.num(data[ky])):.3f}")
                
                # EXPORT 
                if config['saveOutPut']:
                    NEVTS_SELECTED= len(data)
                    oFname=prefix+ fileName.split("/")[-1].replace('.root',f'_{fIdx}_preselected.parquet')
                    print(f"file being saved to : {oFname} ")
                    ak.to_parquet(data,oFname)
                    foutName=oFname.replace(".parquet",".txt")
                    with open(foutName,'w') as f:
                        f.write(f"Input file : {fileName}\n")
                        f.write(f"Total Number of Events in Input : {NEVTS}\n")
                        f.write(f"Total Number of Events after Selections : {NEVTS_SELECTED}\n")
                        f.write(f"Efficiency : {1.0*NEVTS_SELECTED/NEVTS:.3f}")
                    print("\t\t DONE !")
            break    

if __name__=='__main__':
    main()
