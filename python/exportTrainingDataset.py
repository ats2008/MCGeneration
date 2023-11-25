import h5py
import awkward as ak
import numpy as np
import vector
import os
import json,argparse
vector.register_awkward()

def writeHDF5DatasetFromDict(h5File,input_dict,basePath=''):
    for ky in input_dict:
#         print("Processing ",ky)
        if isinstance(input_dict[ky],dict):
            writeHDF5DatasetFromDict(h5File,input_dict[ky],basePath+ky+'/')
        else:
            print("Writing out ",basePath+ky)
            h5File[basePath+ky]=input_dict[ky]

JET_NMAX=16

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',"--fileIn", help="Input file ")
    parser.add_argument('-p',"--parSet", help="Parent collection to store ",default='all')
    parser.add_argument(     "--isSignal", help="Set the signal flag 1",action='store_true')
    parser.add_argument('-o',"--oFname", help="Parent collection to store ",default='workarea/ofile.h5')
    args=parser.parse_args()
    
    filepath=args.fileIn
    data= ak.from_parquet(filepath)
    parentKeySetALL = ['higgs','topQ','wNonTop','wBoson','zBoson']
    parentKeySet = parentKeySetALL
    if args.parSet=='vsZW':
        parentKeySet = ['higgs','wBoson','zBoson']
    elif args.parSet=='vsZWExclusive':
        parentKeySet = ['higgs','wNonTop','zBoson']
    elif args.parSet=='vsTTBar':
        parentKeySet = ['higgs','topQ']
    elif args.parSet=='all':
        pass
    else:
        print("parset possibinities : all , vsZQ , vsZWExclusive , vsTTBar " )
        exit()
    Jet = data.Jet
    alldaughters=[]
    
    dummyParticle=ak.zip({'pt':0.0,'eta':1e3,'phi':1e3,'mass':0.0},with_name='Momentum4D')
    dummyP1 = dummyParticle ; dummyP1['dau1'] = dummyParticle ; dummyP1['dau2'] = dummyParticle
    dummyP2 = dummyParticle ; dummyP2['dau1'] = dummyParticle ; dummyP2['dau2'] = dummyParticle ; dummyP2['dau3'] = dummyParticle
    dummyParents={
        'higgs' : dummyP1,
        'topQ'  : dummyP2,
        'wNonTop' : dummyP1,
        'wBoson' : dummyP1,
        'zBoson' : dummyP1,
    }


    nParentToFill={
        'higgs' : 2,
        'topQ'  : 2,
        'wNonTop' : 2,
        'wBoson' : 3,
        'zBoson' : 2
    }

    exportMap={}
    qIdx=0
    print(f"Exporting dataset for the parent key {parentKeySet=}")
    numGenParents={}
    for parentKey in parentKeySetALL:
        parent = data[parentKey][ np.abs(data[parentKey]['dau1'].PID) < 6 ]
        numGenParents[parentKey] = np.asarray(ak.num(parent),dtype=int)
        if parentKey not in parentKeySet:
            continue
        print(f"{parentKey} Multiplicity :  {np.unique(ak.num(parent),return_counts=True)}")
        parent=ak.fill_none(ak.pad_none(parent,nParentToFill[parentKey]), dummyParents[parentKey])
        for i in range(nParentToFill[parentKey]):
            exportMap[f"{parentKey}{i}"]={}
            exportMap[f"{parentKey}{i}"]["q1"] = qIdx ;  qIdx+=1
            exportMap[f"{parentKey}{i}"]["q2"] = qIdx ;  qIdx+=1
            if parentKey=='topQ':
                exportMap[f"{parentKey}{i}"]["q3"] = qIdx ;  qIdx+=1
        print(f"{parentKey} Multiplicity padding Fill :  {np.unique(ak.num(parent),return_counts=True)}")

        for fld in parent.fields:
            if fld in ['dau1','dau2','dau3']:
                alldaughters.append( parent[fld] )
                print("  Adding ",parentKey,fld)
    daughters = ak.concatenate( alldaughters ,axis=1)
    print("Number of Gen Quarks to match , multiplicity : " ,np.unique(ak.num(daughters) , return_counts=True))
    left,right=ak.unzip(ak.argcartesian([daughters ,Jet]))
    pairsToMatch = ak.zip({'daughters' : daughters[left], 'jet' : Jet[right] } )
    
    deltaR=pairsToMatch.daughters.deltaR(pairsToMatch.jet )
    srtIdx = ak.argsort(deltaR)
    srtIdx = srtIdx[ deltaR[srtIdx] < 0.4 ] 
    # for i in range(3):
    #     print("I = ",i)
    #     print("\t l : ",left[i])
    #     print("\t r : ",right[i])
    #     print("\t srtI : ",srtIdx[i])
    #     print("\t dr : ",deltaR[i])
    #     print("\t dr_srtd : ",deltaR[srtIdx][i])
    
    nMAX= qIdx  #np.max(ak.num(left[srtIdx]))
    nMAX= np.max(ak.num(left[srtIdx]))
    leftSorted=ak.fill_none(ak.pad_none(left[srtIdx],nMAX+1),nMAX)
    parent_dau_idx=ak.to_numpy(leftSorted)
    
    rightSorted=ak.fill_none(ak.pad_none(right[srtIdx],nMAX+1),-1)
    jets_matches=ak.to_numpy(rightSorted)
    
    CONST_SIZE =  qIdx + 1 #np.max(parent_dau_idx) + 1 
    parent_daus=np.ones(( len(parent_dau_idx) , CONST_SIZE ),dtype=int)*-1
    # Making the genQ--> JETs Maping unique . 
    # Improvement needed here by selecting the one with dR smaller rather than the one that comes first
    tMap=np.zeros(nMAX+1,dtype=bool)
    for i in range(len(jets_matches)):
        uq,uid= np.unique(jets_matches[i],return_index=True)
        tMap     =tMap*0
        tMap[uid]=1
        tMap=np.logical_not(tMap)
        jets_matches[i][tMap]=-1
    
    tmp = np.arange(len(parent_dau_idx)).reshape((len(parent_dau_idx),1))
    tmp = np.broadcast_to(np.arange(len(tmp)).reshape((len(tmp),1)),parent_dau_idx.shape)
    evtIdx= np.ndarray.flatten( tmp )
    
    higsDauIdx=np.ndarray.flatten(parent_dau_idx)
    jetsIdx   =np.ndarray.flatten(jets_matches)
    parent_daus[evtIdx,higsDauIdx]=jetsIdx
    
    
    ###   DATA EXPORT
    
    pt=ak.to_numpy(ak.fill_none(ak.pad_none(data.Jet.pt,JET_NMAX),1.0))
    mask= pt  > 0.0 # JET MASK
    
    eta=ak.to_numpy(ak.fill_none(ak.pad_none(data.Jet.eta,JET_NMAX),-8.0))
    phi = ak.to_numpy(ak.fill_none(ak.pad_none(data.Jet.phi,JET_NMAX),0.0))
    s_phi = np.sin(phi)
    c_phi = np.cos(phi)
    
    mass=ak.to_numpy(ak.fill_none(ak.pad_none(data.Jet.mass,JET_NMAX),1.0))
    btag=ak.to_numpy(ak.fill_none(ak.pad_none(data.Jet.BTag,JET_NMAX),-1.0))
    
    scalar_ht=np.ndarray.flatten(ak.to_numpy(data.ScalarHT['pt']))
    
    met=np.ndarray.flatten(ak.to_numpy(data.MissingET.pt))
    _phi=np.ndarray.flatten(ak.to_numpy(data.MissingET.phi))
    met_s_phi= np.sin(_phi)
    met_c_phi= np.cos(_phi)
    
    diphoton_pt =  ak.to_numpy(data.diphotons.pt)
    diphoton_eta =  ak.to_numpy(data.diphotons.eta)
    diphoton_c_phi =  np.cos(ak.to_numpy(data.diphotons.phi))
    diphoton_s_phi =  np.sin(ak.to_numpy(data.diphotons.eta))
    diphoton_mass =  ak.to_numpy(data.diphotons.mass)
    
    data_structure={}
    data_structure['INPUTS']={}
    
    data_structure['INPUTS']['jets']={}
    data_structure['INPUTS']['jets']['MASK']  = mask
    data_structure['INPUTS']['jets']['pt']  = pt
    data_structure['INPUTS']['jets']['eta'] = eta
    data_structure['INPUTS']['jets']['s_phi'] = s_phi
    data_structure['INPUTS']['jets']['c_phi'] = c_phi
    data_structure['INPUTS']['jets']['mass'] = mass
    data_structure['INPUTS']['jets']['btag'] = btag
    
    data_structure['INPUTS']['global']={}
    data_structure['INPUTS']['global']['ht'] = scalar_ht
    data_structure['INPUTS']['global']['met'] = met
    data_structure['INPUTS']['global']['met_s_phi'] = met_s_phi
    data_structure['INPUTS']['global']['met_c_phi'] = met_c_phi
    data_structure['INPUTS']['global']['diphoton_pt'] = diphoton_pt
    data_structure['INPUTS']['global']['diphoton_eta'] = diphoton_eta
    data_structure['INPUTS']['global']['diphoton_s_phi'] = diphoton_s_phi
    data_structure['INPUTS']['global']['diphoton_c_phi'] = diphoton_c_phi
    data_structure['INPUTS']['global']['diphoton_mass'] = diphoton_mass
    data_structure['TARGETS']={}
    # data_structure['TARGETS']["h1"]={}
    # data_structure['TARGETS']["h1"]["b1"]=parent_daus[:,0]
    # data_structure['TARGETS']["h1"]["b2"]=parent_daus[:,1]
    # data_structure['TARGETS']["h2"]={}
    # data_structure['TARGETS']["h2"]["b1"]=parent_daus[:,2]
    # data_structure['TARGETS']["h2"]["b2"]=parent_daus[:,3]
    for ky in exportMap:
        data_structure['TARGETS'][ky]={}
        for dky in exportMap[ky]:
            partIdx = exportMap[ky][dky]
            data_structure['TARGETS'][ky][dky] = parent_daus[:,partIdx]
    
    data_structure['CLASSIFICATIONS']={}
    data_structure['CLASSIFICATIONS']['EVENT']={}
    isSignal=0
    if args.isSignal:
        isSignal=1
    data_structure['CLASSIFICATIONS']['EVENT']['signal']=np.ones(len(data),dtype=np.int64)*isSignal
    
    ## Adding the Miscelaneous branches
    data_structure['MISC'] = {}
    for ky in numGenParents:
        print(" ky " ,ky ,np.unique(numGenParents[ky],return_counts=True))
        data_structure['MISC'][ky]=numGenParents[ky]
    
    ## SETTING PROPER DATATYPES
    for head in data_structure['INPUTS']:
        for ky in data_structure['INPUTS'][head]:
            if ky=='MASK':
                continue
            data_structure['INPUTS'][head][ky]=data_structure['INPUTS'][head][ky].astype(np.float32)
            
    for head in data_structure['TARGETS']:
        for ky in data_structure['TARGETS'][head]:
            data_structure['TARGETS'][head][ky]=data_structure['TARGETS'][head][ky].astype(np.int64)
            
    for head in data_structure['CLASSIFICATIONS']:
        for ky in data_structure['CLASSIFICATIONS'][head]:
            data_structure['CLASSIFICATIONS'][head][ky]=data_structure['CLASSIFICATIONS'][head][ky].astype(np.int64)
    
    print("Exporting file to ",args.oFname)
    oFile = h5py.File(args.oFname,'w')
    writeHDF5DatasetFromDict(oFile,data_structure)
    oFile.close()        



if __name__=='__main__':
    main()

