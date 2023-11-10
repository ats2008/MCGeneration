import awkward as ak
import glob,argparse 

parser = argparse.ArgumentParser()
parser.add_argument('-v',"--version", help="Base trigger set",default='v0')
parser.add_argument("-b", "--base"    , help="Input directory")
parser.add_argument('-o', "--outDir"    , help="Output directory",default=None)
args = parser.parse_args()

baseFs=[
    args.base
]


for base in baseFs:
    base=(base+'/').replace("//","/")
    print("Processing base  = ",base)
    flist=glob.glob(base+"/*.parquet")
 
    tagName = '/'.join(base.split('/')[:-1])
    if not args.outDir:
        tagName = '/'.join(base.split('/')[:-1])
        prefix=''
        foutName=tagName+f'_{args.version}.parquet'
    else:
        tagName = base.split('/')[:-1]
        prefix=args.outDir+'/'
        foutName=prefix+'/'+tagName+f'_{args.version}.parquet'
    print(f"\t Outputs after merging {len(flist)} files will be concatinated to {foutName} !")
    oData=[]
    idx=0
    nEvts=0
    for fl in flist:
        idx+=1
        dta=ak.from_parquet( fl )
        nEvts+=len(dta)
        print(f"\r\t Adding file  {idx} / {len(flist)} : {fl}  , {len(dta)} events  [tot : {nEvts} ]   ",end="")
        oData.append( dta )
            
    
    
    outData=ak.concatenate( oData)
    print()
    print("Number of Events : ",len(outData))
    print(f"\t Output saved to : {foutName}")
    ak.to_parquet(outData,foutName)
