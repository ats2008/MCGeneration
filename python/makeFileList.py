import glob,argparse,os,json 

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--base"    , help="Input directory")
parser.add_argument('-o', "--out"    , help="Output fname",default='out.json')
parser.add_argument('-e', "--ext"    , help="WildCard Extensions for files",default="*.root")
args = parser.parse_args()

base=args.base
print("Processing base  = ",base)
folderList=[]
for fl in glob.glob(base+"/*"):
    if os.path.isdir(fl):
        folderList.append(fl)
print(f"Obtained {len(folderList)} folderlist !")

outDict={}
for folder in folderList:
    searchStr=folder+"/"+args.ext
    tag=folder.split("/")[-1]
    print(f"Processing {tag} with the search str : {searchStr}")
    fList=glob.glob(searchStr)
    outDict[tag]=[]
    for fl in fList:
        outDict[tag].append(os.path.abspath(fl))
    print(f"\t Got {len(outDict[tag])} files for {tag}")

with open(args.out,'w') as f:
    print("Writing out ",args.out)
    json.dump(outDict,f,indent=4)

