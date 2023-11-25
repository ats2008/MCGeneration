import argparse
import glob
import os


parser = argparse.ArgumentParser()
parser.add_argument('-i',"--inputFiles", help="Input File wildcard String ")
parser.add_argument('-p',"--prefix", help="prefix for the output file ")
parser.add_argument("--others", help=" extra flags for the scripts",default="")
parser.add_argument('-e',"--execute", help="Execute the cmd",action='store_true')
#parser.add_argument('-s',"--isSignal", help="add isSignal flag to cmd",action='store_true')
parser.add_argument('-n',"--nfiles", help="Number of files to excute for",default=-1,type=int)
args = parser.parse_args()

input_str=args.inputFiles
input_str = input_str.replace("@","*")
prefix=args.prefix
cmd_str_tpl = 'python python/exportTrainingDataset.py -f @@INFILE -o @@OFILE '
#cmd_str_tpl = 'python python/exportTrainingDataset.py -f @@INFILE -o @@OFILE -p vsZWExclusive '
#if args.isSignal:
#    cmd_str_tpl+=" --isSignal "
miscFlags=args.others.split(",")
for flg in miscFlags:
    cmd_str_tpl+=" "+flg.replace("@",'-')
print("Search string : ",input_str)
files = glob.glob(input_str)
print(f"GOT {len(files)} files")

if 'ggHHH' in cmd_str_tpl:
    if 'isSignal' not in cmd_str_tpl:
        print("isSignal falsg missing !")
        exit()

if not os.path.exists(prefix):
    os.system("mkdir -p "+prefix)
fidx=0
for ifile in files:
    if (args.nfiles > -1 ) and ( fidx >= args.nfiles):
        break
    fidx+=1
    ofile=prefix+'/'+f'spanetD_{fidx}_'+ifile.split('/')[-1].split(".")[0]+'.h5'
    cmd_str = cmd_str_tpl.replace("@@INFILE",ifile).replace("@@OFILE",ofile)
    if args.execute:
        print("== "*10,f"Doing {fidx} / {len(files)} [ n jobs max : {args.nfiles}]","== "*10)
        print(cmd_str)
        os.system(cmd_str)
        print("\n "," - * - *"*15," \n")
    else:
        print(cmd_str)
