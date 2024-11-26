print("starting the reweighting processor ")
#standard python libraries
import re, ast
import awkward as ak
import numpy as np
import uproot
import glob as glob
import matplotlib.pyplot as plt
import argparse
import os,sys
#HEP python libraries
from coffea.nanoevents import NanoEventsFactory
from coffea.processor import AccumulatorABC, ProcessorABC, Runner, FuturesExecutor
import mplhep
#CMS libraries
from ecalautoctrl import RunCtrl
from ecalphisym import EcalPhiSymSchema
from EcalPhiSymProcessor import EcalPhiSymAccumulator, Processor




#%matplotlib inline
plt.rc('font',**{'family':'sans-serif', 'sans-serif':['Helvetica'], 'size':'20'})
plt.rcParams['figure.figsize'] = '10,10'
plt.rcParams['xaxis.labellocation'] = 'right'
plt.rcParams['yaxis.labellocation'] = 'top'
plt.style.use(mplhep.style.ROOT)



def findFilesEop(rctrl,  era, selected_filelist_era=[]):
    selected_filelist_era=[]
    ecalelf = set()

    if "2022C" in era:
        ecalelf.update(rctrl.getOutput(era = era,  process='ecalelf-ntuples'))

    elif "2022D" in era:
        ecalelf.update(rctrl.getOutput(era = era,  process='ecalelf-ntuples'))

    else:
        ecalelf.update(rctrl.getOutput(era = era, process='ecalelf-ntuples-wskim'))
        ecalelf.update(rctrl.getOutput(era = era, process='ecalelf-ntuples-zskim'))

    for ntupla in ecalelf:
        for f in ntupla.split(','):
            if "ntuple" in f:
                selected_filelist_era.append(f)

    return selected_filelist_era


def findFilesEopruns(runs, rctrl):
    ecalelf = set()
    selected_list = []
    ecalelf.update(rctrl.getOutput(runs = runs, process='ecalelf-ntuples-wskim'))
    ecalelf.update(rctrl.getOutput(runs = runs, process='ecalelf-ntuples-zskim'))
    ecalelf.update(rctrl.getOutput(runs = runs,  process='ecalelf-ntuples'))
    for ntuples in ecalelf:
        for ntuple in ntuples.split(','):
            if "/ntuple_" in ntuple:
                selected_list.append(ntuple)
    return selected_list 

def find_run_numbers(min_run, max_run):
    run_numbers = []
    files_pattern = r"/eos/cms/store/group/dpg_ecal/alca_ecalcalib/automation*/phisym/*/phisymreco_nano_*.root"
    files = glob.glob(files_pattern)
    for i, file in enumerate(files):
        file = file.split(r"/")
        try: 
            run_number = int(file[10])
        except:
            run_number = int(file[9])
        run_numbers.append(run_number)
    runs = np.array(run_numbers)
    runs = runs[(runs <= max_run) & (runs >= min_run)]
    return np.unique(runs).tolist()


#parse arguments
parser = parser = argparse.ArgumentParser()
parser.add_argument("--dbname",          action="store",      type=str,                         help="campaing:: to retrieve files from ecalutomation, see doc")
parser.add_argument("--campaign",        action="store",      type=str,                         help="campaing:: to retrieve files from ecalutomation, see doc")
parser.add_argument("-o", "--outputdir", action="store",      type=str,   default="./",         help="output directory")
parser.add_argument("-l", "--label",     action="store",      type=str,   default="",           help="label to identify the weight file")
parser.add_argument("-v", "--verbosity", action="store",      type=int,   default=1,            help="verbosity level")
parser.add_argument("--run_list",          action="store",    type=str,   default=None,          help="minimum then maximum run numbers:: to retrieve files from ecalutomation, see doc (not yet)")
parser.add_argument("--files_list",      action="store",      type=str,   default=None,          help="list of input file directories:: to retrieve files from ecalutomation, see doc (not yet)" )
args = parser.parse_args()
dbname = args.dbname 
campaign = args.campaign
outputdir =args.outputdir 
label = args.label
run_list = args.run_list
files_list = args.files_list

if run_list is None and files_list is None:
    print("must have either the run number list or a list of phisym nanoADO files")
    assert False
if run_list is not None and files_list is not None:
    print("""only provide one from run_list and files_list:
          make your mind up""")
    assert False

"""
    Get files from ecalautomation.
    Need to define a processor (EcalPhiSymProcessor.py)
    and an accumulator function (EcalPhiSymProcessor.py)
    to read multiple files in coffea and sum above the interesting infos
"""
rctrl = RunCtrl(dbname=dbname, campaign=campaign)




if run_list is not None:
    runs = ast.literal_eval(run_list)
    print("run numbers:", runs)

    phisym_files = {'PhiSym': list(rctrl.getOutput(process='phisym-reco', runs = runs))}
    eop_files = findFilesEopruns(runs, rctrl)
else:
    #using a list of physim files to find the eop files
    phisym_list = files_list.split(',')
    phisym_files = {'Phisym': phisym_list}
    phisym_run_numbers = [phisym_file.split(r"/")[11] for phisym_file in phisym_list]
    ope_files = set()
    for num in phisym_run_numbers:
          #typically multiple ope files for the same phisym file
          file_pattern = f"/eos/cms/store/group/dpg_ecal/alca_ecalcalib/automation_prompt/prompt/ecalelf/*skim/num/ntuple_*.root"
          here = set(glob.glob(file_pattern))
          for ope_file in here: 
              ope_files.add(eop_file)
          
          
          
          
    

if args.verbosity >= 1:
    print("Running on these PhiSym files: ")
    print(phisym_files)
    print("----------------------------------------------------------------------")
    print("Running on these Eop files: ")
    print(eop_files)

#for era in eras:
#    eop_files_conc.append( [eop_file+':selected' for eop_file in eop_files[runs]])
eop = uproot.concatenate(eop_files, filter_name = "/charge|eta/" )

print("done merging")

# reading files through a txt, just for testing purposes
#with open('eopFilesC.txt', 'r') as f:
#    eop_files = f.read().splitlines()
#eop = uproot.concatenate(eop_files, filter_name = "/charge|eta/" )
#
#
#
#import ast
#with open("dictionaryFilesC.txt", "r") as data:
#    phisym_files = ast.literal_eval(data.read())


iterative_run = Runner(
    executor = FuturesExecutor(compression=None, workers=4),
    schema=EcalPhiSymSchema
)
out = iterative_run(
    phisym_files,
    treename="Runs",
    processor_instance=Processor(),
)

idx = ak.argsort(out.run)
ebhits = ak.unflatten(out.ebhits[idx], [len(out.run)], axis=0, behavior=out.ebhits.behavior).sum(axis=1)


if args.verbosity >= 1 : print  ("First doing phisym  histos...")
hist_phisym,bins_phisym =np.histogram(ak.to_numpy(ebhits.ieta[0,:]), weights=ak.to_numpy(ebhits.sumet[0,:]),
           bins=171, range=[-85.5, 85.5])
hist_phisym[85] = 0
plt.bar(np.delete(bins_phisym,171)+0.5, height = hist_phisym, width = 1)
plt.xlabel('i$\eta$')
plt.savefig("%s/phisymEta_%s.png"%(outputdir,label))
plt.clf()
plt.close()


hist_phisym = np.delete(hist_phisym, 85)

if args.verbosity >= 2 :
    print (hist_phisym,bins_phisym)
    print(np.sum(hist_phisym))


if args.verbosity >= 1 : print  ("... then doing Eop histos")
hist_eop,bins_eop =np.histogram(ak.to_numpy(ak.mask(eop.etaEle[:,0], abs(eop.chargeEle[:,0])  == 1  , valid_when=True))/0.0175, bins=171, range=[-85.5, 85.5])

hist_eop[85] = 0
plt.bar(np.delete(bins_eop,171)+0.5, height = hist_eop, width = 1)
plt.xlabel('i$\eta$')
plt.savefig("%s/eopEta_%s.png"%(outputdir,label))
plt.clf()
plt.close()


hist_eop = np.delete(hist_eop, 85)

if args.verbosity >= 2 :
    print (hist_eop,bins_eop)
    print (np.sum(hist_eop))




weights = hist_eop/hist_phisym * (np.sum(hist_phisym) / np.sum(hist_eop) )
plt.bar(np.delete(np.delete(bins_eop,85),170)+0.5, height = weights, width = 1)
plt.ylabel('weights')
plt.xlabel('i$\eta$')
plt.savefig("%s/weightsEta_%s.png"%(outputdir,label))
plt.clf()
plt.close()

if args.verbosity >= 1 : print ("--> saving weights in .txt.")
# save weights ordered from -85 to 85
np.savetxt('%s/weights_%s.txt'%(outputdir,label), weights, delimiter=',')   



plt.hist( weights)
plt.xlabel('weights')
plt.savefig('%s/weights_%s.png'%(outputdir,label))
plt.clf()
plt.close()


