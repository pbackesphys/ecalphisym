import awkward as ak
import numpy as np
import uproot
from ecalautoctrl import RunCtrl
from coffea.nanoevents import NanoEventsFactory
from ecalphisym import EcalPhiSymSchema
from coffea.processor import AccumulatorABC, ProcessorABC, Runner, FuturesExecutor
from EcalPhiSymProcessor import EcalPhiSymAccumulator, Processor
import argparse
import os,sys


import matplotlib.pyplot as plt
import mplhep
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






#parse arguments
parser = parser = argparse.ArgumentParser()
parser.add_argument("--dbname",          action="store",      type=str,                         help="campaing:: to retrieve files from ecalutomation, see doc")
parser.add_argument("--campaign",        action="store",      type=str,                         help="campaing:: to retrieve files from ecalutomation, see doc")
parser.add_argument("--eras",            action="store",      type=str,                         help="eras:: to retrieve files from ecalutomation, see doc")
parser.add_argument("-o", "--outputdir", action="store",      type=str,   default="./",         help="output directory")
parser.add_argument("-l", "--label",     action="store",      type=str,   default="",           help="label to identify the weight file")
parser.add_argument("-v", "--verbosity", action="store",      type=int,   default=1,            help="verbosity level")

args = parser.parse_args()
dbname = args.dbname 
campaign = args.campaign 
eras = args.eras.split(',')
outputdir =args.outputdir 
label = args.label

"""
    Get files from ecalautomation.
    Need to define a processor (EcalPhiSymProcessor.py)
    and an accumulator function (EcalPhiSymProcessor.py)
    to read multiple files in coffea and sum above the interesting infos
"""
rctrl = RunCtrl(dbname=dbname, campaign=campaign)

phisym_files =  {}
eop_files = {}
eop_files_conc = []

for era in eras:
    phisym_files[era] = [f for f in rctrl.getOutput(era=era, process='phisym-reco')]
    eop_files [era] = findFilesEop(rctrl,  era)
if args.verbosity >= 1:
    print("Running on these PhiSym files: ")
    print(phisym_files)
    print("----------------------------------------------------------------------")
    print("Running on these Eop files: ")
    print(eop_files)

for era in eras:
    eop_files_conc.append( [eop_file+':selected' for eop_file in eop_files[era]])
eop = uproot.concatenate(eop_files_conc, filter_name = "/charge|eta/" )



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


