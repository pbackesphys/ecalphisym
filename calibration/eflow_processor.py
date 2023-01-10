import awkward as ak
import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep
from coffea.nanoevents import NanoEventsFactory
from ecalphisym import EcalPhiSymSchema
from ecalautoctrl import RunCtrl
from coffea.processor import AccumulatorABC, ProcessorABC, Runner, FuturesExecutor
import argparse
import os,sys
import time
from EcalPhiSymProcessor import EcalPhiSymAccumulator, Processor
#%matplotlib inline
plt.rc('font',**{'family':'sans-serif', 'sans-serif':['Helvetica'], 'size':'20'})
plt.rcParams['figure.figsize'] = '10,10'
plt.rcParams['xaxis.labellocation'] = 'right'
plt.rcParams['yaxis.labellocation'] = 'top'
plt.style.use(mplhep.style.ROOT)



#parse arguments
parser = parser = argparse.ArgumentParser()
parser.add_argument("--dbname",          action="store",      type=str,                         help="dbname:: to retrieve files from ecalutomation, see doc")  
parser.add_argument("--campaign",        action="store",      type=str,                         help="campaing:: to retrieve files from ecalutomation, see doc")  
parser.add_argument("--eras",            action="store",      type=str,                         help="eras:: to retrieve files from ecalutomation, see doc")  
parser.add_argument("-w", "--weights",   action="store",      type=str,                         help="input weights: .txt format")
parser.add_argument("--iovref",          action="store",      type=int,   default=0,            help="reference iov to be normalized")
parser.add_argument("--nhits",           action="store",      type=int,   default=3000000000,   help="value of minimum nhits per iov")
parser.add_argument("-o", "--outputdir", action="store",      type=str,   default="../results", help="output directory for ics")
parser.add_argument("--savePlots",       action="store_true",                                   help="save some plots")
parser.add_argument("-v", "--verbosity", action="store",      type=int,   default=1,            help="verbosity level")

args = parser.parse_args()

dbname = args.dbname # '../Run2018D_test.root'
campaign = args.campaign # '../Run2018D_test.root'
eras = args.eras.split(',')
outputdir =args.outputdir #"../ICs_testing"
weights = args.weights
iovref = args.iovref
nhitsmin = args.nhits

# create the outputdir
os.makedirs(outputdir, exist_ok=True)   


"""
    Get files from ecalautomation.
    Need to define a processor (EcalPhiSymProcessor.py)
    and an accumulator function (EcalPhiSymProcessor.py)
    to read multiple files in coffea and sum above the interesting infos
"""

fileset =  {}
rctrl = RunCtrl(dbname=dbname, campaign=campaign)


for era in eras:
    fileset[era] = [f for f in rctrl.getOutput(era=era, process='phisym-reco')]

if args.verbosity >= 1:
   print("Running on these files: ")
   print(fileset)

#Temporary needed to avoid crashing due to not found files
if "Run2022D" in eras:
    fileset["Run2022D"].remove('/eos/cms/store/group/dpg_ecal/alca_ecalcalib/automation_prompt/phisym/357889/phisymreco_nano_0.root')

# parsing data from a txt (just for testing)
#import ast
#with open("dictionaryFilesC.txt", "r") as data:
#    fileset = ast.literal_eval(data.read())
#fileset = {'PhiSym': ['/eos/cms/store/group/dpg_ecal/alca_ecalcalib/automation_prompt/phisym/356489/phisymreco_nano_1.root', 
#                      '/eos/cms/store/group/dpg_ecal/alca_ecalcalib/automation_prompt/phisym/356951/phisymreco_nano_0.root', 
#                      '/eos/cms/store/group/dpg_ecal/alca_ecalcalib/automation_prompt/phisym/356969/phisymreco_nano_0.root',
#                      '/eos/cms/store/group/dpg_ecal/alca_ecalcalib/automation_prompt/phisym/356469/phisymreco_nano_0.root']}
#print( fileset)

iterative_run = Runner(
    executor = FuturesExecutor(compression=None, workers=4),
    schema=EcalPhiSymSchema
)
out = iterative_run(
    fileset,
    treename="Runs",
    processor_instance=Processor(),
)

runs = out.run
hitseb = ak.flatten(out.info.hitseb)
fills = ak.flatten(out.info.fill)
idx = ak.argsort(runs)

# requiring at least nhistmin for each IOV
runs_merged_i = []
runs_merged_f = []
counts_merged = [0]
nhits = 0
print (len(hitseb[idx]))
for i, hits in enumerate(hitseb[idx]):
    if nhits == 0 : runs_merged_i.append(runs[idx][i])
    if args.verbosity >= 2: print( "Fill: %i and Run: %i"%( fills[idx][i], runs[idx][i]),  "\n hits here %i -- n hits before %i"%(hits , nhits))
    nhits += hits

    if nhits >= nhitsmin:
       if args.verbosity >= 2: print (">>>>> At idx %i --> min hits achieved == %i"%(i, nhits))
       nhits = 0

       if args.verbosity >= 2:  print ("Now closing IOV and resetting to -->= ", nhits )
       counts_merged.append(i+1)
       runs_merged_f.append(runs[idx][i])

counts_merged = np.array(counts_merged)
splits_merged = np.diff(np.concatenate([ counts_merged, [len(runs)]]))
splits_merged = splits_merged [splits_merged != 0]

if args.verbosity >= 2:  print ("So this are the counts merged:: \n  ", counts_merged , "\n", "These are the final splits: \n ", splits_merged)


splits = splits_merged
niovs = len (splits_merged)

if len(runs_merged_f) != len(runs_merged_i): runs_merged_f.append(runs[idx][-1])

if args.verbosity >= 1: 
    print ("Diveded in %i IOVS --> These are the initial runs of the IOV: \n"%len(runs_merged_i), runs_merged_i)
    print ("Diveded in %i IOVS --> These are the final runs of the IOV: \n"%len(runs_merged_f), runs_merged_f)

# Merging the info, ebhits, eehits per  
info = ak.unflatten(out.info[idx], splits, axis=0, behavior=out.info.behavior).sum(axis=1)
ebhits = ak.unflatten(out.ebhits[idx], splits, axis=0, behavior=out.ebhits.behavior).sum(axis=1)
eehits = ak.unflatten(out.eehits[idx], splits, axis=0, behavior=out.eehits.behavior).sum(axis=1)


if args.verbosity >= 2: print ("Number of hits for each IOV: ", info.hitseb)

# Calculating the k - factor and plotting
k = ak.linear_fit(info.miscalibs_eb, ebhits.sumet_v, axis=2)

if args.savePlots:
    plt.hist(ak.flatten(k.slope), bins=1000, range=[0,4])
    plt.xlabel('k-factor Slope')
    plt.savefig("%s/kSlope.png"%outputdir)
    plt.clf()
    plt.close()
    
    plt.hist(ak.flatten(k.intercept), bins=1000, range=[-0.1, 0.1])
    plt.savefig("%s/kIntercept.png"%outputdir)
    plt.xlabel('k-factor Intercept')
    plt.clf()
    plt.close()
    
    
    plt.scatter(runs_merged_i[iovref:], k.slope[:,100])
    plt.xlabel('Runs')
    plt.xlabel('K factor')
    plt.savefig("%s/kslopeVSrun.png"%outputdir)
    plt.clf()
    plt.close()
    
    plt.hist2d(ak.to_numpy(ebhits.iphi[iovref,:]), ak.to_numpy(ebhits.ieta[iovref,:]), weights=ak.to_numpy(k.slope[iovref]), 
               bins=[360, 171], range=[[0.5,360.5], [-85.5, 85.5]], 
               cmap='viridis', cmin=2, cmax=3)
    
    plt.colorbar()
    plt.savefig("%s/map_kslope.png"%outputdir)
    plt.clf()
    plt.close()


def boundaryCrystals(data):
    """
    Flag crystals on module boundaries:
    - first and last crystals in a SM along phi (iphi % 20 == 0|1)
    - first and last crystals in a module along eta (|ieta| = 1, 25,26, 45,46, 65,66, 85)
    """

    bounds = ak.zeros_like(data.ieta)
    for idx in [1, 25, 26, 45, 46, 65, 66, 85]:
        bounds = bounds + (abs(data.ieta) == idx)
        
    return (data.iphi % 20 == 0) | (data.iphi % 20 == 1) | (bounds > 0)

# Opening the files with the eta-weights
with open(weights) as file:
    weights = np.loadtxt(file)

# repeated for 360 xstals and for the number of iovs
ws = ak.flatten(ak.Array([[w] * 360 for w in weights]))
ws = ak.Array([ws] * niovs) 


#Deriving the corrections
if args.verbosity >= 1: print ("Deriving the corrections...") 

######## IMPLEMENTIG normalization per Ring
#sumEtRing = []
#for iRing in range(-85,86):
#    sumEtRing.append( ak.sum(ak.mask(ebhits.sumet, ebhits.ieta == iRing, valid_when=True), axis=1)) # sum 
#
#sumEtRing.pop(85)
#print()
#
#print (sumEtRing)
##print (sumEtRing)
#print (len(sumEtRing))
#
#
#sys.exit()
#norm_ring = ak.Array(np.repeat([ebhits.sumet[iovref]/sumEtRing[iovref]], len(ebhits.sumet), axis=0))
#phisym = (((ebhits.sumet/sumEtRing)/norm_ring)-1)/k.slope+1
#
##        icChEB[index] = 1/((ebXstals_[index].GetSumEt(0)/
##                            ebRingsSumEt_[currentRing][0]-1)/GetChannelKfactor(index, 0).first+1);
#
#norm = ak.Array(np.repeat([ebhits.sumet[iovref]/sumEtEB[iovref]], len(ebhits.sumet), axis=0))
#eflow = (((ebhits.sumet/sumEtEB)/norm)-1)/k.slope+1



sumEtEB = ak.sum(ak.mask(ebhits.sumet, boundaryCrystals(ebhits), valid_when=False), axis=1) # sum 
norm = ak.Array(np.repeat([ebhits.sumet[iovref]/sumEtEB[iovref]], len(ebhits.sumet), axis=0))
eflow = (((ebhits.sumet/sumEtEB)/norm)-1)/k.slope+1


sumEtEB_w = ak.mean(ak.mask(ebhits.sumet, boundaryCrystals(ebhits), valid_when=False), weight = ws ,  axis=1) # weighted average
norm_w = ak.Array(np.repeat([ebhits.sumet[iovref]/sumEtEB_w[iovref]], niovs, axis=0)) 
eflow_w = (((ebhits.sumet/sumEtEB_w)/norm_w)-1)/k.slope+1

start = time.time()
if args.verbosity >= 1: print ("... saving histories in .root ...") 

######### Saving in a ROOT file
with  uproot.recreate("%s/history.root"%outputdir) as file:

    for iov, run in enumerate(runs_merged_i):
        file["history"] = {
            "ieta"    : [ebhits.ieta[iov,:]   for iov in range(0,niovs)],
            "iphi"    : [ebhits.iphi[iov,:]   for iov in range(0,niovs)],
            "run_i"   : [runs_merged_i[iov]   for iov in range(0,niovs)],
            "run_f"   : [runs_merged_f[iov]   for iov in range(0,niovs)],
            "eflow"   : [eflow[iov,:]         for iov in range(0,niovs)],
            #"eflow_w" : [eflow_w[iov,:]       for iov in range(0,niovs)],
    
    }

end = time.time()

if args.verbosity >= 2: print ("How much does it takes???", (end - start) , "s")


#Plotting for some crystals
if args.savePlots:
    if args.verbosity >= 1: print ("... plotting ...") 
    for icry in range(0,61200,100):
        #icry=500
        if (ebhits.status[iovref,icry] != 0): 
            if args.verbosity >= 1: print ("Bad xstal i$\eta$ = %i i$\phi$ = %i ---> Skipping "%(ebhits.ieta[iovref,icry],ebhits.iphi[iovref,icry]) )
            continue
        plt.scatter(runs_merged_i[iovref:], eflow[iovref:,icry], label='EFlow')
        plt.scatter(runs_merged_i[iovref:], eflow_w[iovref:,icry], label='EFlow weighted')
        plt.scatter(runs_merged_i[iovref:], (ebhits.sumlc[iovref,icry]/ebhits.nhits[iovref,icry])/(ebhits.sumlc[iovref:,icry]/ebhits.nhits[iovref:,icry]), 
                    label='1/Laser correction')
        plt.legend()
        plt.grid()
        plt.xlabel('Intial run number of IOV')
        plt.ylabel('Relative response')
        plt.title('Crystal i$\eta$ = %i i$\phi$ = %i'%(ebhits.ieta[iovref,icry],ebhits.iphi[iovref,icry]))
        plt.savefig("%s/history_IEta_%i_IPhi_%i.png"%(outputdir,ebhits.ieta[iovref,icry], ebhits.iphi[iovref,icry]))
        plt.clf()
        plt.close()
    
    
    
    # plot the map of the last run
    plt.hist2d(ak.to_numpy(ebhits.iphi[-1,:]), ak.to_numpy(ebhits.ieta[-1,:]), weights=ak.to_numpy(eflow[-1]), 
               bins=[360, 171], range=[[0.5,360.5], [-85.5, 85.5]], 
               cmap='viridis', cmin=0.9, cmax=1.1)
    plt.colorbar()
    plt.xlabel('i $\phi$')
    plt.ylabel('i $\eta$')
    plt.title('EFlow ICs')
    plt.savefig("%s/ICmap2d_Run%i.png"%(outputdir,runs_merged_i[-1]))
    plt.clf()
    plt.close()
    
    # plot the map of the last run
    plt.hist2d(ak.to_numpy(ebhits.iphi[-1,:]), ak.to_numpy(ebhits.ieta[-1,:]), weights=ak.to_numpy(eflow_w[-1]), 
               bins=[360, 171], range=[[0.5,360.5], [-85.5, 85.5]], 
               cmap='viridis', cmin=0.9, cmax=1.1)
    plt.colorbar()
    plt.xlabel('i $\phi$')
    plt.ylabel('i $\eta$')
    plt.title('EFlow ICs weighted')
    plt.savefig("%s/ICweighted_map2d_Run%i.png"%(outputdir,runs_merged_i[-1]))
    plt.clf()
    plt.close()




if args.verbosity >= 1: print ("... saved plots and now ICs in %s"%outputdir) 

### ICs saved in the correct format
import pandas as pd
# FORMAT: ieta - iphi - zside - eflow - error NB now error is 0
# One folder is created for each IOV
# Passing trough PD
for i,irun in enumerate(runs_merged_i):
     tosave = pd.concat ([ak.to_pandas(ebhits[i].ieta).astype(int), ak.to_pandas(ebhits[i].iphi), ak.to_pandas(ebhits[i].zside()).astype(int), ak.to_pandas(eflow_w[i])], axis = 1)
     tosave[''] = 0 
     tosave.fillna(1, inplace=True)
     os.makedirs('%s/%s'%(outputdir, str(irun).replace("[","").replace("]","")), exist_ok=True)   
     tosave.to_csv('%s/%s/file.txt'%(outputdir, str(irun).replace("[","").replace("]","")), float_format='%.6f', index=False, header=False)
        

        

