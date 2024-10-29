import os
import awkward as ak
import numpy as np
import pandas as pd
import uproot
from coffea.nanoevents import NanoEventsFactory
from ecalphisym import EcalPhiSymSchema
import argparse 

import matplotlib.pyplot as plt
import mplhep
#%matplotlib inline
plt.rc('font',**{'family':'sans-serif', 'sans-serif':['Helvetica'], 'size':'20'})
plt.rcParams['figure.figsize'] = '10,10'
plt.rcParams['xaxis.labellocation'] = 'right'
plt.rcParams['yaxis.labellocation'] = 'top'
plt.style.use(mplhep.style.ROOT)

#parse arguments
parser = parser = argparse.ArgumentParser()
parser.add_argument("-n", "--naod",    action="store",      type=str,     help="input minbias naod: .root format")
parser.add_argument("-w", "--weights",    action="store",      type=str,     help="input weights: .txt format")
parser.add_argument("-o", "--outputdir",    action="store",      type=str,   default="../results"  ,   help="output directory for ics")
parser.add_argument("--iovref",    action="store",      type=int,   default=0  ,   help="reference iov to be normalized")
parser.add_argument("--merge",  action='store_true' ,  help="merge fill with less hits than a thr")
parser.add_argument("--savePlots",  action='store_true' ,  help="save or not the plots")
parser.add_argument("--saveICs",   action='store_true' ,   help="save or not the ics") 
parser.add_argument('--fileFormats' , dest='fileFormats'    ,  default='png', help='Output plot file formats (comma-separated png, pdf). Default "png"')

args = parser.parse_args()

naod = args.naod # '../Run2018D_test.root'
outputdir =args.outputdir #"../ICs_testing"
saveICs = args.saveICs
savePlots = args.savePlots
merge = args.merge
weights = args.weights
fileFormats = args.fileFormats.split(',')
iovref = args.iovref

# replace with proper test file location
runs = NanoEventsFactory.from_root(naod,
                               schemaclass=EcalPhiSymSchema,
                               treepath="/Runs").events()


# create the outputdir
os.makedirs(outputdir, exist_ok=True)   

## Group data by fill number
#    - Prepare the grouping 
#    - Perform the grouping using awkward.unflatten separately on each main collection

counts = np.unique(runs.EcalPhiSymInfo.fill, return_index=True)[1]
splits = np.diff(np.concatenate([counts, [len(runs.EcalPhiSymInfo.fill)]]))

info = ak.unflatten(runs.EcalPhiSymInfo, splits, axis=0, behavior=runs.behavior).sum(axis=1)
ebhits = ak.unflatten(runs.EcalPhiSymEB, splits, axis=0, behavior=runs.behavior).sum(axis=1)
eehits = ak.unflatten(runs.EcalPhiSymEE, splits, axis=0, behavior=runs.behavior).sum(axis=1)

niovs = len(splits)

# try to put togheter fills with less stats
if merge:
    nhitsTOT = ak.sum(ebhits.nhits, axis=1)
    
    it = splits
    merged_splits = []
    nthr = 100000000
    
    it = iter(it)
    i = 0
    while True:    
    
        try:
            current = next (it)
        except StopIteration as e:
            print(e)           
            break      
                    
        while nhitsTOT[i] < nthr:
            try:
    
    
                current += next(it)
                i = i +1
    
            except StopIteration as e:
                print(e)           
                break
        
        merged_splits.append(current)   
        i = i +1
    
        
    splits = merged_splits
    niovs = len(splits)
    
    info = ak.unflatten(runs.EcalPhiSymInfo, splits, axis=0, behavior=runs.behavior).sum(axis=1)
    ebhits = ak.unflatten(runs.EcalPhiSymEB, splits, axis=0, behavior=runs.behavior).sum(axis=1)
    eehits = ak.unflatten(runs.EcalPhiSymEE, splits, axis=0, behavior=runs.behavior).sum(axis=1)



## Compute k-factors
#    - keep same definition as in Run2: the miscalibration values are centered at 0.

k = ak.linear_fit(info.miscalibs_eb, ebhits.sumet_v, axis=2)

# K-factors plots examples
#    - all k-factors (slopes from the fits)
#    - all intercepts
#    - 1 channel k-factor history
#    - 1 iov k-factors map

# k-factor histo
plt.hist(ak.flatten(k.slope), bins=1000, range=[0,4])
plt.xlabel('k-factor')
plt.clf()
plt.close()



plt.hist(ak.flatten(k.intercept), bins=1000, range=[-0.1, 0.1])
plt.xlabel('fit intercept')
plt.clf()
plt.close()



# k-factor vs fill
plt.scatter(info.fill, k.slope[:,100])
plt.xlabel('Fill Number')
plt.ylabel('k-factor')
if savePlots: 
    for form in fileFormats:
        plt.savefig("%s/kfactorVSfill.%s"%(outputdir, form))
plt.clf()
plt.close()


# k-factor map
plt.hist2d(ak.to_numpy(ebhits.iphi[iovref,:]), ak.to_numpy(ebhits.ieta[iovref,:]), weights=ak.to_numpy(k.slope[iovref]), 
           bins=[360, 171], range=[[0.5,360.5], [-85.5, 85.5]], 
           cmap='cividis', cmin=2, cmax=3)
plt.xlabel('i $\phi$')
plt.ylabel('i $\eta$')
plt.title('k-factor')
if savePlots: 
    for form in fileFormats:
        plt.savefig("%s/map_kfactor.%s"%(outputdir, form))
plt.clf()
plt.close()


## Compute EFlow
#    - remove boundary channels (SM and module boundaries) from the total EB sum
#    - Note: double check this definition with the Run2 one
#    - Note: there was an issue with the first fill in the test file used for this example, hence the analysis starts from iov 1
    


def boundaryCrystals(data):
    """
    Flag crystals on module boundaries:
    - first and last crystals in a SM along phi (iphi % 20 == 0|1)
    - first and last crystals in a module along eta (|ieta| = 1, 25,26, 45,46, 65,66, 85)
    """

    bounds = ak.zeros_like(data.ieta)
    for idx in [1, 25, 25, 45, 46, 65, 66, 85]:
        bounds = bounds + (abs(data.ieta) == idx)
        
    return (data.iphi % 20 == 0) | (data.iphi % 20 == 1) | (bounds > 0)

### Weights to match eta ele distribution

with open(weights) as file:
    weights = np.loadtxt(file)

# repeated for 360 xstals and for the number of iovs
ws = ak.flatten(ak.Array([[w] * 360 for w in weights]))
ws = ak.Array([ws] * niovs) 

sumEtEB = ak.sum(ak.mask(ebhits.sumet, boundaryCrystals(ebhits), valid_when=False), axis=1)
sumEtEBw = ak.mean(ak.mask(ebhits.sumet, boundaryCrystals(ebhits), valid_when=False), weight = ws ,  axis=1)

### Eflow ICs normalized to ref IOV

norm = ak.Array(np.repeat([ebhits.sumet[iovref]/sumEtEB[iovref]], niovs, axis=0)) 
eflow = (((ebhits.sumet/sumEtEB)/norm)-1)/k.slope+1

normW = ak.Array(np.repeat([ebhits.sumet[iovref]/sumEtEBw[iovref]], niovs, axis=0)) 
eflowW = (((ebhits.sumet/sumEtEBw)/normW)-1)/k.slope+1
print(len(ebhits.sumet[0]))

## EFlow plots examples

# plot of the laser corrections and eflow ics for some xstals
for ixstal in range(0,62100,10000):
    plt.clf()
    plt.scatter(info.fill[iovref:], eflow[iovref:,ixstal], label='ICs EFlow')
    plt.scatter(info.fill[iovref:], eflowW[iovref:,ixstal], label='ICs EFlow W')
    plt.scatter(info.fill[iovref:], (ebhits.sumlc[iovref,ixstal]/ebhits.nhits[iovref,ixstal])/(ebhits.sumlc[iovref:,ixstal]/ebhits.nhits[iovref:,ixstal]), 
            label='1/Laser correction')
    plt.legend()
    plt.xlabel('Fill number')
    plt.ylabel('Relative response')
    plt.title('i$\eta$ = %i i$\phi$ = %i'%(int(ebhits.ieta[iovref,ixstal]), ebhits.iphi[iovref,ixstal]))
    plt.grid()
    if savePlots: 
        for form in fileFormats:
            plt.savefig('%s/monitoring_ieta%i_iphi%i.%s'%(outputdir, ebhits.ieta[iovref,ixstal], ebhits.iphi[iovref,ixstal], form))
    plt.clf()
    plt.close()

# plot the map of the last iov
plt.hist2d(ak.to_numpy(ebhits.iphi[-1,:]), ak.to_numpy(ebhits.ieta[-1,:]), weights=ak.to_numpy(eflow[-1]), 
           bins=[360, 171], range=[[0.5,360.5], [-85.5, 85.5]], 
           cmap='seismic', cmin=0.9, cmax=1.1)
plt.xlabel('i $\phi$')
plt.ylabel('i $\eta$')
plt.title('EFlow ICs')
plt.colorbar()
if savePlots: 
    for form in fileFormats:
        plt.savefig('%s/mapICs_fill%s.%s'%(outputdir, str(info.fill[-1]).replace("[","").replace("]",""), form))
plt.clf()
plt.close()

# plot the map of the last iov eflow Weighted
plt.hist2d(ak.to_numpy(ebhits.iphi[-1,:]), ak.to_numpy(ebhits.ieta[-1,:]), weights=ak.to_numpy(eflowW[-1]), 
           bins=[360, 171], range=[[0.5,360.5], [-85.5, 85.5]], 
           cmap='seismic', cmin=0.9, cmax=1.1)
plt.xlabel('i $\phi$')
plt.ylabel('i $\eta$')
plt.title('EFlow ICs')
plt.colorbar()

if savePlots: 
    fillnum = str(info.fill[-1]).replace("[","").replace("]","")
    for form in fileFormats:
        plt.savefig('%s/mapICsW_fill%s.%s'%(outputdir, fillnum, form))
plt.clf()
plt.close()


### ICs saved in the correct format

# FORMAT: ieta - iphi - zside - eflow - error NB now error is 0
# One folder is created for each IOV
# Passing trough PD
if saveICs:
    for ifill in range(iovref,len(info.fill)):
        tosave = pd.concat ([ak.to_pandas(ebhits[ifill].ieta).astype(int), ak.to_pandas(ebhits[ifill].iphi), ak.to_pandas(ebhits[ifill].zside()).astype(int), ak.to_pandas(eflowW[ifill])], axis = 1)
        tosave[''] = 0 
        tosave.fillna(1, inplace=True)
        os.makedirs('%s/%s'%(outputdir, str(info.fill[ifill]).replace("[","").replace("]","")), exist_ok=True)   
        tosave.to_csv('%s/%s/file.txt'%(outputdir, str(info.fill[ifill]).replace("[","").replace("]","")), float_format='%.6f', index=False, header=False)
        

        

