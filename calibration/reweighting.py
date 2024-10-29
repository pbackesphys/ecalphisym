import awkward as ak
import numpy as np
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
parser.add_argument("-p", "--phisym",    action="store",      type=str,     help="input ntuple for phisymetry")
parser.add_argument("-e", "--eop",    action="store",      type=str,        help="input ntuple for eop")

args = parser.parse_args()

phisymFile = args.phisym
eopFile = args.eop

#phisym = '../Run2018D_test.root'
#eop = '../EGamma-Run2018B-ZSkim-17Sep2018-v1-316998-319312.root'

# path to W/Z E/p - PhiSym 
print  ("Opening Files")
runs = NanoEventsFactory.from_root( phisymFile,
                                   schemaclass=EcalPhiSymSchema,
                                   treepath="/Runs").events()
eop = uproot.open("%s:selected"%eopFile)

counts = np.unique(runs.EcalPhiSymInfo.fill, return_index=True)[1]
splits = np.diff(np.concatenate([counts, [len(runs.EcalPhiSymInfo.fill)]]))

ebhits = ak.unflatten(runs.EcalPhiSymEB, splits, axis=0, behavior=runs.behavior).sum(axis=1)
ak_eop = eop.arrays(filter_name = "/charge|eta/")

print  ("doing phisym  histos")
# ieat histo phisym

print((ebhits.ieta[0, [0, 1, 2, 10000, 10001]]), (ebhits.sumet[0, [0, 1, 2, 10000, 10001]]))
hist_phisym,bins =np.histogram(ak.to_numpy(ebhits.ieta[0,:]), weights=ak.to_numpy(ebhits.sumet[0,:]), 
           bins=171, range=[-85.5, 85.5])
hist_phisym = np.delete(hist_phisym, 85)


print  ("doing Eop histos")

# ieta histo eop
hist_eop,bins =np.histogram(ak.to_numpy(ak.mask(ak_eop.etaEle[:,0], (ak_eop.chargeEle[:,0]  * ak_eop.chargeEle[:,1]  == 0)  , valid_when=False))/0.0175, bins=171, range=[-85.5, 85.5])
hist_eop = np.delete(hist_eop, 85) 

weights = hist_eop/hist_phisym * (np.sum(hist_phisym) / np.sum(hist_eop) )

# save weights ordered from -85 to 85
np.savetxt('weight.txt', weights, delimiter=',')   


