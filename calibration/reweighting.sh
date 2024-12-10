#!/bin/bash

#the initial terminal directory location is <user directory location>/ecalphisym/calibration

# Source the conda initialization script
source /afs/cern.ch/user/p/pbackes/.bashrc

# Activate the phisym conda environment
conda activate /afs/cern.ch/user/p/pbackes/.conda/envs/phisym

# Change directory to the location required by the python script
cd ..

# Run the Python script
python3 calibration/reweighting.py --phisym /eos/cms/store/group/dpg_ecal/alca_ecalcalib/automation_prompt/prompt/phisym/379003/phisymreco_nano_0.root --eop /eos/cms/store/group/dpg_ecal/alca_ecalcalib/automation_prompt/prompt/ecalelf/wskim/379003/ntuple_0.root