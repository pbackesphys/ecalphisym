#!/bin/bash

#the initial terminal directory location is <user directory location>/ecalphisym/calibration

# Source the conda initialization script
source /afs/cern.ch/user/p/pbackes/.bashrc

# Activate the phisym conda environment
conda activate /afs/cern.ch/user/p/pbackes/.conda/envs/phisym

# Change directory to the location required by the python script
cd ../../automation-control/

# Run the Python script
python3 ../ecalphisym/calibration/reweighting_processor.py --dbname ecal_prompt_v2 --campaign prompt --eras Run2024B -o ../ecalphisym/reweighting_output --verbosity 1
