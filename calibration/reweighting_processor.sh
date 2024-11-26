#!/bin/bash

#$RUNS = [374997,  375665,  379356,  380306,  380945,  382225,  382834,  383537]
#the initial terminal directory location is <user directory location>/ecalphisym/calibration

# Source the conda initialization script
source /afs/cern.ch/user/p/pbackes/.bashrc

# Activate the phisym conda environment
conda activate /afs/cern.ch/user/p/pbackes/.conda/envs/phisym

# Change directory to the location required by the python script
cd ../../automation-control/

# Run the Python script
python3 ../ecalphisym/calibration/reweighting_processor.py --dbname ecal_prompt_v2 --campaign prompt --minrun 382834 --maxrun 382840 -o ../ecalphisym/reweighting_output --verbosity 1
