#!/bin/bash

#the initial terminal directory location is <user directory location>/ecalphisym/calibration

# Source the conda initialization script
source /afs/cern.ch/user/p/pbackes/.bashrc

# Activate the phisym conda environment
conda activate /afs/cern.ch/user/p/pbackes/.conda/envs/phisym

# Run the Python script
python3 eflow_processor.py --dbname ecal_rerecos_v1 --campaign test_eflow_v1 --run_list 386642 -w weight.txt -o ../../www/eflow_processor_out --savePlots --verbosity 2
