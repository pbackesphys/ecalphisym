# ECAL PhiSym calibration package

This package provides tools and examples to perform phisym/eflow calibration of the CMS
ECAL detector using ZeroBias events. This version works with reconstructed events 
processed with CMSSW and saved in NanoAOD format. Please read the compiled [documentation](https://spigazzi.web.cern.ch/docs/ecalphisym/)

## Installation
This package can be installed from github using pip:

```
pip install git+https://github.com/simonepigazzini/ecalphisym.git
```
### Conda environment
An easy way to install and test is to create a dedicated conda environment

```
conda create -n phisym python==3.9
conda activate phisym
pip install git+https://github.com/simonepigazzini/ecalphisym.git
```

Moreover, you need the ecal-automation package to read the phisym nAOD with the pytools:
```
pip install git+https://gitlab.cern.ch/cms-ecal-dpg/ECALELFS/automation-control.git
```

Then just get the calibration codes (in this repo calibration folder) and run them.
Below instructions to run the calibration codes.


### Using docker image 
#### Recommended

Clone the github repository to get the calibration code.

```
git clone git@github.com:flaviacetorelli/ecalphisym.git
```

Then set up the docker image and run the code in the image.
To see how to run the calibration code, read below.


```
set -x APPTAINER_CACHEDIR "/tmp/"(whoami)"/singularity"
singularity shell -B /eos -B /cvmfs -B /afs docker://gitlab-registry.cern.ch/cms-ecal-dpg/ecalelfs/automation:dev
source /home/ecalgit/setup.sh
```

## Running calibrations
Codes are in the calibration folder.

First python script reweighting_processor.py is used to calculate the eta-weight to compare and match the eta-distribution of the MinBias events to the eta distribution of the electrons from W-Z.
No needed in case you only want to study the history of individual crystals.

To run:
```
python3 reweighting_processor.py --dbname ecal_prompt_v1 --campaign prompt --eras Run2022C -o OUTPUTDIR -l LABEL
```

where dbname, campaign, and eras are needed to catch files with the pytools, while -o specify the output directory and -l the label for the weight output file.

The actual script computing the calibration and saving them in .txt and .root file is eflow_processor.py.
To run it:
```
python3 eflow_processor.py --dbname ecal_prompt_v1 --campaign prompt --eras Run2022C,Run2022D -w PATH/weights.txt -o OUTPUTDIR 
```
Where dbname, campaign, and eras are needed to catch files with the pytools, -w specify the path weight file from the previous step (if not given, the eflow ics are calculated without the eta-weights), -o the output directory.
Moreover, --savePlots option saves the history plots of some crystals to check everything is fine.



