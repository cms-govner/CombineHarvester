# CombineHarvester

Full documentation: http://cms-analysis.github.io/CombineHarvester

This is a forked repository. To contribute to the upstream repository, use or branch from the 'upstream' branch.

## Setup Instructions

This package requires HiggsAnalysis/CombinedLimit to be in your local CMSSW area. We follow the [release recommendations](https://cms-hcomb.gitbooks.io/combine/content/part1/#for-end-users-that-dont-need-to-commit-or-do-any-development) of the combine developers . The CombineHarvester framework is compatible with the CMSSW 7_4_X and 8_1_X series releases.

A new full release area can be set up and compiled in the following steps:

    export SCRAM_ARCH=slc6_amd64_gcc530
    scram project CMSSW CMSSW_8_1_0
    cd CMSSW_8_1_0/src
    cmsenv
    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    # IMPORTANT: Checkout the recommended tag of CombinedLimit using the instructions here: https://cms-hcomb.gitbooks.io/combine/content/part1/#for-end-users-that-dont-need-to-commit-or-do-any-development
    git clone https://github.com/cms-govner/CombineHarvester.git CombineHarvester
    scram b

## Make the Datacard

The main script is named "TopEFT_DatacardMaker.py". This creates a datacard from our ROOT trees. To run it:

    cd CombineHarvester/TopEFT/test
    python ../scripts/TopEFT_DatacardMaker.py
    
The resulting datacard will be created in the current directory (.../test/) and named "EFT_MultiDim_Datacard.txt"

## Using asimov or fake data

By default, this datacard sets the observation to be fake data for 16 EFT operators set to 1.0. To use fewer operators in the fake data, swap L32-33 of Process_Input.py. To use SM asimov data instead, add 0 as an argument to the TopEFT_DatacardMaker.py call.

## Do some fits

For the full suite of fitting scripts, see https://github.com/cms-govner/EFTFit

To quickly transfer the datacard to that diretory for fitting, use:

    source ../scripts/transfer.csh
