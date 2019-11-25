# CombineHarvester

Full documentation: http://cms-analysis.github.io/CombineHarvester

This is a forked repository. To contribute to the upstream repository, use or branch from the 'upstream' branch. Do not submit a pull request from master!

## Setup Instructions

This package requires HiggsAnalysis/CombinedLimit to be in your local CMSSW area. We follow the [release recommendations](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#for-end-users-that-dont-need-to-commit-or-do-any-development) of the combine developers. The CombineHarvester framework is compatible with the CMSSW 7_4_X, 8_1_X, and 10_2_X series releases. Currently the repository has been designed for 8_1_X, but issues with 10_2_X are expected to be minimal.

A new full release area can be set up and compiled in the following steps:

    export SCRAM_ARCH=slc7_amd64_gcc530
    scram project CMSSW CMSSW_8_1_0
    cd CMSSW_8_1_0/src
    cmsenv
    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    # IMPORTANT: Checkout the recommended tag of CombinedLimit using the instructions here: https://cms-hcomb.gitbooks.io/combine/content/part1/#for-end-users-that-dont-need-to-commit-or-do-any-development
    git clone https://github.com/cms-govner/CombineHarvester.git CombineHarvester
    scram b

## Make the Datacard

The main script is named "DatacardMaker.py". This creates a datacard from a root file of histograms. Currently the histogram is provided with 2lss categories split by lepton flavor, so we must merge these first.

    cd CombineHarvester/TopEFT/test
    python ../scripts/CombineLepFlavors.py # Edit input histogram file
    python ../python/DatacardMaker.py # Edit input (merged) histogram file
    
The resulting datacard will be created in the current directory and named "EFT_MultiDim_Datacard.txt". Renaming is recommended.

## Using asimov or fake data

By default, this script sets the observation to data.

Options:

    --fakedata # By default, uses Asimov dataset. This can be modified to use any WC point in HistogramProcessor.py in __init__\
    --central # Uses central sample histograms for signal expectation

## Other scripts
SpliceCards.py -- Used to combine histograms from different root files.


## Do some fits

For the full suite of fitting scripts, see https://github.com/cms-govner/EFTFit
