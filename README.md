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

The main script is named "Top_EFT_Maker.py". To run it:

    cd CombineHarvester/TopEFT/test
    python ../scripts/Top_EFT_Maker.py
    
The resulting datacard will be created in your current directory (.../test/)

## Do some fits

These datacards can be fit with the usual combine commands, but we have a script to do a simulataneous 4D fit on the four signal processes (ttZ, ttW, ttH, tZq). It first makes a workspace, then does a MultiDimFit. To run it (from "test"):

    source ../scripts/4DFit.csh

Note there are two places toward the end of the output that give the best fit mus. This is because an increase in verbosity is required to spit out the uncertainties.
