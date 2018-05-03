# CombineHarvester

Full documentation: http://cms-analysis.github.io/CombineHarvester

This is a forked repository. To contribute to the upstream repository, use or branch from the 'upstream' branch.

## Setup Instructions

This package requires HiggsAnalysis/CombinedLimit to be in your local CMSSW area. We follow the release recommendations of the combine developers which can be found [here](https://cms-hcomb.gitbooks.io/combine/content/part1/#for-end-users-that-dont-need-to-commit-or-do-any-development). The CombineHarvester framework is compatible with the CMSSW 7_4_X and 8_1_X series releases.

A new full release area can be set up and compiled in the following steps:

    export SCRAM_ARCH=slc6_amd64_gcc530
    scram project CMSSW CMSSW_8_1_0
    cd CMSSW_8_1_0/src
    cmsenv
    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    # IMPORTANT: Checkout the recommended tag of [CombinedLimit](https://cms-hcomb.gitbooks.io/combine/content/part1/#for-end-users-that-dont-need-to-commit-or-do-any-development)
    git clone https://github.com/cms-govner/CombineHarvester.git CombineHarvester
    scram b

