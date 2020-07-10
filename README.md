Description 
============

Installation
------------

Current CMSSW version: ``CMSSW_10_2_16_UL``.

Get a supported CMSSW release area:

```bash
  scram pro -n MyWorkingAreaName CMSSW <CMSSW_VERSION>
  cd MyWorkingAreaName/src
  # Setup your CMSSW environment
  cmsenv
  # SSH agent is optional, but will save you from typing your password many times
  eval `ssh-agent -s`
  ssh-add
  # Run this before doing ANYTHING else in src
  git cms-init
```

Checkout the FinalStateAnalysis repository:

```bash
  git clone --recursive https://github.com/cecilecaillol/B_mumupi_analyzer.git
  cd B_mumupi_analyzer
```


