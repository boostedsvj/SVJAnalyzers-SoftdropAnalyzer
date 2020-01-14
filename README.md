# SVJAnalyzers/SoftdropAnalyzer

Produces subjets using the softdrop algorithm as implemented in CMSSW, and outputs everything to a flat Ntuple.

## Installation

```
cmsrel CMSSW_9_3_15
cd CMSSW_9_3_15/src
git clone https://github.com/tklijnsma/DataFormats-SVJFormats.git DataFormats/SVJFormats
git clone https://github.com/tklijnsma/SVJAnalyzers-SoftdropAnalyzer.git SVJAnalyzers/SoftdropAnalyzer
scram b
```

The package was developed for `CMSSW_9_3_15`, but should work in most versions >`7_4_X`.

## Running

From `CMSSW_9_3_15/src/SVJAnalyzers/SoftdropAnalyzer`, do:

```
cmsRun python/AnalyzeSoftdrop.py inputFiles=/path/to/input/genfile.root
```
