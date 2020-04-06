# SVJAnalyzers/SoftdropAnalyzer

Produces subjets using the softdrop algorithm as implemented in CMSSW, and outputs everything to a flat Ntuple.

## Installation

```
cmsrel CMSSW_X_Y_Z
cd CMSSW_X_Y_Z/src
git clone https://github.com/tklijnsma/DataFormats-SVJFormats.git DataFormats/SVJFormats
git clone https://github.com/tklijnsma/SVJAnalyzers-SoftdropAnalyzer.git SVJAnalyzers/SoftdropAnalyzer
scram b
```

The package was tested for `CMSSW_9_3_15` (2017) and `CMSSW_10_2_18` (2018), but should work in most versions >`7_4_X`.

## Running


# Softdrop analyzer

From `CMSSW_X_Y_Z/src/SVJAnalyzers/SoftdropAnalyzer`, do:

```
cmsRun python/AnalyzeSoftdrop.py inputFiles=file:/path/to/input/genfile.root
```

# Saving detailed information about the Z'

From `CMSSW_X_Y_Z/src/SVJAnalyzers/SoftdropAnalyzer`, do:

```
cmsRun python/MakeParticleHierarchy.py outputFile=output.root inputFiles=file:/path/to/input/genfile.root
```
