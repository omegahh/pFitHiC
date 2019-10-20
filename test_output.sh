#!/bin/bash -ex

set -o nounset
set -o pipefail
set -o errexit 

# git clone https://github.com/ay-lab/fithic.git
# git checkout 2.0.6
# python setup.py install

DATADIR=test
inI=$DATADIR/contactCounts
inF=$DATADIR/fragmentLists

#####  Settings for only chromosome 1 of        #####
#####  human and mouse embryonic stem cell data     #####
#####  from Dixon et al.                #####

distUpThres=5000000
distLowThres=50000
## other parameters described in fit-hi-c.py
mappabilityThres=1
noOfPasses=1
noOfBins=50
resolution=40000

filename="Dixon_hESC_HindIII_hg18_w40000_chr1.gz"

time fithic -f $inF/$filename -i $inI/$filename -o test/fithic_outputs/ -r $resolution -L $distLowThres -U $distUpThres -b $noOfBins -p $noOfPasses -x intraOnly -v

time python runner.py -f $inF/$filename -i $inI/$filename -o test/pfithic_output -r $resolution -L $distLowThres -U $distUpThres -b $noOfBins -p $noOfPasses -x intraOnly

exit