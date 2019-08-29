#!/bin/bash
#script to make FTest tables for all 9 categories. Each will be stored in FTable_x.txt (x=0..8)
cd /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src
#cmsenv
eval `scramv1 runtime -sh`
#./mkTransferFunctions_run2_cat.py --workdir=test2_for_bennett --TF ConstPOL1,ConstPOL1
#./mkSigTemplates_run2_cat.py --workdir=test2_for_bennett
#./mkBkgTemplates_run2_cat.py --workdir=test2_for_bennett  
#./mkDataTemplates_run2_cat.py --workdir=test2_for_bennett --TF ConstPOL1,ConstPOL1
for i in $(seq 1 8)
do
  echo "Running expPol$i"
  ./myBiasTemplates.py --workdir test_for_bennett --function expPol$i --TF ConstPOL1,ConstPOL1 
  #./myBiasTemplates.py --workdir test2_for_bennett --function expx$i --TF ConstPOL1,ConstPOL1 
#  ./myBiasTemplates.py --workdir test2_for_bennett --function sine$i --TF ConstPOL1,ConstPOL1 
done
python FTest.py
