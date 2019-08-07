#!/bin/env python

#Script to make job scripts to run multiple jobs in condor.
#Make toys, then do a fit on them.
import os
import sys
def main():
#Default values
        #How many toys in each combine call
        ntoys=10
        #How many times to call combine in a single job
        ncalls=10
        #How many jobs to run
        njobs = 500
        #function to use to generate toys
        gen_func = "Pol5"
        #function to use to fit toys
        fit_func = "Pol5"
        #What category number to start at
        catstart = 0
        #what category number to end at
        catend = 8
        #parameters for function to make a random seed
        seed_params = [123456, 1738, 679420, 421738, 42]
        #job flavour (how long to let it run).
        flavour = "longlunch"
#        flavour = "workday"  
#parse input arguments to see if any default values are being replaced
        #only -h or --help is a valid option without a value after it
        if len(sys.argv) % 2 != 1 and len(sys.argv) != 2:
                sys.exit("Error! Please specify arguments as --[option] [value], eg --ntoys 15")
	i = 1
	while i < len(sys.argv)-1 or len(sys.argv) == 2:
                if sys.argv[i] == "--ntoys" or sys.argv[i] == "-t":
                        ntoys = int(sys.argv[i+1])
                elif sys.argv[i] == "--ncalls" or sys.argv[i] == "-n":
                        ncalls = int(sys.argv[i+1])
                elif sys.argv[i] == "--njobs" or sys.argv[i] == "-j":
                        njobs = int(sys.argv[i+1])
                elif sys.argv[i] == "--gen_func" or sys.argv[i] == "-g":
                        gen_func = sys.argv[i+1]
                elif sys.argv[i] == "--fit_func" or sys.argv[i] == "-f":
                        fit_func = sys.argv[i+1]
                elif sys.argv[i] == "--CATS" or sys.argv[i] == "--CATs" or sys.argv[i] == "--cats" or sys.argv[i] == "-c":
                        cats = sys.argv[i+1].split(",")
                        if len(cats) != 2:
                                sys.exit("Error! Please use the --CATS option as: --CATS catstart,catend, eg --CATS 0,0 for only category 0.")
                        catstart = int(cats[0])
                        catend = int(cats[1])
                elif sys.argv[i] == "--seed_params" or sys.argv[i] == "-s":
                        params = sys.argv[i+1].split(",")
                        if len(params) > 5:
                                sys.exit("Error! Please specify no more than 5 seed parameters.")
                        for j in range(len(params)):
                                seed_params[j] = int(params[j])
                else:
                        if not (sys.argv[i] == "--help" or sys.argv[i] == "-h"):
                                print("Error! Unrecognized option %s." %sys.argv[i])
                        else:
                                print("Script to generate executable files for combine calls, and make directories to copy them to.")
                        print("Available options: --ntoys (-t), --ncalls (-n), --njobs (-j), --gen_func (-g), --fit_func (-f), --CATS (-c), --seed_params (-s)")
                        print("Usage: ./make_scripts.py [option] [value] [option] [value]...")
                        sys.exit(0)
                i += 2
#
        #print selected options
        #(total number of toys will be ntoys*ncalls*njobs
        totntoys = ntoys*ncalls*njobs
        print("Writing executables for %d total toys." %(totntoys))
        print("Generating function = %s, fitting function = %s" %(gen_func, fit_func))
        print("Categories %d - %d" %(catstart, catend))
        print("Random seed parameters: " + str(seed_params))
        #datacard to generate toys based on
        gen_datacard = "/afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett/datacards/datacard_vbfHbb_bias%s_m125_CAT%d-CAT%d_CATveto.txt" %(gen_func, catstart, catend)
        #if the generating datacard doesn't exist yet, make it.
        if not os.path.exists(gen_datacard):
                os.system("eval `scramv1 runtime -sh`; /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/mkDatacards_run2_cat.py --CATS %d,%d --bias --function %s --TF ConstPOL1,ConstPOL1 --workdir /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett"%(catstart, catend, gen_func))
        #datacard to fit toys based on
        fit_datacard = "/afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett/datacards/datacard_vbfHbb_bias%s_m125_CAT%d-CAT%d_CATveto.txt" %(fit_func, catstart, catend)
        #make sure the model exists also.
        if not os.path.exists("/afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett/root/bias_shapes_workspace_%s.root" %(gen_func)):
                #generate the bias template.
                os.system("eval `scramv1 runtime -sh`; /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/myBiasTemplates.py --function %s --TF ConstPOL1,ConstPOL1 --workdir /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett"%(gen_func))
        if not os.path.exists("/afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett/root/bias_shapes_workspace_%s.root" %(fit_func)):
                #generate the bias template.
                os.system("eval `scramv1 runtime -sh`; /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/myBiasTemplates.py --function %s --TF ConstPOL1,ConstPOL1 --workdir /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett"%(fit_func))
        
        if not os.path.exists(fit_datacard):
                os.system("eval `scramv1 runtime -sh`; /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/mkDatacards_run2_cat.py --CATS %d,%d --bias --function %s --TF ConstPOL1,ConstPOL1 --workdir /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett"%(catstart, catend, fit_func))
        #make directory in eos that the scripts can copy their big (root) files to when they're done.
        new_head = "/eos/user/b/bgreenbe/cat%d_%d/%s_%s" %(catstart, catend, gen_func, fit_func)
        #if the category directory doesn't exist yet, make that first.
        cat_dir = "/eos/user/b/bgreenbe/cat%d_%d"%(catstart, catend)
        if not os.path.exists(cat_dir):
                os.system("mkdir " + cat_dir)
        if os.path.exists(new_head):
                print("WARNING: Directory %s exists. Are you sure you want to overwrite?"%new_head)
        else:
                os.system("mkdir " + new_head )
        my_directory = os.getcwd()
        for job in range(njobs):
                filename = "cat%d_%d_%s_%s_job%d.sh" %(catstart, catend, gen_func, fit_func, job)
                jobfile = open(filename, "w")
                jobfile.write("#!/bin/bash\n")
                #change to the correct directory
                jobfile.write("cd " + my_directory + "\n") 
                #cmsenv
                jobfile.write("eval `scramv1 runtime -sh`\n")
                #random seed for fitting
                #Do 10 different combine calls, each with 10 toys and a different random seed.
                for c in range(ncalls):
                        seed = job ** job % seed_params[0] - seed_params[1]*(c*c) + job*c + seed_params[2]*c + seed_params[3] 
                        if seed > 2147483647 or seed < -2147483648:
                                sys.exit("Error: random seed out of range: %d. Please specify smaller (or larger) parameters so that the seed is between -2147483648 and 2147483647."%seed)
                #generate toys
                        jobfile.write("combine -M GenerateOnly --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --X-rtd ADDNLL_RECURSIVE=0 --X-rtd FITTER_NEW_CROSSING_ALGO " +  \
                        "--toysFrequentist -t " + str(ntoys) + " --expectSignal 1.0 --saveToys  --rMin=-100 --rMax=100 -n _" + str(job) + " --seed=" + str(seed) +   \
                        " " + gen_datacard + "\n") 
                        toyname = "higgsCombine_" + str(job) + ".GenerateOnly.mH120." + str(seed) + ".root"
                        #run fit
                        seed = seed + seed_params[4]
                        if seed > 2147483647 or seed < -2147483648:
                                sys.exit("Error: random seed out of range: %d. Please specify smaller (or larger) parameters so that the seed is between -2147483648 and 2147483647."%seed)
                        jobfile.write("combine -M FitDiagnostics --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --X-rtd ADDNLL_RECURSIVE=0" + \
                        " --X-rtd FITTER_NEW_CROSSING_ALGO -t " + str(ntoys)  + " --seed=" + str(seed) \
                        + " --expectSignal=1.0 --robustFit=1  --rMin=-100 --rMax=100 --toysFile " + str(toyname) + " -n _" + str(job) + "_" + str(c) +  \
                        " --setRobustFitTolerance 0.01 -d " + fit_datacard + "\n")
                #Now move everything to eos (so afs storage space doesn't run out).
                newdir = new_head + "/job" + str(job)
                jobfile.write("mkdir " + newdir + "\n") 
                jobfile.write("mv fitDiagnostics_" + str(job) + "_* " + newdir + "\n")
                jobfile.write("mv higgsCombine_" + str(job) + "_* " + newdir + "\n")
                jobfile.write("mv higgsCombine_" + str(job) + ".* " + newdir)
        #now make .sub file for condor submission.
        subname = "cat%d_%d_%s_%s.sub"%(catstart,catend,gen_func,fit_func)
        subfile = open(subname, "w")
        subfile.write("executable            = cat%d_%d_%s_%s_job$(ProcId).sh\n" %(catstart, catend, gen_func, fit_func)) 
        subfile.write("arguments             = $(ClusterId) $(ProcId)\n")
        subfile.write("output                = cat%d_%d_%s_%s_$(ProcId).out\n" %(catstart, catend, gen_func, fit_func)) 
        subfile.write("error                 = cat%d_%d_%s_%s_$(ProcId).err\n" %(catstart, catend, gen_func, fit_func)) 
        subfile.write("log                   = cat%d_%d_%s_%s_$(ProcId).log\n" %(catstart, catend, gen_func, fit_func)) 
        subfile.write("+JobFlavour           = \"%s\"\n" %(flavour) )
        subfile.write("queue %d" %(njobs))
        subfile.close()
       #now make script to move the lagging files to the correct directory (if the job gets killed before it can move them)
        mvfile = open("mv_files.sh", "w")
        for job in range(njobs):
                newdir = "/eos/user/b/bgreenbe/cat%d_%d/%s_%s/job%d" %(catstart, catend, gen_func, fit_func, job)
                mvfile.write("mkdir " + newdir + "\n")
                mvfile.write("mv fitDiagnostics_" + str(job) + "* " + newdir + "\n")
                mvfile.write("mv higgsCombine_" + str(job) + ".* " + newdir + "\n")
                mvfile.write("mv higgsCombine_" + str(job) + "_* " + newdir + "\n")
        mvfile.close()
        #make it executable
        os.system("chmod +x mv_files.sh")
        
        print("Done! Submit job with:")
        print("condor_submit %s"%subname)
if __name__ == "__main__":
        main()
