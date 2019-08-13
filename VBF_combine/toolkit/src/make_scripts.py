#!/bin/env python

#Script to make job scripts to run multiple jobs in condor.
#Make toys, then do a fit on them.
import os
import sys
def main():
#Default values
        #How many toys in each combine call
        ntoys=5
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
        flavour = "microcentury" #"longlunch"
        #verbosity for combine calls
        verbosity = 0
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
                elif sys.argv[i] == "--verbosity" or sys.argv[i] == "-v":
                        verbosity = int(sys.argv[i+1])
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
        #gencard is just the name without the path.
        gencard = "datacards/datacard_vbfHbb_bias%s_m125_CAT%d-CAT%d_CATveto.txt" %(gen_func, catstart, catend)
        #if the generating datacard doesn't exist yet, make it.
        if not os.path.exists(gen_datacard):
                os.system("eval `scramv1 runtime -sh`; /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/mkDatacards_run2_cat.py --CATS %d,%d --bias --function %s --TF ConstPOL1,ConstPOL1 --workdir /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett"%(catstart, catend, gen_func))
        #datacard to fit toys based on
        fit_datacard = "/afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett/datacards/datacard_vbfHbb_bias%s_m125_CAT%d-CAT%d_CATveto.txt" %(fit_func, catstart, catend)
        #fitcard is just the name without the path.
        fitcard = "datacards/datacard_vbfHbb_bias%s_m125_CAT%d-CAT%d_CATveto.txt" %(fit_func, catstart, catend)
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
      #make sure random seeds aren't repeated!
        used_seeds = dict()
        for job in range(njobs):
                filename = "cat%d_%d_%s_%s_job%d.sh" %(catstart, catend, gen_func, fit_func, job)
                jobfile = open(filename, "w")
                jobfile.write("#!/bin/bash\n")
                jobfile.write("echo \"Starting in $PWD\"\n")
                #change to the correct directory JUST to cmsenv, then switch back. 
                jobfile.write("cd " + my_directory + "\n") 
                #cmsenv
                jobfile.write("eval `scramv1 runtime -sh`\n")
                #come back
                jobfile.write("cd -\n")
                jobfile.write("echo \"now in $PWD\"\n")
                #copy datacards to here.
                jobfile.write("mkdir datacards\n")
                jobfile.write("cp %s datacards\n" %(gen_datacard))
                #don't need to copy the same datacard twice.
                if not fit_func == gen_func:
                        jobfile.write("cp %s datacards\n"%(fit_datacard))
                #now copy the necessary templates.
                jobfile.write("mkdir root\n")
                jobfile.write("cp /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett/root/data_shapes_workspace.root root\n")
                jobfile.write("cp /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett/root/sig_shapes_workspace.root root\n")
                jobfile.write("cp /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett/root/bias_shapes_workspace_%s.root root\n"%(gen_func))
                if not gen_func == fit_func:
                        jobfile.write("cp /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett/root/bias_shapes_workspace_%s.root root\n"%(fit_func))
                newdir = new_head + "/job" + str(job)
                jobfile.write("mkdir " + newdir + "\n") 
                #random seed for fitting
                #Do 10 different combine calls, each with 10 toys and a different random seed.
                for c in range(ncalls):
                        seed = job ** job % seed_params[0] - seed_params[1]*(c*c) + job*c + seed_params[2]*c + seed_params[3] 
                        while seed > 2147483647: 
                                seed /= 2
                        while seed < -2147483648:
                                seed += 20000000
                        while str(seed) in used_seeds:
                                seed += 1
                        used_seeds[str(seed)] = 1
                #generate toys
                        jobfile.write("combine -M GenerateOnly --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --X-rtd ADDNLL_RECURSIVE=0 --X-rtd FITTER_NEW_CROSSING_ALGO " +  \
                        "--toysFrequentist -t " + str(ntoys) + " --expectSignal 1.0 --saveToys  --rMin=-100 --rMax=100 -n _" + str(job) + " --seed=" + str(seed) +   \
                        " " + gencard + "\n") 
                        toyname = "higgsCombine_" + str(job) + ".GenerateOnly.mH120." + str(seed) + ".root"
                        #run fit
                        seedf = seed + seed_params[4]
                        while seedf > 2147483647:
                                seedf /= 2
                        while seedf < -2147483648:
                                seed += 20000000
                        while str(seedf) in used_seeds:
                                seedf += 1
                        used_seeds[str(seedf)] = 1
#                                sys.exit("Error: random seed out of range: %d. Please specify smaller (or larger) parameters so that the seed is between -2147483648 and 2147483647."%seed)
                        jobfile.write("combine -M FitDiagnostics --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --X-rtd ADDNLL_RECURSIVE=0" + \
                        " --X-rtd FITTER_NEW_CROSSING_ALGO -t " + str(ntoys)  + " --seed=" + str(seedf) + " -v " + str(verbosity) \
                        + " --expectSignal=1.0 --robustFit=1  --rMin=-100 --rMax=100 --toysFile " + str(toyname) + " -n _" + str(job) + "_" + str(c) +  \
                        " --setRobustFitTolerance 0.01 -d " + fitcard + "\n")
                       #mv the files from this call to eos.
                        jobfile.write("mv higgsCombine_%d.GenerateOnly.mH120.%d.root %s\n"%(job, seed, newdir))
                        jobfile.write("mv higgsCombine_%d_%d.FitDiagnostics.mH120.%d.root %s\n"%(job, c, seedf, newdir))
                        jobfile.write("mv fitDiagnostics_%d_%d.root %s\n"%(job, seedf, newdir))
                #Now move everything to eos (so afs storage space doesn't run out).
#                jobfile.write("mv fitDiagnostics_" + str(job) + "_* " + newdir + "\n")
#                jobfile.write("mv higgsCombine_" + str(job) + "_* " + newdir + "\n")
#                jobfile.write("mv higgsCombine_" + str(job) + ".* " + newdir)
        #don't need those copies of datacards and templates anymore.
        jobfile.write("rm -rf datacards\n")
        jobfile.write("rm -rf root\n")
        jobfile.close()
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
#       #now make script to move the lagging files to the correct directory (if the job gets killed before it can move them)
        mvfile = open("mv_files.sh", "w")
        for job in range(njobs):
                newdir = "/eos/user/b/bgreenbe/cat%d_%d/%s_%s/job%d" %(catstart, catend, gen_func, fit_func, job)
                mvfile.write("mkdir " + newdir + "\n")
                mvfile.write("mv fitDiagnostics_" + str(job) + "_* " + newdir + "\n")
                mvfile.write("mv higgsCombine_" + str(job) + ".* " + newdir + "\n")
                mvfile.write("mv higgsCombine_" + str(job) + "_* " + newdir + "\n")
        #now all that remains are the error .root files. Move these to the head.
        mvfile.write("mv *.root " + new_head + "\n")
        mvfile.write("mv *.dot " + new_head + "\n")
        #rm any remaining roostats files (they're useless anyway)
        mvfile.write("rm roostats*\n")
        #only keep one copy of the executable files (to save space)
        mvfile.write("mv cat%d_%d_%s_%s_job0.sh temp\n"%(catstart, catend, gen_func, fit_func))
        mvfile.write("rm cat*.sh\n")
        mvfile.write("mv temp cat%d_%d_%s_%s_job0.sh\n"%(catstart, catend, gen_func, fit_func))
        #make bias plot
        mvfile.write("mkdir plots\n")
        mvfile.write("cp ~/public/plot_bias.py plots\n")
        mvfile.write("cd plots\n")
        mvfile.write("./plot_bias.py --CATS %d,%d --gen_func %s --fit_func %s --njobs %d --ntoys %d --ncalls %d\n"%(catstart, catend, gen_func, fit_func, njobs, ntoys, ncalls))
        mvfile.close()
        #make it executable
        os.system("chmod +x mv_files.sh")
        
        print("Done! Submit job with:")
        print("condor_submit %s"%subname)
if __name__ == "__main__":
        main()
