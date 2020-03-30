#!/bin/env python
import ROOT
from ROOT import *
from ROOT import gStyle
import sys
import os

gROOT.SetBatch(True)
#gStyle.SetPadTopMargin(0.06)
#gStyle.SetPadRightMargin(0.04)
#gStyle.SetPadLeftMargin(0.15)
#how many jobs were ran
njobs = 500
#how many calls to combine in each job
ncalls = 10
#How many toys in each call to combine
ntoys = 5
#total number of toys
totntoys = njobs*ncalls*ntoys
print("%d toys total"%(totntoys))
r = [0.0 for x in range(totntoys)]
rErr = [0.0 for x in range(totntoys)]
#textbox min and max (y-value)
tbmin = 50
tbmax = 400
#how many bins in the histogram
nbins = 100
#name of generating model
gen_func = "Pol3"
#prefix of .out files (before number)
fit_func = "Pol3"
#starting category
catstart = 1
#ending category
catend = 1
#Expected mean
mu = "1"
#special name of directory (past catn_n)
spec_name = ""
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
  elif sys.argv[i] == "--dir" or sys.argv[i] == "-d":
    spec_name = sys.argv[i+1]
    if spec_name[0] != "/":
      spec_name = "/" + spec_name
    if spec_name[len(spec_name)-1] != "/":
      spec_name = spec_name + "/"
  elif sys.argv[i] == "--expectSignal" or sys.argv[i] == "-e":
    mu = sys.argv[i+1]
  else:
    if not (sys.argv[i] == "--help" or sys.argv[i] == "-h"):
      print("Error! Unrecognized option %s." %sys.argv[i])
    else:
      print("Script to generate executable files for combine calls, and make directories to copy them to.")
      print("Available options: --ntoys (-t), --ncalls (-n), --njobs (-j), --gen_func (-g), --fit_func (-f), --CATS (-c), --seed_params (-s)")
      print("Usage: ./make_scripts.py [option] [value] [option] [value]...")
      sys.exit(0)
  i += 2
hist_all = ROOT.TH1F("hist_all","",nbins,-5,5)
#hist_all = ROOT.TH1F("hist_all","",nbins,0,20)
#hist_all = ROOT.TH2D("hist_all","",nbins,-40,40, 100, 10, 20)
canv = ROOT.TCanvas("canv","canv",900,900)
canv.cd()
mu_inj_avg = 0.0
nevents = 0
if os.path.exists("r.txt") and os.path.exists("rErr.txt"):
  print("Reading from input file r.txt and rErr.txt")
  rfile = open("r.txt", "r")
  rErrfile = open("rErr.txt", "r")
  for line in rfile:
    r[nevents] = float(line)
    rErrline = rErrfile.readline()
    if not rErrline:
      sys.exit("Error! Not the same number of entries in r.txt and rErr.txt")
    rErr[nevents] = float(rErrline)
    nevents += 1  
else:
  rfile = open("r.txt", "w")
  rErrfile = open("rErr.txt", "w")
  #1000 toys separated into 100 jobs labeled 0-99, 10 combine calls
  for num in range(njobs):
#    outname = "../cat%d_%d_%s_%s_%d.out"%(catstart, catend, gen_func, fit_func, num)
    #check to make sure .out file exists (if not there was a problem so skip this job)
  #  if not os.path.exists(outname):
  #    print("Error! %s (output file) does not exist! Skipping job %d" %(outname, num))
  #    continue
  #  if os.path.exists(outname):
#    outfile = open(outname)
    #10 different combine calls
    #close then reopen the files every so often in case the process gets killed after 2 hours (fml)
    if num % 50 == 49:
        rfile.close()
        rErrfile.close()
        rfile = open("r.txt", "a")
        rErrfile = open("rErr.txt", "a")
    for c in range(ncalls):
  #    filename = "../fitDiagnostics_" + str(num) + "_" + str(c) + ".root"
      filename = "/eos/user/b/bgreenbe/cat%d_%d%s/%s_%s/job%d/fitDiagnostics_%d_%d.root" %(catstart, catend, spec_name, gen_func, fit_func, num, num, c)
      tfile = ROOT.TFile.Open(filename)
      if not tfile:
        print("Error! Could not read file %s. Skipping job %d call %d" %(filename, num, c))
        continue
      tree = tfile.Get("tree_fit_sb")
      print("scanning tree for job %d call %d"%(num, c))
      if not tree:
        print("Error! Could not read fit in tree for job %d call %d" %(num, c))
        continue
      #Fill histogram with value (r-2)/rErr for each of the 10 toys.
  #    r[num] = 0
  #    rErr[num] = 0
      #goto r branch 
      # tree.SetBranchAddress("r", r)
      # tree.SetBranchAddress("rErr", rErr)
     #end of file reached?
#      eof = False
     #Read beginning of outfile until toy 0 starts, at --FitDiagnostics--.
#      out_text = ""
#      while " --- FitDiagnostics ---" not in out_text: 
#        out_text = outfile.readline()
#        if not out_text:
    #      print("Error: reached end of output file %s. Skipping calls %d-9 in job %d" %(outname, c, num))
#  	  eof = True
#          break
     # if eof: break
    #using try/except to catch the stupid 'new version of root' error.
      try:
          for i, event in enumerate(tree):
            if i > ntoys: 
              print("Error! more than %d events"%(ntoys))
              break
            #First make sure this is a valid toy.
    #        use_event = True
      #     # print("Reading toy %d" %i)
    #        if not eof: out_text = outfile.readline() #next line after --FitDiagnostics--
      #   # Do preliminary check for "Fit failed."
      #      if "Fit failed." in out_text:
      #        use_event = False
      #        print("Error: Fit failed. Not using job %d toy %d" %(num, i))
      #   # Now out_text is not --FitDiagnostics-- anymore.
      #   # When --FitDiagnostics-- is found again, that's the start of the next event.
      #   # Also need to account for eof: "mean" will be start of one of the last lines in the file--after every call!
    #        while " --- FitDiagnostics ---" not in out_text and out_text != "" and not eof: # and "mean" not in out_text: 
    #          out_text = outfile.readline()
      #       # print(out_text)
      #   #  If there's an error mssg or warning with this toy, discard it.
    #         ONLY if it says Error: Fit failed. or Error: Search did not converge. discard the toy.
    #          if "Error: search did not converge" in out_text or "Fit failed." in out_text:
    #            use_event = False
      #        if "Error" in out_text: 
      #        #only print error message again if it hasn't already been printed.
      #          if use_event: 
      #            print("Error! Not using job %d call %d toy %d" %(num, c, i))
      #          use_event = False
      #        elif "[WARNING]" in out_text:
      #          if use_event:
      #            print("WARNING! Not using job %d call %d toy %d" %(num, c, i))
      #          use_event = False
      #        elif "Fit failed." in out_text:
      #          if use_event:
      #            print("Error: Fit failed. Not using job %d call %d toy %d" %(num, c, i))
      #          use_event = False
      #      #If there was any type of error or warning, goto next event.
    #        if not use_event:
    #          print("Search did not converge. Skipping job %d call %d toy %d."%(num, c, i))
    #          continue
            #get r value of this toy
            #If it gets past here, the event is valid!
            r[nevents] = event.r
            rErr[nevents] = event.rErr
         #  if rErr==0, throw away event
            if rErr[nevents] == 0:
              print("Error: rErr = 0, r = %f, num = %d, i = %d"%(r[nevents], num, i))
      #        rErr = 0.0001
              continue
           #Ok at this point, the fit is successful and bias is finite, so toy is valid.
           #write to output file to save time later.
            rfile.write(str(r[nevents]) + "\n")
            rErrfile.write(str(rErr[nevents]) + "\n")
            nevents += 1
      except:
          print("Error: new version of root or something stupid like that.")
rfile.close()
rErrfile.close()

if nevents == 0:
  sys.exit("Error! 0 valid events!")
#how many events are actually un-cut
real_nevents = 0
print("nevents = " + str(nevents))
#make dict so that each value only occurs once
pulls = dict()
#now the arrays are all filled so we can go through them and fill the histogram.
for i in range(nevents):
   #only accept toys with rErr > 11.5
  pull =  (r[i]-float(mu)) / rErr[i]
#  if abs( pull -  -2.20987595716) < 0.00001:
#  if not rErr[i] > 12:
#    print("Not using event %d out of %d: pull = %f" %(i, nevents, pull))
#    continue
  if str(pull) in pulls: # and pulls[str(pull)] > 5:
    print("Not using event %d out of %d: pull = %f" %(i, nevents, pull))
    continue
#  if not str(pull) in pulls:
  pulls[str(pull)] = 0
#  pulls[str(pull)] += 1
  if pull == 0.0:
    print("Not using event %d out of %d: pull = %f" %(i, nevents, pull))
    continue
  real_nevents += 1
  mu_inj_avg += r[i]
  hist_all.Fill((r[i]-float(mu))/rErr[i])
#  hist_all.Fill(rErr[i])
 #hist_all.Fill(r[i])
 #hist_all.Fill(r[i], rErr[i])
   #close outfile once the job is done.
#  outfile.close() #close at end of job (10 calls)
if real_nevents == 0:
  sys.exit("Error! 0 valid (uncut) events!")
print("real nevents = " + str(real_nevents))
mu_inj_avg /= real_nevents
 #tree.Print()
# tree.Draw("(r-%s)/rErr"%mu)
# tree.Draw("(r-%s)/rErr>>hist_all"%mu)
#canv.Draw()
#mean = hist_all.GetMean()
#rms =  hist_all.GetRMS()
hist_all.SetStats(True)
#hist_all.SetLineColor(1)
ymax=hist_all.GetMaximum()*1.1
hist_all.GetYaxis().SetRangeUser(0.,ymax)
#hist_all.Fit("gaus")
#hist_final.GetFunction("gaus").SetParLimits(2,0,1);
#hist_all.GetFunction("gaus").SetLineWidth(3);

hist_all.GetXaxis().SetTitle("(#mu-#mu_{inj})/#sigma")
#hist_all.GetXaxis().SetTitle("#sigma")
#hist_all.GetXaxis().SetTitle("(#mu)")
hist_all.GetYaxis().SetTitle("Toys")
#hist_all.GetYaxis().SetTitle("rErr")
pave22 = ROOT.TPaveText(3.50,tbmin,9.0,tbmax);
pave22.AddText("#mu_{inj}=%s"%mu);
pave22.AddText("Generating model = %s" %gen_func)
pave22.AddText("Fitting model = %s" %fit_func)
pave22.AddText("Mean r for toys = %f" %mu_inj_avg)
pave22.SetFillColor(0)
pave22.SetTextFont(42)
#pave22.SetTextColor(ROOT.kBlue)
pave22.SetTextSize(0.03)
#pave22.SetBorderSize(0)
#gPad.Update()
#
#
#pave2 = ROOT.TPaveText(0.2,0.65,0.4,0.9,"NDC");
#pave2.AddText(toy_func);
#pave2.AddText("Fitted with %s"%fit_func);
#pave2.AddText("#mu_{inj}=%s"%mu);
#pave2.AddText("mean=%0.2f"%mean);
#pave2.AddText("rms=%0.2f"%rms);
#pave2.AddText("# toys =%0.2f"%n_toys[num]);
#pave2.SetFillColor(0)
#pave2.SetTextFont(42)
#pave2.SetTextColor(ROOT.kBlue)
#pave2.SetTextSize(0.04)
#pave2.SetBorderSize(0)
#pave2.Draw();
#gPad.Update()
#right,top   = gStyle.GetPadRightMargin(),gStyle.GetPadTopMargin()
#pave2.SetY1NDC(pave2.GetY2NDC()-top*0.9*pave2.GetListOfLines().GetSize())
#pave2.Draw();
#archive+=[pave2]
#archive+=[pave22]
hist_all.Draw()
#pave22.Draw();
print("Average mu_inj = %f" %mu_inj_avg)
#make display not immediately close
#s = raw_input("Hit enter to exit.")
canv.SaveAs("cat%d_%s_%s_%s.pdf"%(catstart, gen_func, fit_func, spec_name[1:len(spec_name)-1]))
#don't use the first and last chars of spec_name: these are '/'.
#canv.SaveAs("plots/plot_inj%s_%s_cat%s.pdf"%(mu,fit_funcs[0],cat
#canv.SaveAs("plot_inj%s_%s_cat%s.pdf"%(mu,fit_funcs[0],cat))
