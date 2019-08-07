#!/bin/env python
import ROOT
from ROOT import *
from ROOT import gStyle
import sys
import os

#gROOT.SetBatch(True)
#gStyle.SetPadTopMargin(0.06)
#gStyle.SetPadRightMargin(0.04)
#gStyle.SetPadLeftMargin(0.15)
#how many jobs were ran
njobs = 100
#how many calls to combine in each job
ncalls = 10
#How many toys in each call to combine
ntoys = 10
#textbox min and max (y-value)
tbmin = 50
tbmax = 200
#how many bins in the histogram
nbins = 80
#name of generating model
gen_func = "Pol5"
#prefix of .out files (before number)
fit_func = "Pol5"
#starting category
catstart = 0
#ending category
catend = 0
#Expected mean
mu = "1"
hist_all = ROOT.TH1F("hist_all","",nbins,-10,10)
#hist_all = ROOT.TH2D("hist_all","",nbins,-40,40, 100, 10, 20)
canv = ROOT.TCanvas("canv","canv",900,900)
canv.cd()
mu_inj_avg = 0.0
nevents = 0
#1000 toys separated into 100 jobs labeled 0-99, 10 combine calls
for num in range(njobs):
  outname = "../cat%d_%d_%s_%s_%d.out"%(catstart, catend, gen_func, fit_func, num)
  #check to make sure .out file exists (if not there was a problem so skip this job)
  if not os.path.exists(outname):
    print("Error! %s (output file) does not exist! Skipping job %d" %(outname, num))
    continue
  outfile = open(outname)
  #10 different combine calls
  for c in range(ncalls):
#    filename = "../fitDiagnostics_" + str(num) + "_" + str(c) + ".root"
    filename = "/eos/user/b/bgreenbe/cat%d_%d/%s_%s/job%d/fitDiagnostics_%d_%d.root" %(catstart, catend, gen_func, fit_func, num, num, c)
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
    r = 0
    rErr = 0
    #goto r branch 
    # tree.SetBranchAddress("r", r)
    # tree.SetBranchAddress("rErr", rErr)
   #end of file reached?
    eof = False
   #Read beginning of outfile until toy 0 starts, at --FitDiagnostics--.
    out_text = ""
    while " --- FitDiagnostics ---" not in out_text: 
      out_text = outfile.readline()
      if not out_text:
        print("Error: reached end of output file %s. Skipping calls %d-9 in job %d" %(outname, c, num))
	eof = True
        break
    if eof: break
    for i, event in enumerate(tree):
      if i > ntoys: 
        print("Error! more than %d events"%(ntoys))
        break
      #First make sure this is a valid toy.
      use_event = True
     # print("Reading toy %d" %i)
      out_text = outfile.readline() #next line after --FitDiagnostics--
   # Do preliminary check for "Fit failed."
      if "Fit failed." in out_text:
        use_event = False
        print("Error: Fit failed. Not using job %d toy %d" %(num, i))
   # Now out_text is not --FitDiagnostics-- anymore.
   # When --FitDiagnostics-- is found again, that's the start of the next event.
   # Also need to account for eof: "mean" will be start of one of the last lines in the file--after every call!
      while " --- FitDiagnostics ---" not in out_text and out_text != "": # and "mean" not in out_text: 
        out_text = outfile.readline()
       # print(out_text)
   #  If there's an error mssg or warning with this toy, discard it.
        if "Error" in out_text: 
        #only print error message again if it hasn't already been printed.
          if use_event: 
            print("Error! Not using job %d call %d toy %d" %(num, c, i))
          use_event = False
        elif "[WARNING]" in out_text:
          if use_event:
            print("WARNING! Not using job %d call %d toy %d" %(num, c, i))
          use_event = False
        elif "Fit failed." in out_text:
          if use_event:
            print("Error: Fit failed. Not using job %d call %d toy %d" %(num, c, i))
          use_event = False
      #If there was any type of error or warning, goto next event.
      if not use_event:
        continue
      #get r value of this toy
      #If it gets past here, the event is valid!
      r = event.r
      rErr = event.rErr
   #  if rErr==0, throw away event
      if rErr == 0:
        print("Error: rErr = 0, r = %f, num = %d, i = %d"%(r, num, i))
        rErr = 0.0001
        continue
     #only accept toys with rErr > 10
      if not rErr > 10:
        continue
     #Ok at this point, the fit is successful and bias is finite, so toy is valid.
      nevents += 1
      mu_inj_avg += r
      hist_all.Fill((r-1.0)/rErr)
      #hist_all.Fill(r)
      #hist_all.Fill(r, rErr)
   #close outfile once the job is done.
  outfile.close() #close at end of job (10 calls)
if nevents == 0:
  sys.exit("Error! 0 valid events!")
mu_inj_avg /= nevents
 #tree.Print()
# tree.Draw("(r-%s)/rErr"%mu)
# tree.Draw("(r-%s)/rErr>>hist_all"%mu)
#canv.Draw()
#mean = hist_all.GetMean()
#rms =  hist_all.GetRMS()
hist_all.SetStats(True)
#hist_all.SetLineColor(1)
#ymax=hist_all.GetMaximum()*1.1
#hi"/afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/workdir/datacards/datacard_CAT0_Pol6.root"st_all.GetYaxis().SetRangeUser(0.,ymax)
#hist_all.Fit("gaus")
#hist_final.GetFunction("gaus").SetParLimits(2,0,1);
#hist_all.GetFunction("gaus").SetLineWidth(3);

hist_all.GetXaxis().SetTitle("(#mu-#mu_{inj})/#sigma")
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
s = raw_input("Hit enter to exit.")
#c.SaveAs("plot_inj1_%s_%s.pdf"%(toy_func,fit_func))
#canv.SaveAs("plots/plot_inj%s_%s_cat%s.pdf"%(mu,fit_funcs[0],cat
#canv.SaveAs("plot_inj%s_%s_cat%s.pdf"%(mu,fit_funcs[0],cat))
