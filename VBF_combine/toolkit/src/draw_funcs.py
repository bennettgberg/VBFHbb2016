#!/usr/bin/env python

import os,sys,re,json,datetime
from glob import glob
from array import array
from optparse import OptionParser,OptionGroup
import warnings

from  pdf_param_cfi import *
from myGenPdf import *
from generateFormula import *
from signal import signal, SIGINT

warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='.*class stack<RooAbsArg\*,deque<RooAbsArg\*> >' )
basepath=os.path.split(os.path.abspath(__file__))[0]+"/../"
sys.path.append(basepath+"common/")

tempargv = sys.argv[:]
sys.argv = []
import ROOT
from ROOT import *
sys.argv = tempargv

from toolkit import *

colours = ["\033[m"] + ["\033[%d;%dm"%(x,y) for (x,y) in [(0,z) for z in range(31,37)]+[(1,z) for z in range(31,37)]]

gcs={} # for memory issues
gcs_aux={} # for memory issues

workdir = "test_for_bennett"
xmin = 80
xmax = 200
nbins = 120
dX = (xmax - xmin) / nbins
MASS = 125 #higgs mass
TF = "ConstPOL1"
#boundaries between the categories ([min, max] bdt val for each cat)
boundaries = [ [-1, 0.0875],         #cat0
               [0.0875, 0.5775],     #cat1
               [0.5775, 0.7875],     #cat2
               [0.7875, 1.0],        #cat3
               [-1.0, 0.0775],       #cat4
               [0.0775, 0.5775],     #cat5
               [0.5775, 0.7775],     #cat6
               [0.7775, 0.8775],     #cat7
               [0.8775, 1.0]      ]   #cat8

#style options for the canvas
def style():
    gROOT.SetBatch(1)
    gROOT.ProcessLineSync(".x %s/styleCMSTDR.C"%basepath)
    gROOT.ForceStyle()
    gStyle.SetPadTopMargin(0.06)
    gStyle.SetPadRightMargin(0.04)
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetCanvasColor(0)
    gStyle.SetPadColor(0)
    gStyle.SetPadBorderMode(0)
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetFrameBorderMode(0)
    gROOT.ProcessLineSync("gErrorIgnoreLevel = kWarning;")
    RooMsgService.instance().setSilentMode(kTRUE)
    for i in range(2): RooMsgService.instance().setStreamStatus(i,kFALSE)

#note: didn't delete the big blocks of commented code because with new method of getting Pol4/6 fits,
#       have been unable to get error with it (well I get the error but it's nonsensical)
def do_fits(my_cat=0, gen_func="Pol6", fit_func="Pol4", TF="ConstPol1"):
#    style()
    #fBKG = TFile.Open("%s/root/bkg_shapes_workspace.root"%(workdir),"read")
    fBKG = TFile.Open("%s/root/data_shapes_workspace.root"%(workdir),"read")
    fqcd4 = TFile.Open("%s/root/bias_shapes_workspace_%s.root"%(workdir, fit_func), "read") #bpgballin
    #fqcd4 = TFile.Open("%s/plot/draw_funcs/bias_shapes_workspace_Pol4.root"%(workdir), "read")
    fqcd6 = TFile.Open("%s/root/bias_shapes_workspace_%s.root"%(workdir, gen_func), "read")
    wBKG = fBKG.Get("w")
    wqcd4 = fqcd4.Get("w")
    wqcd6 = fqcd6.Get("w")
    w    = RooWorkspace("w","workspace")

#    wqcd4.Print()
#    wqcd6.Print()

## CMS info
   # left,right,top,bottom = gStyle.GetPadLeftMargin(),gStyle.GetPadRightMargin(),gStyle.GetPadTopMargin(),gStyle.GetPadBottomMargin()
   # pCMS1 = TPaveText(left,1.-top,0.4,1.,"NDC")
   # pCMS1.SetTextFont(62)
   # pCMS1.SetTextSize(top*0.75)
   # pCMS1.SetTextAlign(12)
   # pCMS1.SetFillStyle(-1)
   # pCMS1.SetBorderSize(0)
   # pCMS1.AddText("CMS")
   # pCMS2 = TPaveText(0.5,1.-top,1.-right*0.5,1.,"NDC")
   # pCMS2.SetTextFont(42)
   # pCMS2.SetTextSize(top*0.75)
   # pCMS2.SetTextAlign(32)
   # pCMS2.SetFillStyle(-1)
   # pCMS2.SetBorderSize(0)
   # pCMS2.AddText("L = 35.9 fb^{-1} (13 TeV)")
    c0 = TCanvas("can","can",600,600)

    #doubleB is cats 0-3, singleB is cats 4-8
    tag = "double"
## Load tree
    if my_cat < 4:
        fin = TFile.Open("/afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/VBFHbb2016/inputs/FitVBF_BTagCSV_analysis_double_trignone_v25_VBF_newReg.root","read")
    else:
        fin = TFile.Open("/afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/VBFHbb2016/inputs/FitVBF_BTagCSV_analysis_single_trignone_v25_VBF_newReg.root","read")
        tag = "single"
        
    T = fin.Get("VBF/events")
## Containers
    params   = {}
    qcd_pdf      = {}
    model        = {}
    x_name = "mbbReg_CAT%d"%my_cat
    x = RooRealVar(x_name, x_name, xmin, xmax)
    res = {}
    yQ     =  {}
    #do fit separately for orders 3-7
## Category loop (not a loop anymore)
    C = my_cat
#### Start of RooFit part 
### QCD part
### Containers
    catstr = "CAT%d"%my_cat
    h   = TH1F("hMbb_%s"%catstr,"hMbb_%s"%catstr,nbins,xmin,xmax)
    hb  = TH1F("hMbb_blind_%s"%catstr,"hMbb_blind_%s"%catstr,nbins,xmin,xmax)
### Define cut
    cut    = TCut("(bdt_VBF>%1.4f && bdt_VBF<=%1.4f)"%(boundaries[C][0],boundaries[C][1]))
    cutB   = TCut("(bdt_VBF>%1.4f && bdt_VBF<=%1.4f) && mbbRegFSR>100 && mbbRegFSR<150"%(boundaries[C][0],boundaries[C][1]))
### Draw
    c0.cd()
                #c0.SetLogy()
    T.Draw("mbbRegFSR>>%s"%(h.GetName()),cut)
    c0.SaveAs(catstr + "test.png")
    T.Draw("mbbRegFSR>>%s"%(hb.GetName()),cutB)
    c0.SaveAs(catstr + "Btest.png")
### Yields
    #Y = RooRealVar("yield_data_CAT%d"%C,"yield_data_CAT%d"%C,h.Integral())
### Histograms
    rh  = RooDataHist("data_hist_CAT%d"%C,"data_hist_CAT%d"%C,RooArgList(x),h)
    rhb = RooDataHist("data_hist_blind_CAT%d"%C,"data_hist_blind_CAT%d"%C,RooArgList(x),hb)
  ### Yields
    yZ = wBKG.var("yield_ZJets_CAT%d"%C)
    yT = wBKG.var("yield_Top_CAT%d"%C)
    #yq_est = h.Integral()
    yZ.setConstant(kTRUE)
    yT.setConstant(kTRUE)
    yQ4=wqcd4.var("yield_QCD_CAT%d"%(C))
    yQ6=wqcd6.var("yield_QCD_CAT%d"%(C))
  ### PDFs
    zPDF = wBKG.pdf("Z_model_CAT%d"%C)
    tPDF = wBKG.pdf("Top_model_CAT%d"%C)
    #print("literally boutta get the qPDF4 ")
    if TF == "ConstPol1":
        qPDF4= wqcd4.pdf("qcd_model_%s_CAT%d"%(fit_func, C))
        qPDF6= wqcd6.pdf("qcd_model_%s_CAT%d"%(gen_func, C))
    else:
        qPDF4= wqcd4.pdf("qcd_model_%s%s_CAT%d"%(TF, TF, C))
        qPDF6= wqcd6.pdf("qcd_model_%s%s_CAT%d"%(TF, TF, C))
    #print("got 6 too")
    #qPDF4.Print()
    #qPDF6.Print()

### Model
    #for order in range(3, 7):
    #    N = "Pol%d"%order       
    #    n_param = order 
    #    func = "Pol%d"%order
    #    print("starting func %s"%func)
    #    gcs_aux[N] = []
    #    [qcd_pdf[N],params[N]] = generate_pdf(x, pdf_name=func,x_name=x_name,selection=double,gcs=gcs_aux[N],real_cat=my_cat) #, val140=val140) 
 #  #     #print("Just after call to generate_pdf, gcs_aux: %s"%gcs_aux)
    #    if qcd_pdf[N] == None : 
    #        sys.exit("qcd_pdf = None !!!!!!!!!!!!!!!!!")


    #    yQ[N]    = RooRealVar("yield_QCD_CAT%d"%C,"yield_QCD_CAT%d"%C,yq_est,yq_est/2.,2*yq_est)
### #Combined model
    #    yQ[N].setConstant(kTRUE)
    #    #print("before: yZ=%s yT=%s yQ=%s"%(yZ[N], yT[N], yQ[N]))
    #    #print 'before: Yields Z,Top,QCD: ', yZ[N].getVal(),yT[N].getVal(),yQ[N].getVal()
    #    model[N] = RooAddPdf("bkg_model_%s_CAT%d"%(TF,C),"bkg_model_%s_CAT%d"%(TF,C),RooArgList(zPDF,tPDF,qcd_pdf[N]),RooArgList(yZ,yT,yQ[N])) 
### Fit
    #    model[N].Print()
    #    res[N]   = model[N].fitTo(rh,RooFit.Save(),RooFit.Warnings(ROOT.kTRUE))
        #print 'kTrue: Yields Z,Top,QCD: ', yZ[N].getVal(),yT[N].getVal(),yQ[N].getVal()
    #    yQ[N].setConstant(kFALSE)
    #    res[N]   = model[N].fitTo(rh,RooFit.Save(),RooFit.Warnings(ROOT.kTRUE))
        #print("kfalse:yZ=%s yT=%s yQ=%s"%(yZ[N], yT[N], yQ[N]))
        #print 'kFalse:Yields Z,Top,QCD: ', yZ[N].getVal(),yT[N].getVal(),yQ[N].getVal()

    #    chi2 = RooChi2Var("chi2", "chi2", model[N], rh)
    #    chi2_val = chi2.getVal()
        #print 'Yields Z,Top,QCD: ', yZ[N].getVal(),yT[N].getVal(),yQ[N].getVal()
      #  total_fitted_yield = yZ.getVal()+yT.getVal()+yQ[N].getVal()
#        h_data = rh.createHistogram('h_data',x)
        #print "h_data_integral=%.3f"%(h_data.Integral())

#        rh  = RooDataHist("data_hist_CAT%d"%C,"data_hist_CAT%d"%C,RooArgList(x),h.Rebin(10))            
    #print("models: %s"%model)
# part 1
  #  for i in range(3, 7):
  #      poln = "Pol%d"%i
  #      chi2 = RooChi2Var("chi2", "chi2", model[poln], rh)
        #print("AFTER CAT LOOP, CHI2 VALUE FOR %s IS %f"%(poln, chi2.getVal()))
    #need to get error for each of the fits, also plot just the qcd and not the other backgrounds.
  #  model["Pol6"].Print()
    #print("fit result:")
  #  res["Pol6"].Print()

    #histogram of the Pol4 fit
  #  pol4_fitHist = model["Pol4"].createHistogram(x_name, nbins)
    #histogram of the Pol6 fit
  #  pol6_fitHist = model["Pol6"].createHistogram(x_name, nbins)
    #histogram of the Z background data
    Z_hist = zPDF.createHistogram(x_name, nbins)
    #histogram of the top background data
    top_hist = tPDF.createHistogram(x_name, nbins)
    #histogram of the data (just need to convert from RooDataHist to TH1F)
    dat_hist = rh.createHistogram(x_name, nbins)
    #qcd histograms from root files (so not doing the whole ass fit again)
    qcd4_hist = qPDF4.createHistogram(x_name, nbins)
    qcd6_hist = qPDF6.createHistogram(x_name, nbins)

  #  pol4_fitHist.Print()
  #  pol6_fitHist.Print()
    Z_hist.Print()
    top_hist.Print()
    dat_hist.Print()
    x_vals = [0 for i in range(nbins)]
    p4_dat = [0 for i in range(nbins)]
    p4_err = [0 for i in range(nbins)]
    p6_dat = [0 for i in range(nbins)]
    p6_err = [0 for i in range(nbins)]
    datapt = [0 for i in range(nbins)]
    daterr = [0 for i in range(nbins)]
    z_vals = [0 for i in range(nbins)]
    t_vals = [0 for i in range(nbins)]
    q4vals = [0 for i in range(nbins)]
    q6vals = [0 for i in range(nbins)]
    
    print("qcd4 yield: %f"%yQ4.getVal())
    print("qcd6 yield: %f"%yQ6.getVal())
    for i in range(1, nbins+1):
        #x_vals[i-1] = pol4_fitHist.GetBinCenter(i)
        x_vals[i-1] = dat_hist.GetBinCenter(i)
     #   pol4_val = pol4_fitHist.GetBinContent(i)
     #   pol6_val = pol6_fitHist.GetBinContent(i)
        z_vals[i-1] = Z_hist.GetBinContent(i)*yZ.getVal()
        t_vals[i-1] = top_hist.GetBinContent(i)*yT.getVal()
        q4vals[i-1] = qcd4_hist.GetBinContent(i)*yQ4.getVal()
        q6vals[i-1] = qcd6_hist.GetBinContent(i)*yQ6.getVal()
        data_val = dat_hist.GetBinContent(i)
        #subtract Z and top to get actual values for qcd. or just use the actual qcd bro.
        p6_dat[i-1] = q6vals[i-1] #pol6_val #- z_val*yZ.getVal() - top_val*yT.getVal()
        datapt[i-1] = data_val #- z_val*yZ.getVal() - top_val*yT.getVal()
        p4_dat[i-1] = q4vals[i-1] #pol4_val #- z_val*yZ.getVal() - top_val*yT.getVal()
        #also need errors
        daterr[i-1] = dat_hist.GetBinError(i)
      #  p4_err[i-1] = pol4_fitHist.GetBinError(i)
      #  p6_err[i-1] = pol6_fitHist.GetBinError(i)
        p4_err[i-1] = qcd4_hist.GetBinError(i) *p4_dat[i-1]#*yQ4.getVal()
        p6_err[i-1] = qcd6_hist.GetBinError(i) *p6_dat[i-1]#*yQ6.getVal()

    return x_vals, p4_dat, p4_err, p6_dat, p6_err, datapt, daterr, z_vals, t_vals 

#get the signal fit data from the signal template.
# vbf: include the VBF signal or nah? gf: include ggH signal or nah?
def get_sig(my_cat=0, vbf=True, gf=True):
    if not vbf and not gf:
        sys.exit("Error: neither vbf nor gf signal requested.")
    fSIG = TFile.Open("%s/root/sig_shapes_workspace.root"%(workdir), "read")
    wSIG = fSIG.Get("w")
    sPDF = wSIG.pdf("signal_model_m%d_CAT%d"%(MASS,my_cat))
    yG = wSIG.var("yield_signalGF_mass125_CAT%d"%(my_cat)) #ggH yield
    yV = wSIG.var("yield_signalVBF_mass125_CAT%d"%(my_cat))
    x_name = "mbbReg_CAT%d"%my_cat
    sig_hist = sPDF.createHistogram(x_name, nbins)
    sig_vals = [0 for i in range(nbins)]
    gfness, vbfness = 0, 0
    if vbf:
        vbfness = 1
    if gf:
        gfness = 1
    for i in range(1, nbins+1):
        sig_vals[i-1] = sig_hist.GetBinContent(i) * (yG.getVal() * gfness + yV.getVal() * vbfness)

    return sig_vals

#function to plot differences between (data-Z-top), pol6 fit, and pol4 fit.
def plot_diffs(my_cat, spec_name="", gen_func="Pol6", fit_func="Pol4"):
    #style()
    #first run fitting method to get data to plot
    x_vals, p4_dat, p4_err, p6_dat, p6_err, datapt, daterr, z_vals, t_vals = do_fits(my_cat, gen_func=gen_func, fit_func=fit_func)
    #subtract the pol6 values to make this plot
        #Note: to get the full values (not differences) need to subtract z and top values also!!
    for i in range(nbins):
        p4_dat[i] = p4_dat[i] - p6_dat[i]
        datapt[i] = datapt[i] - z_vals[i] - t_vals[i] - p6_dat[i]
        p6_dat[i] = p6_dat[i] - p6_dat[i] #0 of course

    can_name = "cat%d_funcs_diff_%s"%(my_cat, spec_name)
    points = [datapt]
    errs = [daterr]
    curves = [p4_dat, p6_dat]
    point_names = [ "data-Z-top-%s"%(gen_func) ]
    curve_names = [ "%s - %s"%(fit_func, gen_func),
                    "%s - %s"%(gen_func, gen_func) ]
    return draw_plots(can_name, x_vals, curves=curves, points=points, pt_errs=errs, curve_names=curve_names, point_names=point_names)

#plot ratio of fits of pol6 and pol4 (as of yet untested as independent function!)
# note: currently useless!
def plot_64ratio(my_cat=0):
    #PLOTTING RATIO OF POL6 / POL4 NOW 
    #style()
    canrat = TCanvas("cat%d_Pol6_4_ratio"%my_cat,"cat%d_Pol6_4_ratio"%my_cat,600,600)
    canrat.cd()
#    pCMS1.Draw()
#    pCMS2.Draw()
    canrat.SetTitle("Pol3,4,5,6")
    canrat.SetFillColor(0)
    canrat.SetFillStyle(-1)
    nevents = 13142089 #from root file opened earlier
    #print("*****************Pol6 fit values:")
#    print("res: %s"%res["Pol6"].printArgs())
    coefs = res["Pol6"].floatParsFinal()
    varNarrow = RooRealVar("varNar", "varNar", 100, 99, 101)
    setNarrow = RooArgSet()
    setNarrow.add(varNarrow)
    valNarrow = model["Pol6"].getVal(setNarrow)
    #print("***GetVal (narrow range) result: %f" %valNarrow)
    x6name = x_name
    fitHist6 = model["Pol6"].createHistogram(x6name, nbins)
    fitHist6.Print()
    x4name = x_name
    fitHist4 = model["Pol4"].createHistogram(x4name, nbins)
    fitHist4.Print()
    ratio64 = fitHist6.Divide(fitHist4) #ratio64 is only a bool for some reason (ratio his is still in fitHist6)
    
    frameboi = x.frame()
    #draw ratio histogram with axes, curve style, with error bars.
    frameboi.GetXaxis().SetTitle("m_{bb} (GeV)")
    frameboi.GetYaxis().SetTitle("Ratio Pol6 / Pol4")
    frameboi.SetAxisRange(0.95, 1.05, "Y")  
    frameboi.Draw()
    fitHist6.Draw("same AC") #add "E" for error bars (they'll prolly be hella big)
    canrat.Modified()
    canrat.Update()

    canrat.SaveAs("%s/plot/draw_funcs/%s.pdf"%(workdir, canrat.GetName()))
    return

#helper function to write avg,err in each bin to file
#  format: bin_avg \t bin_err \n (for each of the numbins bins)
#  used in get_toys()
def write_avgs(filename, numbins, numtoys, totals, errs, zjets=[], top=[], ggH=[], qqH=[]):
    if numtoys == 0:
        sys.exit("Error: 0 valid toys.")
    store_file = open(filename, "w")
    for i in range(nbins):
        fin_avg = totals[i]/numtoys
        if zjets == []:
            fin_err = errs[i]/numtoys
            store_file.write("%f\t%f\n"%(fin_avg, fin_err))
        else:
            zjet_avg = zjets[i] / numtoys
            top_avg = top[i] / numtoys
            ggH_avg = ggH[i] / numtoys
            qqH_avg = qqH[i] / numtoys
            store_file.write("%f\t%f\t%f\t%f\t%f\n"%(fin_avg, zjet_avg, top_avg, ggH_avg, qqH_avg))
    store_file.close()

#to be called from within a try block
def end_program(message):
    print(message)
    sys.exit(0)

#bkg true means use total_background instead of qcd
def get_toys(my_cat=0, spec_name="noSys_shapes", raw=False, bkg=False, sig=False, error=False, fit_func="Pol4", gen_func="Pol6", vbf=False, ggH=False):
    if bkg and sig:
        sys.exit("Please specify either bkg or sig (or neither), not both.")
    if raw and (bkg or sig):
        sys.exit("Error: the raw toy data has no concept of bkg or signal.")
    #NOW PLOTTING TOY RATIO THINGY
    #style()
    #need to use all ~10000 toys and take their average in each bin.
    #average of and error on all toys, for each bin.
    toy_bin_avgs = [0 for i in range(nbins)]
    #for error need to get error on the fit--but how to do this??
    toy_bin_errs = [0 for i in range(nbins)]
    extraname = ""
    if bkg:
        extraname = "totbkg_"
    elif sig:
        extraname = "totsig_"
    elif vbf:
        extraname = "sigvbf_"
    elif ggH:
        extraname = "sigggH_"
    store_name = "total_toy_bins_%s%scat%d_%s_%dbins_%s_%s.txt"%(extraname, "raw" if raw else "", my_cat, spec_name,nbins, gen_func, fit_func)
    #if file storing averages exists, just read from it.
    if os.path.exists(store_name):
        print("Reading from file %s"%store_name)
        store_file = open(store_name, "read")
        for i,line in enumerate(store_file):
            if i > nbins:
                sys.exit("Error: more than nbins=%d lines in store file!"%nbins)
            words = line.split()
            toy_bin_avgs[i] = float(words[0])
            toy_bin_errs[i] = float(words[1])
    #####end store file exists
    #if the file doesn't exist, need to read from the root files and write the file.
    else:
        #list to hold totals for each bin (will be divided by tot_valid at end to get avg)
        bin_tots = [0 for i in range(nbins)]
        #total sum of squares of errors (will be divided by tot_valid at end to get avg error)
        bin_err2 = [0 for i in range(nbins)] 
        #total number of valid toys
        tot_valid = 0
        #number of different jobs to make/fit toys
        njobs = 500
        #number of times combine is called in each job (also number of toys per job)
        ncalls = 25
        for job in range(njobs):
            #write to and close the file every so often just in case job gets kilt
            if job % 50 == 49:
                write_avgs(store_name, nbins, tot_valid, bin_tots, bin_err2)
            for call in range(ncalls):
                toy_path = "/eos/user/b/bgreenbe/cat%d_%d/%s/%s_%s/job%d/fitDiagnostics_%d_%d.root"%(my_cat, my_cat, spec_name, gen_func, fit_func, job, job, call)
                if raw:
                            #toy file needs to be renamed to the call number instead of the random seed.
                    toy_path = "/eos/user/b/bgreenbe/cat%d_%d/%s/Pol6_Pol4/job%d/higgsCombine_%d.GenerateOnly.mH120.%d.root"%(my_cat, my_cat, spec_name, job, job, call)
                    
                try:
                    toy = TFile.Open(toy_path, "read")
                    if bkg:
                        qcd = toy.Get("shapes_fit_s/CAT%d/total_background"%my_cat)
                    elif sig:
                        qcd = toy.Get("shapes_fit_s/CAT%d/total_signal"%my_cat)
                    elif vbf:
                        qcd = toy.Get("shapes_fit_s/CAT%d/qqH_hbb"%my_cat)
                    elif ggH:
                        qcd = toy.Get("shapes_fit_s/CAT%d/ggH_hbb"%my_cat)
                    elif raw:
                        toydat = toy.Get("toys/toy_1")
                        #not so simple in this sh*tty version of root.
                        #qcd = toydat.createHistogram("mbbReg_CAT%d"%my_cat, nbins)
                        hist = toydat.binnedClone()
                 #       hist.Print()
                        qcd = hist.createHistogram("mbbReg_CAT%d"%my_cat, nbins)
                    #    qcd.Print()
                    #    print("bin 1: %f bin 50: %f bin 120: %f"%(qcd.GetBinContent(1), qcd.GetBinContent(50), qcd.GetBinContent(120)))
                   # print("obtained qcd: toy %s"%(toy_path))
                    else:
                        #qcd = toy.Get("shapes_fit_s/CAT%d/qcd"%my_cat)
                        qcd = toy.Get("shapes_fit_s/CAT%d/total"%my_cat)
                    nent_toy = int(qcd.GetEntries())
                    if bkg or sig or vbf or ggH:
                        nent_toy = 1200 #don't believe those f*ckers when they say it's 3600--the last 2400 entries are all 0
                    if raw:
                        nent_toy = nbins  #again, they don't want you to know the truth
                    for i in range(1, nbins+1):
                        toybinN = (i-1)*(nent_toy/nbins) + (nent_toy/nbins)/2 + 1
                        #print("nent_toy: %d toybinN: %d center of toybin: %f"%(nent_toy, toybinN, qcd.GetBinCenter(toybinN)))
                        toybin = qcd.GetBinContent(toybinN)
                        toyerr = qcd.GetBinError(toybinN)
                        bin_tots[i-1] += toybin
                        bin_err2[i-1] += (toyerr ** 2)
                    tot_valid += 1
                except:
                    print("Error: could not retrieve toy %s"%(toy_path))
                    signal(SIGINT, handler)
                    continue
            ###end of call loop
            print("done reading job %d"%job) 
        ###end of job loop
        #write final results to the file
        write_avgs(store_name, nbins, tot_valid, bin_tots, bin_err2)
        for i in range(nbins):
            toy_bin_avgs[i] = bin_tots[i] / tot_valid
            toy_bin_errs[i] = (bin_err2[i])**0.5 / tot_valid
    ###end store file does not yet exist
    print("done filling bins")
    #might need to subtract Z and top from the toy data,
    # but this function just returns the raw toy data.
    print("done getting toys")
    if not error:
        return toy_bin_avgs
    return toy_bin_avgs, toy_bin_errs

#compare the fitted signal of the signal model with the toys.
def plot_sigdiff(my_cat=0, spec_name="noSys_shapes", fit_func="Pol4", gen_func="Pol6"):
    #get the toy data: total_background.
    toy_bin_avgs = get_toys(my_cat, spec_name, bkg=False, vbf=True, fit_func=fit_func, gen_func=gen_func) #bpgballin
    #run fitting function to get the pol4 and pol6 fits
    x_vals, p4_dat, p4_err, p6_dat, p6_err, datapt, daterr, z_vals, t_vals = do_fits(my_cat, gen_func=gen_func, fit_func=fit_func)
    sig_vals = get_sig(my_cat, vbf=True, gf=False)
    
    #print("finished doing the fits.")
#    print("before doing shit: datapt=%s; toy_bin_avgs=%s"%(datapt, toy_bin_avgs))
    #now make the things I actually want to plot
    #difference in pol4_fit and toy avg; diff in pol6_fit and toy avg; diff in data and toy_avg
#    p4_diff = [0 for i in range(nbins)]
#    p6_diff = [0 for i in range(nbins)]
#    datdiff = [0 for i in range(nbins)]
#    #toydiff will always be 0
#    toydiff = [0 for i in range(nbins)]
#    for i in range(nbins):
#        #subtract z and top data from toys to get just qcd; add them to models to make up for datapts including z and top.
#        toy_bin_avgs[i] = toy_bin_avgs[i]
##        p4_diff[i] = p4_dat[i]# - toy_bin_avgs[i]
##        p6_diff[i] = p6_dat[i]# - toy_bin_avgs[i] 
##        datdiff[i] = datapt[i]# - z_vals[i] - t_vals[i] - toy_bin_avgs[i]
#        toydiff[i] = toy_bin_avgs[i] #- toy_bin_avgs[i]

    #group the like lists for plotting.
    #lists to be plotted as data points (with error bars)
    points = [sig_vals]
    #errors associated with each of the lists in points
#    errs = [daterr]
    #curves with no error bars
    curves = [toy_bin_avgs] 
    #curves = []
    #name for the canvas that we're gonna draw on.
    can_name = "cat%d_vbf_diff_%s_%s_%s"%(my_cat, spec_name, gen_func, fit_func)
    #list of names of each function (for the legend)
    curve_names = ["average toy fitted vbf, %s_%s CAT%d"%(gen_func, fit_func, my_cat)]
    point_names = [ "vbf model CAT%d"%my_cat ]
    #now draw the plots on a canvas with provided name.
    #print("boutta call the plotting function boiiiii")
    return draw_plots(can_name, x_vals, curves=curves, points=points, curve_names=curve_names, point_names=point_names, y_range=[0,100])
    

#plot the toys compared to the fits and data points.
def plot_toydiff(my_cat=0, spec_name="noSys_shapes", fit_func="Pol4", gen_func="Pol6"):
    #get the toy data: total_background.
    toy_bin_avgs = get_toys(my_cat, spec_name, bkg=True, sig=False, fit_func=fit_func, gen_func=gen_func) #bpgballin
    #run fitting function to get the pol4 and pol6 fits
    x_vals, p4_dat, p4_err, p6_dat, p6_err, datapt, daterr, z_vals, t_vals = do_fits(my_cat, gen_func=gen_func, fit_func=fit_func)
    
    #print("finished doing the fits.")
#    print("before doing shit: datapt=%s; toy_bin_avgs=%s"%(datapt, toy_bin_avgs))
    #now make the things I actually want to plot
    #difference in pol4_fit and toy avg; diff in pol6_fit and toy avg; diff in data and toy_avg
    p4_diff = [0 for i in range(nbins)]
    p6_diff = [0 for i in range(nbins)]
    datdiff = [0 for i in range(nbins)]
    #toydiff will always be 0
    toydiff = [0 for i in range(nbins)]
    for i in range(nbins):
        #subtract z and top data from toys to get just qcd; add them to models to make up for datapts including z and top.
        toy_bin_avgs[i] = toy_bin_avgs[i] - t_vals[i] - z_vals[i]
        p4_diff[i] = p4_dat[i] - toy_bin_avgs[i]
        p6_diff[i] = p6_dat[i] - toy_bin_avgs[i] 
        datdiff[i] = datapt[i] - z_vals[i] - t_vals[i] - toy_bin_avgs[i]
        toydiff[i] = toy_bin_avgs[i] - toy_bin_avgs[i]

    #group the like lists for plotting.
    #lists to be plotted as data points (with error bars)
    points = [datdiff]
    #errors associated with each of the lists in points
    errs = [daterr]
    #curves with no error bars
    curves = [p4_diff, p6_diff, toydiff] 
    #curves = []
    #name for the canvas that we're gonna draw on.
    can_name = "cat%d_toy_diff_%s"%(my_cat, spec_name)
    #list of names of each function (for the legend)
    curve_names = ["qcd_model_%s-(toy avg-Z-top)"%(fit_func),
                   "qcd_model_%s-(toy avg-Z-top)"%(gen_func),
                   "toy avg-Z-top-(toy avg-Z-top)" ]
    point_names = [ "data-Z-top-(toy avg-Z-top)" ]
    #now draw the plots on a canvas with provided name.
    #print("boutta call the plotting function boiiiii")
    return draw_plots(can_name, x_vals, curves=curves, points=points, pt_errs=errs, curve_names=curve_names, point_names=point_names)

#function to draw the actual plots.
#name is string (name for the canvas); x is list of x values.
#curve_names is list of names for curves and line_names is same for point sets (for the legend)
# all other inputs are lists of lists (and optional).
def draw_plots(name, x_vals, curves=[], cv_errs=[], points=[], pt_errs=[], curve_names=[], point_names=[], x_range=[], y_range=[]):
    print("items passed to draw_plots:\ncurves: %s \n points: %s \n names: %s%s"%( str(curves) , str(points) , str(curve_names) , str(point_names)) )
    #sequence of colors for each different dataset (may not all be used)
    colors = [kBlack, kBlue, kRed, kGreen, kMagenta, kOrange, kCyan, kViolet, kPink, kSpring]
    #convert to arrays for easy plotting 
    #list of curve arrays
    curve_arrs = []
    #list of curve error arrays
    curveEarrs = []
    for i in range(len(curves)):
        curve_arrs.append(array("f", curves[i]))
        if not cv_errs == []:
            curveEarrs.append(array("f", cv_errs[i]))
    #Now same for pointy bois
    point_arrs = []
    pointEarrs = []
    for i in range(len(points)):
        point_arrs.append(array("f", points[i]))
        if not pt_errs == []:
            pointEarrs.append(array("f", pt_errs[i]))
    #array of x values
    x_varr = array("f", x_vals)

    #xerr is 0 but is a required argument for the TGraph constructor for some stupid ass reason
    xerr = array('f', [0 for i in range(nbins)])

    #create canvas to draw on
    can2 = TCanvas(name,name,600,600)
    can2.cd(0)

    #draw cms official text at top
   # pCMS1.Draw()
   # pCMS2.Draw()
#    can2.SetTitle("Pol3,4,5,6")
    can2.SetFillColor(0)
    can2.SetFillStyle(-1)
    gPad.SetBottomMargin(0.3)
    #get the pad ready
    pad = TPad("pad","pad",0.,0.,1.,1.)
#    pad.SetTopMargin(0.7)
    pad.SetFillStyle(-1)
    pad.SetBorderSize(0)
    pad.Draw()
    pad.cd(0)
    #Make legend
    L00 = TLegend(0.50,0.85,1.0,1.0)
    L00.SetFillColor(0)
    L00.SetLineColor(1)
    L00.SetBorderSize(10)
    L00.SetTextFont(42)
    L00.SetTextSize(gStyle.GetPadTopMargin()*0.50)
    #Make and draw the necessary TGraphs (or TGraphErrors)
    tgraphs = []
    c = 0 #counter to keep track of how many tgraphs there are.
#Now same for points
    for i in range(len(point_arrs)):
        #if there are errors on some of the points, and this particular points set has errors, TGraphErrors.
        if pt_errs != [] and pt_errs[i] != []:
            tgraphs.append( TGraphErrors(nbins, x_varr, point_arrs[i], xerr, pointEarrs[i]) )
        #otherwise it's a normal TGraph.
        else:
            tgraphs.append( TGraph(nbins, x_varr, point_arrs[i]) )
        tgraphs[c].SetMarkerStyle(7)
        tgraphs[c].SetLineColor(colors[c])
        tgraphs[c].SetFillColor(colors[c])
        tgraphs[c].SetMarkerColor(colors[c])
        #can only use 'a' for axis for first tgraph, not the rest (root is not very smart)
        if c == 0:
            tgraphs[c].Draw("same APE0")
        else:
            tgraphs[c].Draw("same PE0")
        if point_names != []:
            L00.AddEntry(tgraphs[c], point_names[i], 'p')
        c += 1

    for i in range(len(curve_arrs)):
        #if there are errors on some of the curves, and this particular curve has errors, TGraphErrors.
        if cv_errs != [] and cv_errs[i] != []:
            tgraphs.append( TGraphErrors(nbins, x_varr, curve_arrs[i], xerr, curveEarrs[i]) )
        #otherwise it's a normal TGraph.
        else:
            tgraphs.append( TGraph(nbins, x_varr, curve_arrs[i]) )
        tgraphs[c].SetFillStyle(3005)
        tgraphs[c].SetLineColor(colors[c])
        tgraphs[c].SetFillColor(colors[c])
        #do the following only once
        if c == 0:
            tgraphs[c].Draw("same ace4")
        else:
            tgraphs[c].Draw("same ce4")
        if curve_names != []:
            L00.AddEntry(tgraphs[c], curve_names[i], 'l')
        c += 1

    #label axes
    newtg = tgraphs[0]
    print("Labelling axes!")
    newtg.GetXaxis().SetTitleSize(.055)
    newtg.GetXaxis().SetLabelSize(.035)
#zoom in on the beginning of the x-axis!
    if x_range != []:
        newtg.GetXaxis().SetRangeUser(x_range[0], x_range[1])
    if y_range != []:
        newtg.GetYaxis().SetRangeUser(y_range[0], y_range[1])

    newtg.GetYaxis().SetLabelSize(0.035)
    newtg.GetYaxis().SetTitleSize(0.055)
    newtg.GetYaxis().SetTitleOffset(1.5);
    #grdt.GetYaxis().SetTitle("Events / %.1f GeV"%(opts.dX[0]*10))
    newtg.GetYaxis().SetTitle("Events / %.1f GeV"%(dX))
    newtg.GetXaxis().SetTitle("M_{bb} (GeV)")
    L00.Draw("same")
    gPad.Update()
    gPad.Modified()
    #now save canvas to file so can view.
    path = "%s/plot/draw_funcs/%s.pdf"%(workdir,can2.GetName())
    can2.SaveAs(path)
#    can2.SaveAs("%s/plot/draw_funcs/%s.png"%(workdir,can2.GetName()))
    print("Saved canvas %s"%path)
    return path

#plot differences in Z and top models
def plot_zt_diff(my_cat=0, spec_name="noSys_shapes"):

    toy_bin_avgs = get_toys(my_cat, spec_name)
    #run fitting function to get the pol4 and pol6 fits
    x_vals, p4_dat, p4_err, p6_dat, p6_err, datapt, daterr, z_vals, t_vals = do_fits(my_cat)
    
#    print("before doing shit: datapt=%s; toy_bin_avgs=%s"%(datapt, toy_bin_avgs))
    #now make the things I actually want to plot
    #difference in pol4_fit and toy avg; diff in pol6_fit and toy avg; diff in data and toy_avg
    p4_datdiff = [0 for i in range(nbins)]
    p6_datdiff = [0 for i in range(nbins)]
    p4_toydiff = [0 for i in range(nbins)]
    p6_toydiff = [0 for i in range(nbins)]
    zt_dat = [0 for i in range(nbins)]
    #toydiff will always be 0
#    toydiff = [0 for i in range(nbins)]
    for i in range(nbins):
        #subtract z and top data from toys to get just qcd; add them to models to make up for datapts including z and top.
        toy_bin_avgs[i] = toy_bin_avgs[i] #- z_vals[i] - t_vals[i]
        p4_datdiff[i] = datapt[i] - p4_dat[i] #p4_dat[i] - toy_bin_avgs[i]
        p6_datdiff[i] = datapt[i] - p6_dat[i] #p6_dat[i] - toy_bin_avgs[i] 
        p4_toydiff[i] = toy_bin_avgs[i] - p4_dat[i]  #what's left is Z+top #- z_vals[i] - t_vals[i] - toy_bin_avgs[i]
        p6_toydiff[i] = toy_bin_avgs[i] - p6_dat[i]  #what's left is Z+top #- z_vals[i] - t_vals[i] - toy_bin_avgs[i]
        zt_dat[i] = z_vals[i] + t_vals[i]
#        toydiff[i] = toy_bin_avgs[i] - toy_bin_avgs[i]
    
    points = [p4_datdiff, p6_datdiff]
    errs = [daterr, daterr]
    curves = [p4_toydiff, p6_toydiff, zt_dat]
    point_names = [ "data-qcd_model_Pol4_CAT%d"%(my_cat),
                    "data-qcd_model_Pol6_CAT%d"%(my_cat) ]
    curve_names = [ "toy avg-qcd__model_Pol4_CAT%d"%(my_cat),
                    "toy avg-qcd_model_Pol6_CAT%d"%(my_cat),
                    "Z + top model (CAT%d)"%(my_cat) ]
    can_name = "cat%d_zt_diff_%s"%(my_cat, spec_name)
    return draw_plots(can_name, x_vals, points=points, curves=curves, pt_errs=errs, curve_names=curve_names, point_names=point_names)

####################################################################################################
def plot_topdiff(my_cat=0, spec_name="noSys_shapes"):
    #get toys
    toy_bin_avgs = get_toys(my_cat, spec_name)
    #run fitting function to get the pol4 and pol6 fits
    x_vals, p4_dat, p4_err, p6_dat, p6_err, datapt, daterr, z_vals, t_vals = do_fits(my_cat)
    
#    print("before doing shit: datapt=%s; toy_bin_avgs=%s"%(datapt, toy_bin_avgs))
    #now make the things I actually want to plot
    #difference in pol4_fit and toy avg; diff in pol6_fit and toy avg; diff in data and toy_avg
    p4_datdiff = [0 for i in range(nbins)]
    p6_datdiff = [0 for i in range(nbins)]
    p4_toydiff = [0 for i in range(nbins)]
    p6_toydiff = [0 for i in range(nbins)]
    zt_dat = [0 for i in range(nbins)]
    #toydiff will always be 0
#    toydiff = [0 for i in range(nbins)]
    for i in range(nbins):
        #subtract z and top data from toys to get just qcd; add them to models to make up for datapts including z and top.
        toy_bin_avgs[i] = toy_bin_avgs[i] #- z_vals[i] - t_vals[i]
        p4_datdiff[i] = datapt[i] - p4_dat[i] -z_vals[i]#p4_dat[i] - toy_bin_avgs[i]
        p6_datdiff[i] = datapt[i] - p6_dat[i] -z_vals[i]#p6_dat[i] - toy_bin_avgs[i] 
        p4_toydiff[i] = toy_bin_avgs[i] - p4_dat[i] - z_vals[i]  #what's left is top #- z_vals[i] - t_vals[i] - toy_bin_avgs[i]
        p6_toydiff[i] = toy_bin_avgs[i] - p6_dat[i] - z_vals[i]  #what's left is top #- z_vals[i] - t_vals[i] - toy_bin_avgs[i]
        zt_dat[i] = t_vals[i]
#        toydiff[i] = toy_bin_avgs[i] - toy_bin_avgs[i]

    points = [p4_datdiff, p6_datdiff]
    errs = [daterr, daterr]
    curves = [p4_toydiff, p6_toydiff, zt_dat]
    point_names = [ "data-Z-qcd_model_Pol4_CAT%d"%(my_cat),
                    "data-Z-qcd_model_Pol6_CAT%d"%(my_cat) ]
    curve_names = [ "toy avg-Z-qcd__model_Pol4_CAT%d"%(my_cat),
                    "toy avg-Z-qcd_model_Pol6_CAT%d"%(my_cat),
                    "top model (CAT%d)"%(my_cat) ]
    can_name = "cat%d_topdiff_%s"%(my_cat, spec_name)
    return draw_plots(can_name, x_vals, points=points, curves=curves, pt_errs=errs, curve_names=curve_names, point_names=point_names)

def plot_Z_diff(my_cat=0, spec_name="noSys_shapes"):
    toy_bin_avgs = get_toys(my_cat, spec_name)
    #run fitting function to get the pol4 and pol6 fits
    x_vals, p4_dat, p4_err, p6_dat, p6_err, datapt, daterr, z_vals, t_vals = do_fits(my_cat)
    
#    print("before doing shit: datapt=%s; toy_bin_avgs=%s"%(datapt, toy_bin_avgs))
    #now make the things I actually want to plot
    #difference in pol4_fit and toy avg; diff in pol6_fit and toy avg; diff in data and toy_avg
    p4_datdiff = [0 for i in range(nbins)]
    p6_datdiff = [0 for i in range(nbins)]
    p4_toydiff = [0 for i in range(nbins)]
    p6_toydiff = [0 for i in range(nbins)]
    zt_dat = [0 for i in range(nbins)]
    #toydiff will always be 0
#    toydiff = [0 for i in range(nbins)]
    for i in range(nbins):
        #subtract z and top data from toys to get just qcd; add them to models to make up for datapts including z and top.
        toy_bin_avgs[i] = toy_bin_avgs[i] #- z_vals[i] - t_vals[i]
        p4_datdiff[i] = datapt[i] - p4_dat[i] -t_vals[i]#p4_dat[i] - toy_bin_avgs[i]
        p6_datdiff[i] = datapt[i] - p6_dat[i] -t_vals[i]#p6_dat[i] - toy_bin_avgs[i] 
        p4_toydiff[i] = toy_bin_avgs[i] - p4_dat[i] - t_vals[i]  #what's left is top #- z_vals[i] - t_vals[i] - toy_bin_avgs[i]
        p6_toydiff[i] = toy_bin_avgs[i] - p6_dat[i] - t_vals[i]  #what's left is top #- z_vals[i] - t_vals[i] - toy_bin_avgs[i]
        zt_dat[i] = z_vals[i]
#        toydiff[i] = toy_bin_avgs[i] - toy_bin_avgs[i]

    points = [p4_datdiff, p6_datdiff]
    errs = [daterr, daterr]
    curves = [p4_toydiff, p6_toydiff, zt_dat]
    point_names = [ "data-top-qcd_model_Pol4_CAT%d"%(my_cat),
                    "data-top-qcd_model_Pol6_CAT%d"%(my_cat) ]
    curve_names = [ "toy avg-top-qcd__model_Pol4_CAT%d"%(my_cat),
                    "toy avg-top-qcd_model_Pol6_CAT%d"%(my_cat),
                    "Z model (CAT%d)"%(my_cat) ]
    can_name = "cat%d_Z_diff_%s"%(my_cat, spec_name)
    return draw_plots(can_name, x_vals, points=points, curves=curves, pt_errs=errs, curve_names=curve_names, point_names=point_names)

#function to find the difference between the raw (unfitted) toy data and the actual data (which it should be very close to, I think.)
def toy_comp(my_cat=0, spec_name = "noSys_shapes"):
    #first get qcd fits and data data.
    x_vals, p4_dat, p4_err, p6_dat, p6_err, datapt, daterr, z_vals, t_vals = do_fits(my_cat)
    #now get toy (raw) data.
    raw_toys, raw_errs = get_toys(my_cat, spec_name, raw=True, error=True)
    fit_toys = get_toys(my_cat, spec_name, raw=False)

    raw_toys[nbins-1] = raw_toys[nbins-2] #small mistake in reading the data
    for i in range(nbins):
        fit_toys[i] = fit_toys[i] + z_vals[i] + t_vals[i] #- raw_toys[i]
        datapt[i] = datapt[i] #- raw_toys[i]
        raw_toys[i] = raw_toys[i] #- raw_toys[i]
    points = [datapt, raw_toys]
    errs = [daterr, []] #raw_errs]
    curves = [fit_toys]
    point_names = ["data-raw_toys", "raw_toys(gen w/ Pol6)-raw_toys"] 
    curve_names = ["Fit toys(fit w/ Pol4)-raw_toys"]
    can_name = "cat%d_toycomp_%s"%(my_cat, spec_name)
    return draw_plots(can_name, x_vals, points=points, curves=curves, pt_errs=errs, curve_names=curve_names, point_names=point_names) #y_range=[-1800,500])

#get all the individual backgrounds and signal from one toy
def get_toy_all(my_cat=0, spec_name="", gen_func="Pol6", fit_func="Pol4", job=0, call=0):
   # tf = TFile.Open("%s/root/data_shapes_workspace.root"%(workdir),"read")
   # w = tf.Get("w")
   # dset = w.data("data_hist_CAT0")
#    print("my_cat=%d, spec_name=%s, gen_func=%s, fit_func=%s, job=%d, call=%d"%(my_cat, spec_name, gen_func, fit_func, job, call))
    toy = TFile.Open("/eos/user/b/bgreenbe/cat%d_%d/%s/%s_%s/job%d/fitDiagnostics_%d_%d.root"%(my_cat, my_cat, spec_name, gen_func, fit_func, job, job, call), "read")
    toy_qcd = toy.Get("shapes_fit_s/CAT%d/qcd"%my_cat)
    toy_Z = toy.Get("shapes_fit_s/CAT%d/zjets"%my_cat)
    toy_top = toy.Get("shapes_fit_s/CAT%d/top"%my_cat)
    toy_ggH = toy.Get("shapes_fit_s/CAT%d/ggH_hbb"%my_cat)
    toy_qqH = toy.Get("shapes_fit_s/CAT%d/qqH_hbb"%my_cat)
    
    t_qcd = [0 for i in range(nbins)]
    t_Z = [0 for i in range(nbins)]
    t_top = [0 for i in range(nbins)]
    t_ggH = [0 for i in range(nbins)]
    t_qqH = [0 for i in range(nbins)]
    nent_toy = 1200
    for i in range(1, nbins+1):
        n = (i-1)*(nent_toy/nbins) + (nent_toy/nbins)/2 + 1
        t_qcd[i-1] = toy_qcd.GetBinContent(n)
        t_Z[i-1] = toy_Z.GetBinContent(n)
        t_top[i-1] = toy_top.GetBinContent(n)
        t_ggH[i-1] = toy_ggH.GetBinContent(n)
        t_qqH[i-1] = toy_qqH.GetBinContent(n)
#    print("done reading from toy file.")
    return t_qcd, t_Z, t_top, t_ggH, t_qqH #added t_qqH recently!

def handler():
    sys.exit(0)

#get the averages over all toys for all of the processes
def obtain_toy_all(my_cat=0, spec_name="", gen_func="Pol6", fit_func="Pol4", njobs=500, ncalls=25):
    #filename to read from or write to
    filename = "cat%d_%s_%s_%s_toy_all.txt"%(my_cat, spec_name, gen_func, fit_func)
    qcd = [0.0 for i in range(nbins)]
    zjets= [0.0 for i in range(nbins)]
    top = [0.0 for i in range(nbins)]
    ggH = [0.0 for i in range(nbins)]
    sig = [0.0 for i in range(nbins)]
    if os.path.exists(filename):
        print("reading from %s"%filename)
        #read from file
        fin = open(filename, "r")
        for n in range(nbins):
            line = fin.readline()
            datums = line.split()
            qcd[n] = float(datums[0])
            zjets[n] = float(datums[1])
            top[n] = float(datums[2])
            ggH[n] = float(datums[3])
            sig[n] = float(datums[4])
    else:
        print("reading from root files")
        fout = open(filename, "w")
        numtoys = 0
        for job in range(njobs):
            for call in range(ncalls):
                try:
                    t_qcd, t_Z, t_top, t_ggH, t_qqH = get_toy_all(my_cat=my_cat, spec_name=spec_name, job=job, call=call, gen_func=gen_func, fit_func=fit_func)
                    for n in range(nbins):
                        qcd[n] += t_qcd[n]
                        zjets[n] += t_Z[n]
                        top[n] += t_top[n]
                        ggH[n] += t_ggH[n]
                        sig[n] += t_qqH[n]
                    numtoys += 1
                except:
                    signal(SIGINT, handler)
                    print("Error: could not read file for job %d, call %d"%(job, call))
            #periodically write results to output file in case of crash
            print("finished reading job %d out of %d"%(job+1, njobs))
            try:
#            if 0 == 0:
                write_avgs(filename, nbins, numtoys, qcd, [], zjets=zjets, top=top, ggH=ggH, qqH=sig)
            except:
#            else:
                print("Warning: 0 valid toys in job %d"%job)
        #now that all the results are read, can do the final division by numtoys to get the averages.
        if numtoys == 0:
            os.system("rm %s"%filename)
            sys.exit("Error: 0 valid toys")
        else:
            for n in range(nbins):
                qcd[n] /= numtoys
                zjets[n] /= numtoys
                top[n] /= numtoys
                ggH[n] /= numtoys
                sig[n] /= numtoys
    return qcd, zjets, top, ggH, sig
    
#this function plots the average Z and top for all toys vs the Z and top models.
def plot_toy_all(my_cat=0, spec_name="", gen_func="Pol6", fit_func="Pol4", TF="ConstPol1"):
    t_qcd, t_Z, t_top, t_ggH, t_qqH = obtain_toy_all(my_cat=my_cat, spec_name=spec_name, gen_func=gen_func, fit_func=fit_func)
    x_vals, p4_dat, p4_err, p6_dat, p6_err, datapt, daterr, z_vals, t_vals = do_fits(my_cat, gen_func=gen_func, fit_func=fit_func, TF=TF)
    #curves = [p4_dat, p6_dat, z_vals, t_vals, t_qcd, t_Z, t_top]
    curves = [ z_vals, t_vals, t_Z, t_top]
    curve_names = [ "Z model", "top model", "avg toy Z", "avg toy top"]
    can_name = "toy_Z_top_CAT%d_%s"%(my_cat, spec_name)
    return draw_plots(can_name, x_vals, curves=curves, curve_names=curve_names) #, y_range=[0, 62000])

#this function plots the differences bt pol4,5,6 models
def plot_toy_full(my_cat=0, spec_name="", gen_func="Pol6", fit_func="Pol4", sig_strength=1.0, TF="ConstPol1"):
    #get toy data
    t_qcd, t_Z, t_top, t_ggH, t_qqH = obtain_toy_all(my_cat=my_cat, spec_name=spec_name, gen_func=gen_func, fit_func=fit_func)
    #get model data
    x_vals, p4_dat, p4_err, p6_dat, p6_err, datapt, daterr, z_vals, t_vals = do_fits(my_cat, gen_func="Pol6", fit_func="Pol4", TF=TF)
    x_vals, p4_dat, p4_err, p5_dat, p5_err, datapt, daterr, z_vals, t_vals = do_fits(my_cat, gen_func="Pol5", fit_func="Pol4", TF=TF)
    #get signal models for gF and qqH
    g_vals = get_sig(my_cat, gf=True, vbf=False)
    s_vals = get_sig(my_cat, gf=False, vbf=True)
    tot_toy = [0 for i in range(nbins)]
    totPol4 = [0 for i in range(nbins)]
    totPol5 = [0 for i in range(nbins)]
    totPol6 = [0 for i in range(nbins)]
    for i in range(nbins):
        tot_toy[i] = t_qcd[i] #+ t_Z[i] + t_top[i] + t_ggH[i] + t_qqH[i]
        totPol4[i] = p4_dat[i] #+ z_vals[i] + t_vals[i] + g_vals[i] + s_vals[i]
        totPol5[i] = p5_dat[i] #+ z_vals[i] + t_vals[i] + g_vals[i] + s_vals[i]
        totPol6[i] = p6_dat[i] #+ z_vals[i] + t_vals[i] + g_vals[i] + s_vals[i]
    
        tot_toy[i] -= totPol4[i]
        totPol5[i] -= totPol4[i]
        totPol6[i] -= totPol4[i]
        datapt[i] -= totPol4[i]
        totPol4[i] = 0

        s_vals[i] = sig_strength*(s_vals[i] + g_vals[i])
    
    #print("tot_toy: %s"%tot_toy)
    #print("totPol4: %s"%totPol4)
    points = [] #[s_vals] #[datapt]
    errs = [] #[daterr]
    curves = [tot_toy, totPol4, totPol5, totPol6]
    #point_names = ["data-Pol4model"]
    point_names = ["signal model *%f"%sig_strength]
    curve_names = ["toy fitted qcd(%s_%s)-Pol4model"%(gen_func, fit_func), "Pol4model-Pol4model", "Pol5model-Pol4model", "Pol6model-Pol4model"]
    can_name = "tot_toy_%s_%s_%s_v_mod_CAT%d"%(spec_name, gen_func, fit_func, my_cat)
    return draw_plots(can_name, x_vals, points=points, pt_errs=errs, curves=curves, curve_names=curve_names, y_range=[-15, 15])

def plot_sigbois(my_cat=0, spec_name="", gen_func="Pol6", fit_func="Pol4", sig_strength=1.0, TF="ConstPol1"):
    #get all models
   # t_qcd, t_Z, t_top, t_ggH, t_qqH = obtain_toy_all(my_cat=my_cat, spec_name=spec_name, gen_func=gen_func, fit_func=fit_func)
    x_vals, p4_dat, p4_err, p6_dat, p6_err, datapt, daterr, z_vals, t_vals = do_fits(my_cat, gen_func="Pol6", fit_func=fit_func, TF=TF)

#whoops forgot to include signal in the toy_all thing haha----NOT HAHA YOU IDIOT YOU JUST WASTED A WHOLE HOUR OF MY TIME A$$HOLE - future you
#    t_vbf = get_toys(my_cat, spec_name, sig=True, fit_func=fit_func, gen_func=gen_func)
    t_qcd, t_Z, t_top, t_ggH, t_qqH = obtain_toy_all(my_cat=my_cat, spec_name=spec_name, gen_func=gen_func, fit_func=fit_func)

    #get full signal (vbf+gf)
    s_vals = get_sig(my_cat, vbf=True, gf=False)
    totsig = [0 for i in range(nbins)]
    for i in range(nbins):
        totsig[i] = t_qqH[i]  #t_vbf[i] #t_ggH[i] + 
        
        s_vals[i] *= sig_strength

    points = [s_vals]
    curves = [totsig]                      
    point_names = ["signal model"]
    curve_names = ["toy fitted signal (%s_%s)"%(gen_func, fit_func)]
    can_name = "sig_%s_%s_%s_CAT%d"%(spec_name, gen_func, fit_func, my_cat)
    return draw_plots(can_name, x_vals, points=points, point_names=point_names, curves=curves, curve_names=curve_names) #, y_range=[-1, 400])
    

#plot difference between Pol6,Pol4; along with diff bt Pol5,Pol4 (along with signal model)
def plot_qcddiff(my_cat=0):
    x_vals, q4_dat, q4_err, q5_dat, q5_err, datapt, daterr, z_vals, t_vals = do_fits(my_cat, gen_func="Pol5", fit_func="Pol4")
    x_vals, q4_dat, q4_err, q6_dat, q6_err, datapt, daterr, z_vals, t_vals = do_fits(my_cat, gen_func="Pol6", fit_func="Pol4")
    sig_vals = get_sig(my_cat)
#    toy_sig = get_toys(my_cat=my_cat, sig=True, fit_func
    #curves = [p4_dat, p6_dat, z_vals, t_vals, t_qcd, t_Z, t_top]
    p4_diff = [0 for i in range(nbins)]
    p5_diff = [0 for i in range(nbins)]
    p6_diff = [0 for i in range(nbins)]
    datdiff = [0 for i in range(nbins)]
    daterr_shrunk = [0 for i in range(nbins)]
    for i in range(nbins):
        p6_diff[i] = q6_dat[i] - q4_dat[i]
        p5_diff[i] = q5_dat[i] - q4_dat[i] 
        datdiff[i] = datapt[i] - z_vals[i] - t_vals[i] - q4_dat[i]
        daterr_shrunk[i] = daterr[i] / 50
        
    #points = [sig_vals]
    points = [datdiff]
    errs = [daterr_shrunk]
    curves = [sig_vals, p4_diff, p5_diff, p6_diff]
    point_names = ["data-Z-top-qcdModel_Pol4"]
    curve_names = ["signal model CAT2", "qcdModel_Pol4-qcdModel_Pol4", "qcdModel_Pol5-qcdModel_Pol4", "qcdModel_Pol6-qcdModel_Pol4"]
    can_name = "qcd_diff_dat_CAT%d"%my_cat
    return draw_plots(can_name, x_vals, points=points, pt_errs=errs, curves=curves, curve_names=curve_names, point_names=point_names, y_range=[-20, 20])

def main():
    style()
    my_cat = 5
    #special name of the directory to find the toy shapes
   # spec_name = "noSys_shapes"
    #spec_name = "cminD0"
    #spec_name = "preciseRate"
    #spec_name = "perfNoSys"
    #spec_name = "expectSig50"
    #spec_name = "noSys50"
    spec_name = "TFPOL1_shapes"
    #spec_name = "useBkg" 
    #path = plot_diffs(my_cat) #, "useBkg")
    #plot_zt_diff(my_cat, spec_name)
    #plot_topdiff(my_cat, spec_name)
    #path = plot_Z_diff(my_cat, spec_name)
    gen_func = "Pol5"
    fit_func = "Pol4"
    #path = plot_toydiff(my_cat, spec_name, gen_func=gen_func, fit_func=fit_func)
    #path = plot_qcddiff(2)
    #path = plot_sigdiff(my_cat, spec_name, gen_func=gen_func, fit_func=fit_func)
#    path = toy_comp(my_cat, spec_name)
    #open new pdf created
    #path = plot_toy_all(my_cat, spec_name, gen_func=gen_func, fit_func=fit_func, TF="POL1")
    #path = plot_toy_full(my_cat, spec_name, gen_func=gen_func, fit_func=fit_func, sig_strength=1.0, TF="POL1")
    path = plot_sigbois(my_cat, spec_name, gen_func=gen_func, fit_func=fit_func, sig_strength=1.0, TF="POL1")
    os.system("evince %s"%path)

if __name__=='__main__':
    main()
