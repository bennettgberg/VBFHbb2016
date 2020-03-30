#!/usr/bin/env python

import os,sys,re,json,datetime
from glob import glob
from array import array
from optparse import OptionParser,OptionGroup
import warnings

from  pdf_param_cfi import *
from myGenPdf import *
from generateFormula import *

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

def do_fits(my_cat=0):
    fBKG = TFile.Open("%s/root/bkg_shapes_workspace.root"%(workdir),"read")
    fSIG = TFile.Open("%s/root/sig_shapes_workspace.root"%(workdir), "read")
    fqcd4 = TFile.Open("%s/root/bias_shapes_workspace_Pol4.root"%(workdir), "read")
    fqcd6 = TFile.Open("%s/root/bias_shapes_workspace_Pol6.root"%(workdir), "read")
    wBKG = fBKG.Get("w")
    wSIG = fSIG.Get("w")
    wqcd4 = fqcd4.Get("w")
    wqcd6 = fqcd6.Get("w")
    w    = RooWorkspace("w","workspace")

    wqcd4.Print()
    wqcd6.Print()

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
    sPDF = wSIG.pdf("signal_model_mH%d_sel%s_CAT%d"%(MASS,tag,C))
    #print("literally boutta get the qPDF4 ")
    qPDF4= wqcd4.pdf("qcd_model_Pol4_CAT%d"%(C))
    qPDF6= wqcd6.pdf("qcd_model_Pol6_CAT%d"%(C))
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
        p6_dat[i-1] = q6vals[i-1]#pol6_val #- z_val*yZ.getVal() - top_val*yT.getVal()
        datapt[i-1] = data_val #- z_val*yZ.getVal() - top_val*yT.getVal()
        p4_dat[i-1] = q4vals[i-1] #pol4_val #- z_val*yZ.getVal() - top_val*yT.getVal()
        #also need errors
        daterr[i-1] = dat_hist.GetBinError(i)
      #  p4_err[i-1] = pol4_fitHist.GetBinError(i)
      #  p6_err[i-1] = pol6_fitHist.GetBinError(i)
        p4_err[i-1] = qcd4_hist.GetBinError(i) *p4_dat[i-1]#*yQ4.getVal()
        p6_err[i-1] = qcd6_hist.GetBinError(i) *p6_dat[i-1]#*yQ6.getVal()

    return x_vals, p4_dat, p4_err, p6_dat, p6_err, datapt, daterr, z_vals, t_vals 

#function to plot differences between (data-Z-top), pol6 fit, and pol4 fit.
def plot_diffs(my_cat):
    style()
    #first run fitting method to get data to plot
    x_vals, p4_dat, p4_err, p6_dat, p6_err, datapt, daterr, z_vals, t_vals = do_fits(my_cat)
    #subtract the pol6 values to make this plot
        #Note: to get the full values (not differences) need to subtract z and top values also!!
    for i in range(nbins):
        p4_dat[i] = p4_dat[i] - p6_dat[i]
        datapt[i] = datapt[i] - p6_dat[i]
        p6_dat[i] = p6_dat[i] - p6_dat[i] #0 of course
    #now convert lists to arrays so they can be plotted easily
    x_varr = array("f", x_vals)
    p4_arr = array("f", p4_dat)
    p4earr = array("f", p4_err)
    p6_arr = array("f", p6_dat)
    p6earr = array("f", p6_err)
    datarr = array("f", datapt)
    dtearr = array("f", daterr)
    #xerr is 0 but is a required argument for the TGraph constructor for some stupid ass reason
    xerr = array('f', [0 for i in range(nbins)])

    #now plot these bois (with error of course)
    can2 = TCanvas("cat%d_funcs_diff"%(my_cat),"cat%d_funcs_diff"%(my_cat),600,600)
    #pol4 graph
    #draw cms official text at top
    can2.cd(0)
   # pCMS1.Draw()
   # pCMS2.Draw()
    can2.SetTitle("Pol3,4,5,6")
    can2.SetFillColor(0)
    can2.SetFillStyle(-1)
    gPad.SetBottomMargin(0.3)
    print("x_varr: %s \n p4_arr: %s \n p4earr:%s \n lengths: %d, %d, %d"%(x_varr, p4_arr, p4earr, len(x_varr), len(p4_arr), len(p4earr)))
    print("x_varr: %s \n p6_arr: %s \n p6earr:%s \n lengths: %d, %d, %d"%(x_varr, p6_arr, p6earr, len(x_varr), len(p6_arr), len(p6earr)))
    print("x_varr: %s \n datarr: %s \n dtearr:%s \n lengths: %d, %d, %d"%(x_varr, datarr, dtearr, len(x_varr), len(datarr), len(dtearr)))
    grp4 = TGraphErrors(nbins, x_varr, p4_arr, xerr, p4earr)
    #pol6 graph
    grp6 = TGraphErrors(nbins, x_varr, p6_arr, xerr, p6earr)
    #data graph
    grdt = TGraphErrors(nbins, x_varr, datarr, xerr, dtearr)
    
    #width of each bin
#    #ndf = int((opts.X[1]-opts.X[0])/(opts.dX[0])) - (n_param)
#    ndf = int((xmax-xmin)/dX) - n_param
#    chi2=chi2_val/ndf
#    prob = ROOT.TMath.Prob(chi2_val,ndf)
# part 2
    pad = TPad("pad","pad",0.,0.,1.,1.)
#    pad.SetTopMargin(0.7)
    pad.SetFillStyle(-1)
    pad.SetBorderSize(0)
    pad.Draw()
    pad.cd(0)
#    i = n_param
    grdt.SetMarkerStyle(7)
#    grdt.GetYaxis().SetRangeUser(41000, 45000) 
#    grdt.GetXaxis().SetRangeUser(120, 130) 
    grdt.GetXaxis().SetTitleSize(.055)
    grdt.GetXaxis().SetLabelSize(.035)
    grdt.GetYaxis().SetLabelSize(0.035)
    grdt.GetYaxis().SetTitleSize(0.055)
    grdt.GetYaxis().SetTitleOffset(1.5);
    #grdt.GetYaxis().SetTitle("Events / %.1f GeV"%(opts.dX[0]*10))
    grdt.GetYaxis().SetTitle("Events / %.1f GeV"%(dX))
    grdt.GetXaxis().SetTitle("M_{bb} (GeV)")
    grdt.Draw("same APE0")
    #grdt.SetLineColor(kBlack)
    grp4.Draw("same ce4")
    grp4.SetLineColor(kBlue)
    grp4.SetFillColor(kBlue)
    grp4.SetFillStyle(3005)
    grp6.Draw("same ce4")
    grp6.SetLineColor(kRed)
    grp6.SetFillColor(kRed)
    grp6.SetFillStyle(3005)

    gPad.Update()

    #L00 = TLegend(0.7,0.6,1.-gStyle.GetPadRightMargin()-gStyle.GetPadTopMargin()*0.3333,1.-gStyle.GetPadTopMargin()*1.3333)
    L00 = TLegend(0.75,0.85,1.0,1.0)
    L00.SetFillColor(0)
    L00.SetLineColor(1)
    L00.SetBorderSize(10)
    L00.SetTextFont(42)
    L00.SetTextSize(gStyle.GetPadTopMargin()*0.50)
#    L00.SetY1(L00.GetY2()-L00.GetNRows()*gStyle.GetPadTopMargin()*0.95)
    L00.AddEntry(grp4, "Pol4 - Pol6", 'l')
    L00.AddEntry(grp6, "Pol6 - Pol6", 'l')
    L00.AddEntry(grdt, "data-Z-top-Pol6", 'p')
#    L00.AddEntry(dat_hist, "data", 'p')
#    dat_hist.SetMarkerStyle(7)
#    dat_hist.SetMarkerColor(kGreen)
    #draw raw data (without subtracting other background)
    #dat_hist.Draw("same pe")
    L00.Draw("same")
    print("drew legend")
    #now save canvas to file so can view.
    can2.SaveAs("%s/plot/draw_funcs/%s.pdf"%(workdir,can2.GetName()))
    can2.SaveAs("%s/plot/draw_funcs/%s.png"%(workdir,can2.GetName()))

    print("saved my boi")
#    sys.exit(0)
    return

#plot ratio of fits of pol6 and pol4 (as of yet untested as independent function!)
def plot_64ratio(my_cat=0):
    #PLOTTING RATIO OF POL6 / POL4 NOW 
    style()
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
#  used in plot_toydiff()
def write_avgs(filename, numbins, numtoys, totals, errs):
    if numtoys == 0:
        sys.exit("Error: 0 valid toys.")
    store_file = open(filename, "w")
    for i in range(nbins):
        fin_avg = totals[i]/numtoys
        fin_err = errs[i]/numtoys
        store_file.write("%f\t%f\n"%(fin_avg, fin_err))
    store_file.close()

def plot_toydiff(my_cat=0, spec_name="noSys_shapes"):
    #NOW PLOTTING TOY RATIO THINGY
    style()
    #need to use all ~10000 toys and take their average in each bin.
    #average of and error on all toys, for each bin.
    toy_bin_avgs = [0 for i in range(nbins)]
    #for error need to get error on the fit--but how to do this??
    toy_bin_errs = [0 for i in range(nbins)]

    store_name = "total_toy_bins_cat%d_%s_%dbins.txt"%(my_cat, spec_name,nbins)
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
                toy_path = "/eos/user/b/bgreenbe/cat%d_%d/%s/Pol6_Pol4/job%d/fitDiagnostics_%d_%d.root"%(my_cat, my_cat, spec_name, job, job, call)
                try:
                    toy = TFile.Open(toy_path, "read")
                    qcd = toy.Get("shapes_fit_s/CAT%d/qcd"%my_cat)
                   # print("obtained qcd: toy %s"%(toy_path))
                    nent_toy = int(qcd.GetEntries())
                    for i in range(1, nbins+1):
                        toybinN = (i-1)*(nent_toy/nbins) + (nent_toy/nbins)/2
                        #print("nent_toy: %d toybinN: %d center of toybin: %f"%(nent_toy, toybinN, qcd.GetBinCenter(toybinN)))
                        toybin = qcd.GetBinContent(toybinN)
                        toyerr = qcd.GetBinError(toybinN)
                        bin_tots[i-1] += toybin
                        bin_err2[i-1] += (toyerr ** 2)
                    tot_valid += 1
                except:
                    print("Error: could not retrieve toy %s"%(toy_path))
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
    #need to subtract Z and top from the toy data.

    #run fitting function to get the pol4 and pol6 fits
    x_vals, p4_dat, p4_err, p6_dat, p6_err, datapt, daterr, z_vals, t_vals = do_fits(my_cat)
    
#    print("before doing shit: datapt=%s; toy_bin_avgs=%s"%(datapt, toy_bin_avgs))
    #now make the things I actually want to plot
    #difference in pol4_fit and toy avg; diff in pol6_fit and toy avg; diff in data and toy_avg
    p4_diff = [0 for i in range(nbins)]
    p6_diff = [0 for i in range(nbins)]
    datdiff = [0 for i in range(nbins)]
    #toydiff will always be 0
    toydiff = [0 for i in range(nbins)]
    for i in range(nbins):
        toy_bin_avgs[i] = toy_bin_avgs[i] - z_vals[i] - t_vals[i]
        p4_diff[i] = p4_dat[i] - toy_bin_avgs[i]
        p6_diff[i] = p6_dat[i] - toy_bin_avgs[i] 
        datdiff[i] = datapt[i] - z_vals[i] - t_vals[i] - toy_bin_avgs[i]
        toydiff[i] = toy_bin_avgs[i] - toy_bin_avgs[i]

    #convert to arrays for easy plotting 
    x_varr = array("f", x_vals)
    p4_arr = array("f", p4_diff)
    p4earr = array("f", p4_err)
    p6_arr = array("f", p6_diff)
    p6earr = array("f", p6_err)
    datarr = array("f", datdiff)
    #error array for data
    dtearr = array("f", daterr)
    toyarr = array("f", toydiff)
    #error array for toys
    tyearr = array("f", toy_bin_errs)
    #xerr is 0 but is a required argument for the TGraph constructor for some stupid ass reason
    xerr = array('f', [0 for i in range(nbins)])

    #now plot these bois (with error of course)
    can2 = TCanvas("cat%d_toy_diff_%s"%(my_cat, spec_name),"cat%d_toy_diff_%s"%(my_cat, spec_name),600,600)
    #pol4 graph
    #draw cms official text at top
    can2.cd(0)
   # pCMS1.Draw()
   # pCMS2.Draw()
#    can2.SetTitle("Pol3,4,5,6")
    can2.SetFillColor(0)
    can2.SetFillStyle(-1)
    gPad.SetBottomMargin(0.3)
    print("x_varr: %s \n p4_arr: %s \n p4earr:%s \n lengths: %d, %d, %d"%(x_varr, p4_arr, p4earr, len(x_varr), len(p4_arr), len(p4earr)))
    print("x_varr: %s \n p6_arr: %s \n p6earr:%s \n lengths: %d, %d, %d"%(x_varr, p6_arr, p6earr, len(x_varr), len(p6_arr), len(p6earr)))
    print("x_varr: %s \n datarr: %s \n dtearr:%s \n lengths: %d, %d, %d"%(x_varr, datarr, dtearr, len(x_varr), len(datarr), len(dtearr)))
    grp4 = TGraphErrors(nbins, x_varr, p4_arr, xerr, p4earr)
    #pol6 graph
    grp6 = TGraphErrors(nbins, x_varr, p6_arr, xerr, p6earr)
    #data graph
    grdt = TGraphErrors(nbins, x_varr, datarr, xerr, dtearr)
    #toy graph
    #grty = TGraphErrors(nbins, x_varr, toyarr, xerr, tyearr)
    grty = TGraph(nbins, x_varr, toyarr)
    
    #width of each bin
#    #ndf = int((opts.X[1]-opts.X[0])/(opts.dX[0])) - (n_param)
#    ndf = int((xmax-xmin)/dX) - n_param
#    chi2=chi2_val/ndf
#    prob = ROOT.TMath.Prob(chi2_val,ndf)
# part 2
    pad = TPad("pad","pad",0.,0.,1.,1.)
#    pad.SetTopMargin(0.7)
    pad.SetFillStyle(-1)
    pad.SetBorderSize(0)
    pad.Draw()
    pad.cd(0)
#    i = n_param
    grdt.SetMarkerStyle(7)
#    grdt.GetYaxis().SetRangeUser(41000, 45000) 
#    grdt.GetXaxis().SetRangeUser(120, 130) 
    grdt.GetXaxis().SetTitleSize(.055)
    grdt.GetXaxis().SetLabelSize(.035)
#zoom in on the beginning of the x-axis!
#    grdt.GetXaxis().SetRangeUser(80, 100)
#    grdt.GetYaxis().SetRangeUser(50000, 60000)

    grdt.GetYaxis().SetLabelSize(0.035)
    grdt.GetYaxis().SetTitleSize(0.055)
    grdt.GetYaxis().SetTitleOffset(1.5);
    #grdt.GetYaxis().SetTitle("Events / %.1f GeV"%(opts.dX[0]*10))
    grdt.GetYaxis().SetTitle("Events / %.1f GeV"%(dX))
    grdt.GetXaxis().SetTitle("M_{bb} (GeV)")
    grdt.Draw("same APE0")
    #grdt.SetLineColor(kBlack)
    #draw pol4 with errors
    grp4.Draw("same ce4")
    grp4.SetLineColor(kBlue)
    grp4.SetFillColor(kBlue)
    grp4.SetFillStyle(3005)
    #draw pol6 with errors
    grp6.Draw("same ce4")
    grp6.SetLineColor(kRed)
    grp6.SetFillColor(kRed)
    grp6.SetFillStyle(3005)
    #draw toy diff with errors
    grty.Draw("same ce4")
    grty.SetLineColor(kGreen)
    grty.SetFillColor(kGreen)
    grty.SetFillStyle(3005)


    #L00 = TLegend(0.7,0.6,1.-gStyle.GetPadRightMargin()-gStyle.GetPadTopMargin()*0.3333,1.-gStyle.GetPadTopMargin()*1.3333)
    L00 = TLegend(0.65,0.85,1.0,1.0)
    L00.SetFillColor(0)
    L00.SetLineColor(1)
    L00.SetBorderSize(10)
    L00.SetTextFont(42)
    L00.SetTextSize(gStyle.GetPadTopMargin()*0.50)
#    L00.SetY1(L00.GetY2()-L00.GetNRows()*gStyle.GetPadTopMargin()*0.95)
#    L00.AddEntry(grp4, "Pol4 - toy avg", 'l')
#    L00.AddEntry(grp6, "Pol6 - toy avg", 'l')
#    L00.AddEntry(grdt, "data-Z-top-toy avg", 'p')
#    L00.AddEntry(grty, "toy avg - toy avg", 'l')
    #L00.AddEntry(grp4, "Pol4-Z-top", 'l')
    #L00.AddEntry(grp6, "Pol6-Z-top", 'l')
    L00.AddEntry(grp4, "qcd_model_Pol4", 'l')
    L00.AddEntry(grp6, "qcd_model_Pol6", 'l')
    L00.AddEntry(grdt, "data-Z-top", 'p')
    L00.AddEntry(grty, "toy avg-Z-top", 'l')
#    L00.AddEntry(dat_hist, "data", 'p')
#    dat_hist.SetMarkerStyle(7)
#    dat_hist.SetMarkerColor(kGreen)
    #draw raw data (without subtracting other background)
    #dat_hist.Draw("same pe")
    L00.Draw("same")
    print("drew legend")
    grty.SetTitle("CAT%d"%my_cat)
    #can2.Draw("same")
    gPad.Update()
    gPad.Modified()
    #now save canvas to file so can view.
    can2.SaveAs("%s/plot/draw_funcs/%s.pdf"%(workdir,can2.GetName()))
    can2.SaveAs("%s/plot/draw_funcs/%s.png"%(workdir,can2.GetName()))

    print("saved my boi")

####################################################################################################
def main():
    my_cat = 0
    #special name of the directory to find the toy shapes
    spec_name = "noSys_shapes"
   # spec_name = "expectSig0" 
    #plot_diffs(my_cat)
    plot_toydiff(my_cat, spec_name)

if __name__=='__main__':
    main()
