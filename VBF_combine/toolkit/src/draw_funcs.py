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

####################################################################################################
def parser(mp=None):
    if not mp: mp = OptionParser()
#
    mp.add_option('--workdir',help=colours[5]+'Case workdir.'+colours[0],default='case0',type='str')
    mp.add_option('--long',help=colours[5]+'Long filenames.'+colours[0],default=False,action='store_true')
    mp.add_option('-v','--verbosity',help=colours[5]+'Verbosity.'+colours[0],default=0,action='count')
    mp.add_option('-q','--quiet',help=colours[5]+'Quiet.'+colours[0],default=True,action='store_false')
    mg1 = OptionGroup(mp,'Selection setup')
    mg1.add_option('--SELCATs',help=colours[5]+'Selection/Category setup.'+colours[0],default='double;DoubleB;-1#0.0875#0.5775#0.7875#1.,single;SingleB;-1#0.0775#0.5775#0.7775#0.8775#1.',type='str',action='callback',callback=SELsetup,dest='SC')# final
    mp.add_option_group(mg1)
#
    mg2 = OptionGroup(mp,'Data template setup')
    mg2.add_option('--bounds',help=colours[5]+'Template boundaries: 80,200 (min,max)'+colours[0],default=[80.,200.],type='str',action='callback',callback=optsplitfloat,dest='X')
    mg2.add_option('--binwidth',help=colours[5]+'Template bin width: 1.,1. (Single,Double)'+colours[0],default=[1.,1.],type='str',action='callback',callback=optsplitfloat,dest='dX')
    mg2.add_option('--lumi',help=colours[5]+'Luminosity: 35900(single,double)'+colours[0],default=[35900,35900.],type='str',action='callback',callback=optsplitfloat,dest='LUMI')
    mg2.add_option('--TF',help=colours[5]+'Transfer function label: POL1,POL1 (Single,Double)'+colours[0],default=['POL2','POL4'],type='str',action='callback',callback=optsplit)
    mg2.add_option('--function',help=colours[5]+'Function for bias study :  expPow'+colours[0],default='expPow',type='str')
    mg2.add_option('--forfit',help=colours[5]+'Freeze/free parameters for fit or not (for toys)'+colours[0],default=True,action='store_false')
    mp.add_option_group(mg2)
#
    return mp

####################################################################################################
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

def hStyle(h,i,sumw2=True):
    if sumw2: h.Sumw2()
    h.SetMarkerStyle([20,20,23,21,22,24][i])
    h.SetMarkerSize(1.2)
    h.SetMarkerColor([kBlack,kBlue,kRed,kGreen+2,kOrange][i])
    h.SetLineColor([kBlack,kBlue,kRed,kGreen+2,kOrange][i])

def gpStyle(g,i):
    g.SetFillColor([kBlack,kBlue,kRed,kGreen+2,kOrange][i])
    g.SetFillStyle(3004)

def gaStyle(g,i):
    g.SetFillColor(kGray)
    g.SetFillStyle(1001)

def main():
# Parse options
    mp = parser()
    opts,args = mp.parse_args()

# Style
    style()

# Create directories if needed
    makeDirs('%s'%opts.workdir)
    makeDirs('%s/plot'%opts.workdir)
    makeDirs('%s/root'%opts.workdir)
#    longtag1 = "_%s"%(func)
    
# Setup
    SC = opts.SC if not type(opts.SC)==str else SELsetup(opts.SC)
#    fTF = json.loads(filecontent("%s/vbfHbb_transfer_2016_run2_linear.json"%basepath))
    LUMI = opts.LUMI
    NBINS = [int((opts.X[1]-opts.X[0])/(x)) for x in opts.dX] 
    print("NBINS: %s"%NBINS)
    MASS=125

# Files
    fBKG = TFile.Open("%s/root/bkg_shapes_workspace%s.root"%(opts.workdir,"" if not opts.long else longtag1),"read")
    fSIG = TFile.Open("%s/root/sig_shapes_workspace%s.root"%(opts.workdir,"" if not opts.long else longtag1),"read")
    fTRF = TFile.Open("%s/root/TransferFunctions%s.root"%(opts.workdir,"" if not opts.long else longtag1),"read")
    wBKG = fBKG.Get("w")
    wSIG = fSIG.Get("w")
    w    = RooWorkspace("w","workspace")

# New histograms, graphs, containers
    archive = []
    C0  = TCanvas("c1","c1",900,750)
    rh  = {}
    rhb = {}
    h   = {}
    hb  = {}
    Y   = {}
## CMS info
    left,right,top,bottom = gStyle.GetPadLeftMargin(),gStyle.GetPadRightMargin(),gStyle.GetPadTopMargin(),gStyle.GetPadBottomMargin()
    pCMS1 = TPaveText(left,1.-top,0.4,1.,"NDC")
    pCMS1.SetTextFont(62)
    pCMS1.SetTextSize(top*0.75)
    pCMS1.SetTextAlign(12)
    pCMS1.SetFillStyle(-1)
    pCMS1.SetBorderSize(0)
    pCMS1.AddText("CMS")
    pCMS2 = TPaveText(0.5,1.-top,1.-right*0.5,1.,"NDC")
    pCMS2.SetTextFont(42)
    pCMS2.SetTextSize(top*0.75)
    pCMS2.SetTextAlign(32)
    pCMS2.SetFillStyle(-1)
    pCMS2.SetBorderSize(0)
    pCMS2.AddText("L = 35.9 fb^{-1} (13 TeV)")
    c0 = TCanvas("can","can",600,600)

### WHICH CATEGORY TO ANALYZE
    my_cat = 7
# Selection loop
#    for iS,S in enumerate(SC.selections):
## Load tree
    if my_cat < 4:
        fin = TFile.Open("/afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/VBFHbb2016/inputs/FitVBF_BTagCSV_analysis_double_trignone_v25_VBF_newReg.root","read")
    else:
        fin = TFile.Open("/afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/VBFHbb2016/inputs/FitVBF_BTagCSV_analysis_single_trignone_v25_VBF_newReg.root","read")
        
    T = fin.Get("VBF/events")
## Containers
    brn          = {}
    params   = {}
    qcd_pdf      = {}
    qcd_pdf_aux  = {}
    qcd_pdf_aux2 = {}
    model        = {}
    nQCD         = {}
    x            = {}
    x_name       = {}
    res = {}
#    transfer     = {}
    zPDF,tPDF,qPDF,sPDF = {},{},{},{}
    yZ,yT,yQ     =  {} ,{},{}
    can2 = TCanvas("cat%d_funcs_diff"%(my_cat),"cat%d_funcs_diff"%(my_cat),600,600)
    iS = 0
    S = SC.selections[0]
    if my_cat > 3:
        S = SC.selections[1]
    for order in range(3, 7):
        N = "Pol%d"%order       
        n_param = order 
        func = "Pol%d"%order
        print("starting func %s"%func)
## Category loop
        #for C in range(S.ncat):
        for C in range(my_cat, my_cat+1):
            Cp = C%4
            if C == 8:
                Cp = 4
#### Start of RooFit part 
  ## x
            for i in range(3, order):
                poln = "Pol%d"%i
                chi2 = RooChi2Var("chi2", "chi2", model[poln], rh[poln])
                #print("beginning of cat loop order %d, chi2 value for %s is %f"%(order, poln, chi2.getVal()))
            x[N] = RooRealVar("mbbReg_CAT%d"%C,"mbbReg_CAT%d"%C,opts.X[0],opts.X[1])
            x_name[N]="mbbReg_CAT%d"%C
            trans_p = {}

  ### QCD part
      ### Containers
            h[N]   = TH1F("hMbb_%s"%N,"hMbb_%s"%N,NBINS[iS],opts.X[0],opts.X[1])
            hb[N]  = TH1F("hMbb_blind_%s"%N,"hMbb_blind_%s"%N,NBINS[iS],opts.X[0],opts.X[1])
      ### Define cut
            cut    = TCut("(bdt_VBF>%1.4f && bdt_VBF<=%1.4f)"%(S.boundaries[Cp],S.boundaries[Cp+1]))
            cutB   = TCut("(bdt_VBF>%1.4f && bdt_VBF<=%1.4f) && mbbRegFSR>100 && mbbRegFSR<150"%(S.boundaries[Cp],S.boundaries[Cp+1]))
            print("C=%d, Cp=%d, cut=%s, cutB=%s"%(C, Cp, cut, cutB))
            print("Cat = %d, Cp = %d, cut = %s, cutB = %s"%(C, Cp, cut, cutB))
      ### Draw
            c0.cd()
                        #c0.SetLogy()
            T.Draw("mbbRegFSR>>%s"%(h[N].GetName()),cut)
            c0.SaveAs(N+"test.png")
            T.Draw("mbbRegFSR>>%s"%(hb[N].GetName()),cutB)
         #   c0.SaveAs(N+"Btest.png")
      ### Yields
            Y[N]   = RooRealVar("yield_data_CAT%d"%C,"yield_data_CAT%d"%C,h[N].Integral())
            #print("h[n] integral: " + str(h[N].Integral()))
      ### Histograms
            rh[N]  = RooDataHist("data_hist_CAT%d"%C,"data_hist_CAT%d"%C,RooArgList(x[N]),h[N])
            rhb[N] = RooDataHist("data_hist_blind_CAT%d"%C,"data_hist_blind_CAT%d"%C,RooArgList(x[N]),hb[N])
#            h_data = rh[N].createHistogram('h_data',x[N])
  ### Model
            for i in range(3, order):
                poln = "Pol%d"%i
                chi2 = RooChi2Var("chi2", "chi2", model[poln], rh[poln])
                #print("middle0 of cat loop order %d, chi2 value for %s is %f"%(order, poln, chi2.getVal()))
#                if Cp==4 :  gcs_aux[:]=[]
            gcs_aux[N] = []
          #  val140 = h[N].GetBinContent(NBINS[iS]/2)
            for i in range(3, order):
                poln = "Pol%d"%i
                chi2 = RooChi2Var("chi2", "chi2", model[poln], rh[poln])
                #print("middle1/4 of cat loop order %d, chi2 value for %s is %f"%(order, poln, chi2.getVal()))
            [qcd_pdf[N],params[N]] = generate_pdf(x[N], pdf_name=func,x_name=x_name[N],selection=double,gcs=gcs_aux[N],real_cat=my_cat) #, val140=val140) 
     #       #print("Just after call to generate_pdf, gcs_aux: %s"%gcs_aux)
            if qcd_pdf[N] == None : 
                sys.exit("qcd_pdf = None !!!!!!!!!!!!!!!!!")
    
            for i in range(3, order):
                poln = "Pol%d"%i
                chi2 = RooChi2Var("chi2", "chi2", model[poln], rh[poln])
                #print("middle1/2 of cat loop order %d, chi2 value for %s is %f"%(order, poln, chi2.getVal()))
      ### PDFs
            zPDF[N] = wBKG.pdf("Z_model_CAT%d"%C)
            tPDF[N] = wBKG.pdf("Top_model_CAT%d"%C)
            sPDF[N] = wSIG.pdf("signal_model_mH%d_sel%s_CAT%d"%(MASS,S.tag,C))

      ### Yields
            yZ[N]    = wBKG.var("yield_ZJets_CAT%d"%C)
            yT[N]    = wBKG.var("yield_Top_CAT%d"%C)
            yq_est = h[N].Integral()
            yQ[N]    = RooRealVar("yield_QCD_CAT%d"%C,"yield_QCD_CAT%d"%C,yq_est,yq_est/2.,2*yq_est)
            yZ[N].setConstant(kTRUE)
            yT[N].setConstant(kTRUE)

            for i in range(3, order):
                poln = "Pol%d"%i
                chi2 = RooChi2Var("chi2", "chi2", model[poln], rh[poln])
                #print("middle1 of cat loop order %d, chi2 value for %s is %f"%(order, poln, chi2.getVal()))

  ### Combined model
#            model[N] = RooAddPdf("bkg_model_CAT%d"%(Cp),"bkg_model_CAT%d"%(Cp),RooArgList(zPDF[N],tPDF[N],qcd_pdf[N]),RooArgList(yZ[N],yT[N],yQ[N]))
            yQ[N].setConstant(kTRUE)
            #print("before: yZ=%s yT=%s yQ=%s"%(yZ[N], yT[N], yQ[N]))
            #print 'before: Yields Z,Top,QCD: ', yZ[N].getVal(),yT[N].getVal(),yQ[N].getVal()
            model[N] = RooAddPdf("bkg_model_%s_CAT%d"%(opts.TF[iS],C),"bkg_model_%s_CAT%d"%(opts.TF[iS],C),RooArgList(zPDF[N],tPDF[N],qcd_pdf[N]),RooArgList(yZ[N],yT[N],yQ[N])) 
  ### Fit
      #      print("Just before call to fitTo, gcs_aux: %s"%gcs_aux)
            res[N]   = model[N].fitTo(rh[N],RooFit.Save(),RooFit.Warnings(ROOT.kTRUE))
            #print 'kTrue: Yields Z,Top,QCD: ', yZ[N].getVal(),yT[N].getVal(),yQ[N].getVal()
            yQ[N].setConstant(kFALSE)
            res[N]   = model[N].fitTo(rh[N],RooFit.Save(),RooFit.Warnings(ROOT.kTRUE))
            #print("kfalse:yZ=%s yT=%s yQ=%s"%(yZ[N], yT[N], yQ[N]))
            #print 'kFalse:Yields Z,Top,QCD: ', yZ[N].getVal(),yT[N].getVal(),yQ[N].getVal()

            #print("model[%s] = %s" %(N, model[N]))
            for i in range(3, order+1):
                poln = "Pol%d"%i
                chi2 = RooChi2Var("chi2", "chi2", model[poln], rh[poln])
                #print("end of cat loop order %d, chi2 value for %s is %f"%(order, poln, chi2.getVal()))
            chi2 = RooChi2Var("chi2", "chi2", model[N], rh[N])
            chi2_val = chi2.getVal()
            #print 'Yields Z,Top,QCD: ', yZ[N].getVal(),yT[N].getVal(),yQ[N].getVal()
            total_fitted_yield = yZ[N].getVal()+yT[N].getVal()+yQ[N].getVal()
            h_data = rh[N].createHistogram('h_data',x[N])
            #print "h_data_integral=%.3f"%(h_data.Integral())
            ks=0
            #print 'KS probability = ',ks

           # chi2 = RooChi2Var("chi2", "chi2", model[N], rh[N])
           # chi2_val = chi2.getVal()
           # print("chi2 of the fit for %s: %f"%(N, chi2_val))
#            rh[N]  = RooDataHist("data_hist_CAT%d"%C,"data_hist_CAT%d"%C,RooArgList(x[N]),h[N].Rebin(10))            
###--- end of CAT loop
    #print("models: %s"%model)
    #how many bins to use for the fitted histogram
    nbins = 120
    #now plot the four curves.
    pCMS1.Draw()
    pCMS2.Draw()
    can2.cd(0)
    can2.SetTitle("Pol3,4,5,6")
    can2.SetFillColor(0)
    can2.SetFillStyle(-1)
    gPad.SetBottomMargin(0.3)
# part 1
    for i in range(3, 7):
        poln = "Pol%d"%i
        chi2 = RooChi2Var("chi2", "chi2", model[poln], rh[poln])
        #print("AFTER CAT LOOP, CHI2 VALUE FOR %s IS %f"%(poln, chi2.getVal()))
    #need to get error for each of the fits, also plot just the qcd and not the other backgrounds.
    frametop = x["Pol3"].frame()
#    framebot = x["Pol3"].frame()
#    rh["Pol3"].plotOn(frametop)
#    model["Pol3"].plotOn(frametop,RooFit.LineWidth(2), RooFit.LineColor(2)) #red
##    rh["Pol4"].plotOn(frametop)
#    model["Pol4"].plotOn(frametop,RooFit.LineWidth(2), RooFit.LineColor(3)) #green
##    rh["Pol5"].plotOn(frametop)
#    model["Pol5"].plotOn(frametop,RooFit.LineWidth(2), RooFit.LineColor(4)) #blue  
##    rh["Pol6"].plotOn(frametop)
#    model["Pol6"].plotOn(frametop,RooFit.LineWidth(2), RooFit.LineColor(5)) #yellow
   # print("model after plotting:")
    model["Pol6"].Print()
    #print("fit result:")
    res["Pol6"].Print()

    #histogram of the Pol4 fit
    pol4_fitHist = model["Pol4"].createHistogram(x["Pol4"].GetName(), nbins)
    #histogram of the Pol6 fit
    pol6_fitHist = model["Pol6"].createHistogram(x["Pol6"].GetName(), nbins)
    #histogram of the Z background data
    Z_hist = zPDF["Pol4"].createHistogram(x["Pol4"].GetName(), nbins)
    #histogram of the top background data
    top_hist = tPDF["Pol6"].createHistogram(x["Pol6"].GetName(), nbins)
    #don't worry, x["Pol4"] and x["Pol6"] have the same name (mbbreg...CAT[my_cat] or whatever)
    #histogram of the data (just need to convert from RooDataHist to TH1F)
    dat_hist = rh["Pol4"].createHistogram(x["Pol4"].GetName(), nbins)
    pol4_fitHist.Print()
    pol6_fitHist.Print()
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
    for i in range(1, nbins+1):
        x_vals[i-1] = pol4_fitHist.GetBinCenter(i)
        pol4_val = pol4_fitHist.GetBinContent(i)
        pol6_val = pol6_fitHist.GetBinContent(i)
        z_val = Z_hist.GetBinContent(i)
        top_val = top_hist.GetBinContent(i)
        data_val = dat_hist.GetBinContent(i)
        #subtract Z and top to get actual values for qcd.
        print("boutta set datapt: data_val=%f, z_val=%f, top_val=%f"%(data_val, z_val, top_val))
        p6_dat[i-1] = pol6_val - z_val*yZ["Pol3"].getVal() - top_val*yT["Pol3"].getVal()
        #subtract off the pol6 value for this plot.
        datapt[i-1] = data_val - z_val*yZ["Pol3"].getVal() - top_val*yT["Pol3"].getVal() - p6_dat[i-1]
        p4_dat[i-1] = pol4_val - z_val*yZ["Pol3"].getVal() - top_val*yT["Pol3"].getVal() - p6_dat[i-1]
        p6_dat[i-1] = 0
        #also need errors
        daterr[i-1] = dat_hist.GetBinError(i)
        p4_err[i-1] = pol4_fitHist.GetBinError(i)
        p6_err[i-1] = pol6_fitHist.GetBinError(i)

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
    #pol4 graph
    print("x_varr: %s \n p4_arr: %s \n p4earr:%s \n lengths: %d, %d, %d"%(x_varr, p4_arr, p4earr, len(x_varr), len(p4_arr), len(p4earr)))
    print("x_varr: %s \n datarr: %s \n dtearr:%s \n lengths: %d, %d, %d"%(x_varr, datarr, dtearr, len(x_varr), len(p4_arr), len(p4earr)))
    grp4 = TGraphErrors(nbins, x_varr, p4_arr, xerr, p4earr)
    #pol6 graph
    grp6 = TGraphErrors(nbins, x_varr, p6_arr, xerr, p6earr)
    #data graph
    grdt = TGraphErrors(nbins, x_varr, datarr, xerr, dtearr)
    
    
#plot qcd also
#plot normalized by yield
 #   qPDF.plotOn(frametop,RooFit.LineWidth(1), RooFit.LineColor(2)) #red
 #   zPDF.plotOn(frametop,RooFit.LineWidth(1), RooFit.LineColor(3)) #green
 #   tPDF.plotOn(frametop,RooFit.LineWidth(1), RooFit.LineColor(5)) #yellow


    #rhres = frametop.pullHist()
   # p1 = RooArgSet(qPDF)
   # p2 = RooArgSet(zPDF)
   # p3 = RooArgSet(tPDF)
    #model.plotOn(frametop,RooFit.Components(p1),RooFit.LineWidth(2),RooFit.LineColor(kBlack),RooFit.LineStyle(kDashed))
    
    ndf = int((opts.X[1]-opts.X[0])/(opts.dX[0])) - (n_param)
    chi2=chi2_val/ndf
    prob = ROOT.TMath.Prob(chi2_val,ndf)
# part 2
    pad = TPad("pad","pad",0.,0.,1.,1.)
#    pad.SetTopMargin(0.7)
    pad.SetFillStyle(-1)
    pad.SetBorderSize(0)
    pad.Draw()
    pad.cd(0)
#    framebot.addPlotable(rhres,"p")
#    framebot.GetYaxis().SetTitle("Pulls")
#    framebot.GetXaxis().SetTitle("M_{bb} (GeV)")
#    framebot.GetYaxis().SetLabelSize(0.03)
#    framebot.GetYaxis().SetTitleSize(0.04)
#    line = ROOT.TLine(opts.X[0],0,opts.X[1],0)
#    line.SetLineStyle(7)
#    line.SetLineWidth(2)
    i = n_param
#    line.SetLineColor(kRed+i)
#    framebot.addObject(line)

#    framebot.addObject(L00)
#    pave = TPaveText(0.75,0.7,0.9,0.9,"NDC")
#    pave.SetTextAlign(11)
#    pave.SetFillStyle(-1)
#    pave.SetBorderSize(0)
#    pave.SetTextFont(42)
#    pave.SetTextSize(gStyle.GetPadTopMargin()*0.60)
#    for i in range(4, 7, 2):
#        pave.AddText("Pol%d"%i) 
#        pave.GetListOfLines().Last().SetTextColor(i-1);
#    pave.Draw()
    #frametop.addObject(pave)
#    grdt.SetFillColor(6);
    grdt.SetMarkerStyle(7)
#    grdt.SetMarkerSize(2)
    #grdt.SetFillStyle(3005);
#    grdt.GetYaxis().SetRangeUser(41000, 45000) 
#    grdt.GetXaxis().SetRangeUser(120, 130) 
    grdt.GetXaxis().SetTitleSize(.055)
    grdt.GetXaxis().SetLabelSize(.035)
    grdt.GetYaxis().SetLabelSize(0.035)
    grdt.GetYaxis().SetTitleSize(0.055)
    grdt.GetYaxis().SetTitleOffset(1.5);
    grdt.GetYaxis().SetTitle("Events / %.1f GeV"%(opts.dX[0]*10))
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

    frametop.GetXaxis().SetTitleSize(.055)
    frametop.GetXaxis().SetLabelSize(.035)
    frametop.GetYaxis().SetLabelSize(0.035)
    frametop.GetYaxis().SetTitleSize(0.055)
    frametop.GetYaxis().SetTitleOffset(1.5);
    frametop.GetYaxis().SetTitle("Events / %.1f GeV"%(opts.dX[0]*10))
    frametop.GetXaxis().SetTitle("M_{bb} (GeV)")
  #  frametop.Draw("same")
    gPad.Update()
    print("drew my bois")

    #L00 = TLegend(0.7,0.6,1.-gStyle.GetPadRightMargin()-gStyle.GetPadTopMargin()*0.3333,1.-gStyle.GetPadTopMargin()*1.3333)
    L00 = TLegend(0.75,0.85,1.0,1.0)
    L00.SetFillColor(0)
    L00.SetLineColor(1)
    L00.SetBorderSize(10)
    L00.SetTextFont(42)
    L00.SetTextSize(gStyle.GetPadTopMargin()*0.50)
#    L00.SetY1(L00.GetY2()-L00.GetNRows()*gStyle.GetPadTopMargin()*0.95)
   # if "6" in title:
    L00.AddEntry(grp4, "Pol4 - Pol6", 'l')
    L00.AddEntry(grp6, "Pol6 - Pol6", 'l')
    L00.AddEntry(grdt, "data-Z-top-Pol6", 'p')
#    L00.AddEntry(dat_hist, "data", 'p')
    dat_hist.SetMarkerStyle(7)
    dat_hist.SetMarkerColor(kGreen)
    #draw raw data (without subtracting other background)
    #dat_hist.Draw("same pe")
    L00.Draw("same")
    print("drew legend")
    #now save canvas to file so can view.
    makeDirs("%s/plot/draw_funcs/"%opts.workdir)
    can2.SaveAs("%s/plot/draw_funcs/%s.pdf"%(opts.workdir,can2.GetName()))
    can2.SaveAs("%s/plot/draw_funcs/%s.png"%(opts.workdir,can2.GetName()))

    print("saved my boi")
    sys.exit(0)

    #PLOTTING RATIO OF POL6 / POL4 NOW
    canrat = TCanvas("cat%d_Pol6_4_ratio"%my_cat,"cat%d_Pol6_4_ratio"%my_cat,600,600)
    canrat.cd()
    pCMS1.Draw()
    pCMS2.Draw()
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
    x6name = x["Pol6"].GetName()
    fitHist6 = model["Pol6"].createHistogram(x6name, nbins)
    fitHist6.Print()
    x4name = x["Pol4"].GetName()
    fitHist4 = model["Pol4"].createHistogram(x4name, nbins)
    fitHist4.Print()
    ratio64 = fitHist6.Divide(fitHist4) #ratio64 is only a bool for some reason (ratio his is still in fitHist6)
    
    frameboi = x["Pol4"].frame()
    #draw ratio histogram with axes, curve style, with error bars.
    frameboi.GetXaxis().SetTitle("m_{bb} (GeV)")
    frameboi.GetYaxis().SetTitle("Ratio Pol6 / Pol4")
    frameboi.SetAxisRange(0.95, 1.05, "Y")  
    frameboi.Draw()
    fitHist6.Draw("same AC") #add "E" for error bars (they'll prolly be hella big)
    canrat.Modified()
    canrat.Update()

    canrat.SaveAs("%s/plot/draw_funcs/%s.pdf"%(opts.workdir, canrat.GetName()))

    #NOW PLOTTING TOY RATIO THINGY
    #now get toy, and divide the fits by it
    toy_path = "/afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/bias_studies/cat%d/noSysFit/Pol6_Pol4/fitDiagnostics_0_0.root"%my_cat
    #if toy doesn't exist yet then create it.
    if not os.path.exists(toy_path):
        os.system("cd /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/bias_studies/cat%d/noSysFit/Pol6_Pol4;"%(my_cat) \
            + " eval `scramv1 runtime -sh`;" \
            + " combine -M GenerateOnly --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --X-rtd ADDNLL_RECURSIVE=0 --X-rtd FITTER_NEW_CROSSING_ALGO --toysNoSystematics  -t 1 --expectSignal 1.0 --saveToys  --rMin=-100 --rMax=100 -n _0 --seed=679 /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett/datacards/datacard_vbfHbb_biasPol6_m125_CAT%d-CAT%d_CATveto.txt;"%(my_cat, my_cat) \
            + " combine -M FitDiagnostics --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --X-rtd ADDNLL_RECURSIVE=0 --X-rtd FITTER_NEW_CROSSING_ALGO -t 1 --saveShapes --seed=-1 -v 0 --expectSignal=1.0 --robustFit=1 --rMin=-100 --rMax=100 --toysFile higgsCombine_0.GenerateOnly.mH120.679.root -n _0_0 --setRobustFitTolerance 0.01  -S 0  -d /afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/lata_code/VBFHbb2016/VBF_combine/toolkit/src/test_for_bennett/datacards/datacard_vbfHbb_biasPol4_m125_CAT%d-CAT%d_CATveto.txt"%(my_cat, my_cat)   )
    toy = TFile.Open(toy_path, "read")
    qcd = toy.Get("shapes_fit_s/CAT%d/qcd"%my_cat)
    print("obtained qcd:")
    qcd.Print()
    fitHist4.Print()
    #Now will divide Pol4 fitHist by the toy hist.
    #but need to rebin first so that the to hists have the same number of bins.
    #argument for rebin method is ngroups instead of new nbins for some reason.
    #also need to normalize by number of entries in each histogram.
    nent_toy = qcd.GetEntries()
    nent_hist = fitHist4.GetEntries()
    #c1 is additional, normalization argument to multiply first histogram, c2 for second hist.
    c1 = 1.0 / nent_hist
    c2 = 1.0 / nent_toy 
    print("nentries toy=%d, nentries fitHist=%d, c1=%f, c2=%f"%(nent_toy, nent_hist, c1, c2))
    #x values to be plotted
    xes = [0.0 for i in range(nbins)]
    #y values (ratios) to be plotted
    toyfit_ratio = [0.0 for i in range(nbins)]
    #bin 0 is the overflow bin for some reason
    for i in range(1, nbins+1):
        fitHist_entry = fitHist4.GetBinContent(i)
        #it's 10 times more bins in the toy histogram than in the data fitting one.
        toyHist_entry = qcd.GetBinContent((i-1)*10+5)
        if toyHist_entry != 0:
            toyfit_ratio[i-1] = fitHist_entry / toyHist_entry
        else:
            toyfit_ratio[i-1] = 0
            print("WARNING: bin %d has value 0 for toy histogram"%i)
        xes[i-1] = fitHist4.GetBinCenter(i)
        x2 = qcd.GetBinCenter((i-1)*10+5)
        print("bin:%d xval: %f, x2: %f, fit val: %f, toy val: %f, ratio:%f"%(i, xes[i-1], x2, fitHist_entry, toyHist_entry, toyfit_ratio[i-1]))
    #now print to file, exactly as before.
    #print("fitHist4 after division:")
    #fitHist4.Print()
    cantoy = TCanvas("cat%d_toyPol6_fitPol4"%my_cat,"cat%d_toyPol6_fitPol4"%my_cat,600,600)
    cantoy.cd()
    pCMS1.Draw()
    pCMS2.Draw()
#    cantoy.SetTitle("toyPol6_fitPol4_ratio")
#    cantoy.SetTitle("toyPol6_fitPol4_ratio")
    cantoy.SetFillColor(0)
    cantoy.SetFillStyle(-1)
    frametoi = x["Pol4"].frame()
    #draw ratio histogram with axes, curve style, with error bars.
    frametoi.GetXaxis().SetTitle("m_{bb} (GeV)")
    frametoi.GetYaxis().SetTitle("fit(Pol4)/toy(Pol6) [normalized]")
    #frametoi.SetAxisRange(0.9, 1.1, "Y")  
    frametoi.SetAxisRange(0, 30000, "Y")  
    frametoi.Draw()

    #fit_int = fitHist4.Integral()
    #toy_int = qcd.Integral()
    #scale the histograms so their integrals are the same
    #fitHist4.Scale(1.0/fit_int)
    #qcd.Scale(1.0/toy_int)    

    fitHist4.SetLineColor(kBlue)
    qcd.SetLineColor(kRed)
    #fitHist4.Draw("same AC") #add "E" for error bars (they'll prolly be hella big)

    qcd.Draw("same AC")
    fitHist4.Draw("same AC")
    
    toyRatio = fitHist4.Clone("toyRatio")
    toyRatio.SetLineColor(kBlack)
    #convert lists to arrays so they can be plotted
    x = array("f", xes)
    tfr = array("f", toyfit_ratio)
    ratGraph = TGraph(nbins, x, tfr)
    ratGraph.Draw("same AC")
    ratGraph.GetYaxis().SetRangeUser(0.95, 1.05)
    ratGraph.GetXaxis().SetTitle("m_{bb} (GeV)")
    ratGraph.GetYaxis().SetTitle("Ratio fit(Pol4)/toy(Pol6)")
    ratGraph.GetYaxis().SetTitleSize(0.055)
    ratGraph.GetYaxis().SetLabelSize(0.035)
    ratGraph.GetXaxis().SetLabelSize(0.035)
    ratGraph.GetXaxis().SetTitleSize(0.055)
#    toyRatio.Divide(qcd)
#    toyRatio.Draw("same AC")

    cantoy.Modified()
    cantoy.Update()

    cantoy.SaveAs("%s/plot/draw_funcs/%s.pdf"%(opts.workdir, cantoy.GetName()))
#
#--- end of SEL loop
#

#    makeDirs("%s/root/"%opts.workdir)
#    w.writeToFile("%s/root/biasCATS_shapes_workspace_%s.root"%(opts.workdir,opts.function))
#    w.writeToFile("%s/root/bias_shapes_workspace_%s.root"%(opts.workdir,func))
    print "Done."
####################################################################################################
if __name__=='__main__':
    main()
