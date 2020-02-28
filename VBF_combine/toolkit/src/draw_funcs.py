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

####################################################################################################
def RooDraw(opts,can,C,S,x,rh,model,qPDF,zPDF,tPDF,archive,chi2_val,n_param,title,ks, iS): #BG added iS
    can.cd(0)
    can.SetTitle(title)
    can.SetFillColor(0)
    can.SetFillStyle(-1)
    print 'plotting'
    gPad.SetBottomMargin(0.3)
    frametop = x.frame()
    framebot = x.frame()
# part 1
    rh.plotOn(frametop)
    model.plotOn(frametop,RooFit.LineWidth(2)) 
#plot qcd also
#plot normalized by yield
 #   qPDF.plotOn(frametop,RooFit.LineWidth(1), RooFit.LineColor(2)) #red
 #   zPDF.plotOn(frametop,RooFit.LineWidth(1), RooFit.LineColor(3)) #green
 #   tPDF.plotOn(frametop,RooFit.LineWidth(1), RooFit.LineColor(5)) #yellow


    rhres = frametop.pullHist()
    p1 = RooArgSet(qPDF)
    p2 = RooArgSet(zPDF)
    p3 = RooArgSet(tPDF)
    model.plotOn(frametop,RooFit.Components(p1),RooFit.LineWidth(2),RooFit.LineColor(kBlack),RooFit.LineStyle(kDashed))
    frametop.GetXaxis().SetTitleSize(0)
    frametop.GetXaxis().SetLabelSize(0)
    frametop.GetYaxis().SetLabelSize(0.035)
    frametop.GetYaxis().SetTitleSize(0.055)
    frametop.GetYaxis().SetTitleOffset(1.5);
    frametop.GetYaxis().SetTitle("Events / %.1f GeV"%(opts.dX[0]*10))
   # if "6" in title:
    frametop.Draw()
    gPad.Update()
    
    ndf = int((opts.X[1]-opts.X[0])/(opts.dX[0])) - (n_param)
    chi2=chi2_val/ndf
    prob = ROOT.TMath.Prob(chi2_val,ndf)
# part 2
    pad = TPad("pad","pad",0.,0.,1.,1.)
    pad.SetTopMargin(0.7)
    pad.SetFillStyle(-1)
    pad.SetBorderSize(0)
   # if "6" in title:
    pad.Draw()
    pad.cd(0)
    framebot.addPlotable(rhres,"p")
    framebot.GetYaxis().SetTitle("Pulls")
    framebot.GetXaxis().SetTitle("M_{bb} (GeV)")
    framebot.GetYaxis().SetLabelSize(0.03)
    framebot.GetYaxis().SetTitleSize(0.04)
    line = ROOT.TLine(opts.X[0],0,opts.X[1],0)
    line.SetLineStyle(7)
    line.SetLineWidth(2)
    i = n_param
    line.SetLineColor(kRed+i)
    framebot.addObject(line)

    L00 = TLegend(0.7,0.6,1.-gStyle.GetPadRightMargin()-gStyle.GetPadTopMargin()*0.3333,1.-gStyle.GetPadTopMargin()*1.3333)
    L00.SetFillColor(-1)
    L00.SetBorderSize(0)
    L00.SetTextFont(42)
    L00.SetTextSize(gStyle.GetPadTopMargin()*0.60)
    L00.SetY1(L00.GetY2()-L00.GetNRows()*gStyle.GetPadTopMargin()*0.95)
   # if "6" in title:
    L00.Draw()
    framebot.addObject(L00)


####################################################################################################
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

# Selection loop
#    for iS,S in enumerate(SC.selections):
## Load tree
    fin = TFile.Open("/afs/cern.ch/user/b/bgreenbe/private/CMSSW_8_1_0/src/HiggsAnalysis/VBFHbb2016/inputs/FitVBF_BTagCSV_analysis_double_trignone_v25_VBF_newReg.root","read")
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
    can2 = TCanvas("can2_funcs","can2_funcs",600,600)
    iS = 0
    S = SC.selections[0]
    for order in range(3, 7):
        N = "Pol%d"%order       
        n_param = order 
        func = "Pol%d"%order
        print("starting func %s"%func)
## Category loop
        #for C in range(S.ncat):
        for C in range(2, 3):
            Cp = 2
#### Start of RooFit part 
  ## x
            for i in range(3, order):
                poln = "Pol%d"%i
                chi2 = RooChi2Var("chi2", "chi2", model[poln], rh[poln])
                print("beginning of cat loop order %d, chi2 value for %s is %f"%(order, poln, chi2.getVal()))
            x[N] = RooRealVar("mbbReg_CAT%d"%Cp,"mbbReg_CAT%d"%Cp,opts.X[0],opts.X[1])
            x_name[N]="mbbReg_CAT%d"%Cp
            trans_p = {}

  ### QCD part
      ### Containers
            h[N]   = TH1F("hMbb_%s"%N,"hMbb_%s"%N,NBINS[iS],opts.X[0],opts.X[1])
            hb[N]  = TH1F("hMbb_blind_%s"%N,"hMbb_blind_%s"%N,NBINS[iS],opts.X[0],opts.X[1])
      ### Define cut
            cut    = TCut("(bdt_VBF>%1.4f && bdt_VBF<=%1.4f)"%(S.boundaries[C],S.boundaries[C+1]))
            cutB   = TCut("(bdt_VBF>%1.4f && bdt_VBF<=%1.4f) && mbbRegFSR>100 && mbbRegFSR<150"%(S.boundaries[C],S.boundaries[C+1]))
      ### Draw
            c0.cd()
                        #c0.SetLogy()
            T.Draw("mbbRegFSR>>%s"%(h[N].GetName()),cut)
            c0.SaveAs(N+"test.png")
            T.Draw("mbbRegFSR>>%s"%(hb[N].GetName()),cutB)
         #   c0.SaveAs(N+"Btest.png")
      ### Yields
            Y[N]   = RooRealVar("yield_data_CAT%d"%Cp,"yield_data_CAT%d"%Cp,h[N].Integral())
            print("h[n] integral: " + str(h[N].Integral()))
      ### Histograms
            rh[N]  = RooDataHist("data_hist_CAT%d"%Cp,"data_hist_CAT%d"%Cp,RooArgList(x[N]),h[N])
            rhb[N] = RooDataHist("data_hist_blind_CAT%d"%Cp,"data_hist_blind_CAT%d"%Cp,RooArgList(x[N]),hb[N])
            h_data = rh[N].createHistogram('h_data',x[N])
  ### Model
            for i in range(3, order):
                poln = "Pol%d"%i
                chi2 = RooChi2Var("chi2", "chi2", model[poln], rh[poln])
                print("middle0 of cat loop order %d, chi2 value for %s is %f"%(order, poln, chi2.getVal()))
#                if Cp==4 :  gcs_aux[:]=[]
            gcs_aux[N] = []
          #  val140 = h[N].GetBinContent(NBINS[iS]/2)
            for i in range(3, order):
                poln = "Pol%d"%i
                chi2 = RooChi2Var("chi2", "chi2", model[poln], rh[poln])
                print("middle1/4 of cat loop order %d, chi2 value for %s is %f"%(order, poln, chi2.getVal()))
            [qcd_pdf[N],params[N]] = generate_pdf(x[N], pdf_name=func,x_name=x_name[N],selection=double,gcs=gcs_aux[N],real_cat=2) #, val140=val140) 
     #       print("Just after call to generate_pdf, gcs_aux: %s"%gcs_aux)
            if qcd_pdf[N] == None : 
                sys.exit("qcd_pdf = None !!!!!!!!!!!!!!!!!")
    
            for i in range(3, order):
                poln = "Pol%d"%i
                chi2 = RooChi2Var("chi2", "chi2", model[poln], rh[poln])
                print("middle1/2 of cat loop order %d, chi2 value for %s is %f"%(order, poln, chi2.getVal()))
      ### PDFs
            zPDF[N] = wBKG.pdf("Z_model_CAT%d"%Cp)
            tPDF[N] = wBKG.pdf("Top_model_CAT%d"%Cp)
            sPDF[N] = wSIG.pdf("signal_model_mH%d_sel%s_CAT%d"%(MASS,S.tag,Cp))

      ### Yields
            yZ[N]    = wBKG.var("yield_ZJets_CAT%d"%Cp)
            yT[N]    = wBKG.var("yield_Top_CAT%d"%Cp)
            yq_est = h[N].Integral()
            yQ[N]    = RooRealVar("yield_QCD_CAT%d"%Cp,"yield_QCD_CAT%d"%Cp,yq_est,yq_est/2.,2*yq_est)
            yZ[N].setConstant(kTRUE)
            yT[N].setConstant(kTRUE)

            for i in range(3, order):
                poln = "Pol%d"%i
                chi2 = RooChi2Var("chi2", "chi2", model[poln], rh[poln])
                print("middle1 of cat loop order %d, chi2 value for %s is %f"%(order, poln, chi2.getVal()))

  ### Combined model
#            model[N] = RooAddPdf("bkg_model_CAT%d"%(Cp),"bkg_model_CAT%d"%(Cp),RooArgList(zPDF[N],tPDF[N],qcd_pdf[N]),RooArgList(yZ[N],yT[N],yQ[N]))
            yQ[N].setConstant(kTRUE)
            #print("before: yZ=%s yT=%s yQ=%s"%(yZ[N], yT[N], yQ[N]))
            print 'before: Yields Z,Top,QCD: ', yZ[N].getVal(),yT[N].getVal(),yQ[N].getVal()
            model[N] = RooAddPdf("bkg_model_%s_CAT%d"%(opts.TF[iS],Cp),"bkg_model_%s_CAT%d"%(opts.TF[iS],Cp),RooArgList(zPDF[N],tPDF[N],qcd_pdf[N]),RooArgList(yZ[N],yT[N],yQ[N])) 
  ### Fit
      #      print("Just before call to fitTo, gcs_aux: %s"%gcs_aux)
            res[N]   = model[N].fitTo(rh[N],RooFit.Save(),RooFit.Warnings(ROOT.kTRUE))
            print 'kTrue: Yields Z,Top,QCD: ', yZ[N].getVal(),yT[N].getVal(),yQ[N].getVal()
            yQ[N].setConstant(kFALSE)
            res[N]   = model[N].fitTo(rh[N],RooFit.Save(),RooFit.Warnings(ROOT.kTRUE))
            #print("kfalse:yZ=%s yT=%s yQ=%s"%(yZ[N], yT[N], yQ[N]))
            print 'kFalse:Yields Z,Top,QCD: ', yZ[N].getVal(),yT[N].getVal(),yQ[N].getVal()

            print("model[%s] = %s" %(N, model[N]))
            for i in range(3, order+1):
                poln = "Pol%d"%i
                chi2 = RooChi2Var("chi2", "chi2", model[poln], rh[poln])
                print("end of cat loop order %d, chi2 value for %s is %f"%(order, poln, chi2.getVal()))
            #chi2 = RooChi2Var("chi2", "chi2", model[N], rh[N])
            #chi2_val = chi2.getVal()
#            print("chi2 of the fit for %s: %f"%(N, chi2_val))
            chi2 = RooChi2Var("chi2", "chi2", model[N], rh[N])
            chi2_val = chi2.getVal()
#            print("chi2 of the fit for %s: %f"%(N, chi2_val))
            print 'Yields Z,Top,QCD: ', yZ[N].getVal(),yT[N].getVal(),yQ[N].getVal()
            total_fitted_yield = yZ[N].getVal()+yT[N].getVal()+yQ[N].getVal()
            h_data = rh[N].createHistogram('h_data',x[N])
            print "h_data_integral=%.3f"%(h_data.Integral())
            ks=0
            print 'KS probability = ',ks

  ### Draw
           # chi2 = RooChi2Var("chi2", "chi2", model[N], rh[N])
           # chi2_val = chi2.getVal()
           # print("chi2 of the fit for %s: %f"%(N, chi2_val))
            rh[N]  = RooDataHist("data_hist_CAT%d"%Cp,"data_hist_CAT%d"%Cp,RooArgList(x[N]),h[N].Rebin(10))            
            #find chi2 again to make sure it didn't change 
#            chi2 = RooChi2Var("chi2", "chi2", model[N], rh[N])
#            chi2_val = chi2.getVal()
#            print("chi2 of the fit for %s: %f"%(N, chi2_val))
            #RooDraw(opts,can2,C,S,x,rh[N],model[N],qcd_pdf[N],zPDF[N],tPDF[N],archive,chi2_val,n_param,func,ks,iS)
#            model[N].plotOn(frametop,RooFit.LineWidth(2), RooFit.LineColor(order-1)) #red
 #           can2.cd(0)

  ### Reset constant settings
           # for ib in range(n_param): 
           #     if C==0 :
           #         gcs_aux[ib].setConstant(kFALSE if opts.forfit else kTRUE)  ## false for final fit, true if fixed)
  ### Saving
#            for o in [rh[N],rhb[N],model[N],Y[N]]:
#                getattr(w,'import')(o,RooFit.RenameConflictNodes("(1)"))
#                if opts.verbosity>0 and not opts.quiet: o.Print()


###
###--- end of CAT loop
    print("models: %s"%model)
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
        print("AFTER CAT LOOP, CHI2 VALUE FOR %s IS %f"%(poln, chi2.getVal()))
    frametop = x["Pol3"].frame()
    framebot = x["Pol3"].frame()
    rh["Pol3"].plotOn(frametop)
    model["Pol3"].plotOn(frametop,RooFit.LineWidth(2), RooFit.LineColor(2)) #red
#    rh["Pol4"].plotOn(frametop)
    model["Pol4"].plotOn(frametop,RooFit.LineWidth(2), RooFit.LineColor(3)) #green
#    rh["Pol5"].plotOn(frametop)
    model["Pol5"].plotOn(frametop,RooFit.LineWidth(2), RooFit.LineColor(4)) #blue  
#    rh["Pol6"].plotOn(frametop)
    model["Pol6"].plotOn(frametop,RooFit.LineWidth(2), RooFit.LineColor(5)) #yellow
    print("model after plotting:")
    model["Pol6"].Print()
    print("fit result:")
    res["Pol6"].Print()
    
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
    frametop.GetXaxis().SetTitleSize(.055)
    frametop.GetXaxis().SetLabelSize(.035)
    frametop.GetYaxis().SetLabelSize(0.035)
    frametop.GetYaxis().SetTitleSize(0.055)
    frametop.GetYaxis().SetTitleOffset(1.5);
    frametop.GetYaxis().SetTitle("Events / %.1f GeV"%(opts.dX[0]*10))
    frametop.GetXaxis().SetTitle("M_{bb} (GeV)")
    frametop.Draw()
    gPad.Update()
    
    ndf = int((opts.X[1]-opts.X[0])/(opts.dX[0])) - (n_param)
    chi2=chi2_val/ndf
    prob = ROOT.TMath.Prob(chi2_val,ndf)
# part 2
    pad = TPad("pad","pad",0.,0.,1.,1.)
    pad.SetTopMargin(0.7)
    pad.SetFillStyle(-1)
    pad.SetBorderSize(0)
    pad.Draw()
    pad.cd(0)
#    framebot.addPlotable(rhres,"p")
    framebot.GetYaxis().SetTitle("Pulls")
    framebot.GetXaxis().SetTitle("M_{bb} (GeV)")
    framebot.GetYaxis().SetLabelSize(0.03)
    framebot.GetYaxis().SetTitleSize(0.04)
    line = ROOT.TLine(opts.X[0],0,opts.X[1],0)
    line.SetLineStyle(7)
    line.SetLineWidth(2)
    i = n_param
    line.SetLineColor(kRed+i)
    framebot.addObject(line)

    L00 = TLegend(0.7,0.6,1.-gStyle.GetPadRightMargin()-gStyle.GetPadTopMargin()*0.3333,1.-gStyle.GetPadTopMargin()*1.3333)
    L00.SetFillColor(-1)
    L00.SetBorderSize(0)
    L00.SetTextFont(42)
    L00.SetTextSize(gStyle.GetPadTopMargin()*0.60)
    L00.SetY1(L00.GetY2()-L00.GetNRows()*gStyle.GetPadTopMargin()*0.95)
   # if "6" in title:
    L00.Draw()
    framebot.addObject(L00)
    pave = TPaveText(0.75,0.7,0.9,0.9,"NDC")
    pave.SetTextAlign(11)
    pave.SetFillStyle(-1)
    pave.SetBorderSize(0)
    pave.SetTextFont(42)
    pave.SetTextSize(gStyle.GetPadTopMargin()*0.60)
    for i in range(3, 7):
        pave.AddText("Pol%d"%i) 
        pave.GetListOfLines().Last().SetTextColor(i-1);
    pave.Draw()
    framebot.addObject(pave)
    #now save canvas to file so can view.
    makeDirs("%s/plot/draw_funcs/"%opts.workdir)
    can2.SaveAs("%s/plot/draw_funcs/%s.pdf"%(opts.workdir,can2.GetName()))
    can2.SaveAs("%s/plot/draw_funcs/%s.png"%(opts.workdir,can2.GetName()))

    #PLOTTING RATIO OF POL6 / POL4 NOW
    canrat = TCanvas("cat2_Pol6_4_ratio","cat2_Pol6_4_ratio",600,600)
    canrat.cd()
    pCMS1.Draw()
    pCMS2.Draw()
    canrat.SetTitle("Pol3,4,5,6")
    canrat.SetFillColor(0)
    canrat.SetFillStyle(-1)
    nevents = 13142089 #from root file opened earlier
    print("*****************Pol6 fit values:")
#    print("res: %s"%res["Pol6"].printArgs())
    coefs = res["Pol6"].floatParsFinal()
    varNarrow = RooRealVar("varNar", "varNar", 100, 99, 101)
    setNarrow = RooArgSet()
    setNarrow.add(varNarrow)
    valNarrow = model["Pol6"].getVal(setNarrow)
    print("***GetVal (narrow range) result: %f" %valNarrow)
    xname = x["Pol6"].GetName()
    fitHist = model["Pol6"].createHistogram(xname, 20)
    fitHist.Print()
    fitHist2 = model["Pol3"].createHistogram(xname, 20)
    fitHist2.Print()
    for i in range(20):
        print("%f \t %f \t %f" %(fitHist.GetBinContent(i), fitHist2.GetBinContent(i), fitHist.GetBinError(i)) )
    print("Integral: " + str(fitHist.Integral()))
    
    #final fit parameters for Pol6, to be used to calculate ratio
    final_params_6 = []
    for i in range(8):
        param = coefs.at(i)
        print("parameter %d is %s"%(i, param))
        final_params_6.append(param.getValV())
    print("FINAL PARAMS: %s"%final_params_6)
    #x values (will be in range 80-220 GeV)
    xvals = []
    #corresponding values of the Pol6
    val6s = []
    #debugging: are these parameters correct?
    mycalcs = []
    #j is to count the number of entries.
    j = 0
    for mass in range(85, 196, 10):
        cheby = [1 for k in range(8)]
        cheby[1] = mass
        mycalcs.append(0)
        for k in range(2, 7):
            cheby[k] = 2*mass*cheby[k-1] - cheby[k-2]
        #value of Pol6 at this x value
        xvals.append(mass)
   #     val6s.append(1)
        #need to use the correct chebychev polynomial
        # defined as T0(x)=1, T1(x)=x, Tn+1(x)=2xTn(x)-Tn-1(x)
#        cheby = [1 for i in range(8)]
        #for some godforsaken reason the coeffs start with 1 instead of 0 (0 has to be 1???)
#        for i in range(1, 8):
#            if i == 1: cheby[i] = mass
#            else: cheby[i] = 2*mass*cheby[i-1] - cheby[i-2]
#            addval = final_params_6[i-1]*cheby[i]
#            print("mass=%f, i=%d, adding value %f"%(mass, i, addval))
#            val6s[j] += addval
        print(final_params_6)
        #val6s.append( ROOT.Math.Chebyshev6(mass, final_params_6[0], final_params_6[1], final_params_6[2],final_params_6[3],final_params_6[4],final_params_6[5],final_params_6[6]))
        #first parameter is 'implicitly' assumed to be 1 for some reason?
        val6s.append( ROOT.Math.Chebyshev6(mass, final_params_6[0], final_params_6[1], final_params_6[2],final_params_6[3],final_params_6[4],final_params_6[5], final_params_6[6]) / final_params_6[7])
        print("mass=%f, val6=%f"%(mass, val6s[j]))
        #last parameter is normalization for the whole func (I think) or maybe not
#        val6s[j] *= final_params_6[6]
        j += 1
    print("final_params_6=%s, j=%d, xvals=%s, val6s=%s"%(str(final_params_6), j, xvals, val6s))
    #now convert lists to 'arrays' so pyroot can understand it with its low IQ
    xs = array('f', xvals)
    y6 = array('f', val6s)
    print("val6s: %s \n y6: %s"%(str(val6s), str(y6)))
#    mc = array('f', mycalcs)
    graphRat = TGraph(j, xs, y6)
#    graphMc = TGraph(j, xs, mc)
    graphRat.GetXaxis().SetTitle("m_{bb} (GeV)")
    graphRat.GetYaxis().SetTitle("Ratio Pol6 / Pol4")
    graphRat.Draw("AC")
    #plot the RooFitResult res instead of the model
#    frameRat = x["Pol6"].frame()
#    res["Pol6"].plotOn(frameRat, ) #pink

#    graphMc.Draw("same")
    graphRat.Paint()
    canrat.Modified()
    canrat.Update()
    canrat.SaveAs("%s/plot/draw_funcs/%s.pdf"%(opts.workdir, canrat.GetName()))
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
