#!/usr/bin/env python

import os,sys,re,json,datetime
from glob import glob
from array import array
from optparse import OptionParser,OptionGroup
import warnings
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

####################################################################################################
def parser(mp=None):
	if not mp: mp = OptionParser()
#	mp.add_option('','',help=colours[5]+''+colours[0],default='',type='',dest='')
#
#	mg1 = OptionGroup(mp,'GroupTitle')
#	mg1.add_option('','',help=colours[5]+''+colours[0],default='',type='',dest='')
#	mp.add_option_group(mg1)
#
	mp.add_option('--workdir',help=colours[5]+'Case workdir.'+colours[0],default='case0',type='str')
	mp.add_option('--long',help=colours[5]+'Long filenames.'+colours[0],default=False,action='store_true')
	mp.add_option('-v','--verbosity',help=colours[5]+'Verbosity.'+colours[0],default=0,action='count')
	mp.add_option('-q','--quiet',help=colours[5]+'Quiet.'+colours[0],default=True,action='store_false')
#	mp.add_option('-q','--quiet',help=colours[5]+'Quiet.'+colours[0],default=False,action='store_false')
#
	mg1 = OptionGroup(mp,'Selection setup')
#	mg1.add_option('--SELCATs',help=colours[5]+'Selection/Category setup: "Double;DoubleB;-1.,-0.0175,0.2925,1.".'+colours[0],default='Double;DoubleB;-1.#-0.0175#0.2925#1.',type='str',action='callback',callback=SELsetup,dest='SC')# final
	mg1.add_option('--SELCATs',help=colours[5]+'Selection/Category setup: "Double;DoubleB;-1.,-0.2275, 0.1125,0.3425,1.".'+colours[0],default='Double;DoubleB;-1.#-0.2275#0.1125#0.3425#1.',type='str',action='callback',callback=SELsetup,dest='SC')# final
	mp.add_option_group(mg1)
#
	mg2 = OptionGroup(mp,'Data template setup')
	mg2.add_option('--bounds',help=colours[5]+'Template boundaries: 60,170 (min,max)'+colours[0],default=[60.,170.],type='str',action='callback',callback=optsplitfloat,dest='X')
	mg2.add_option('--binwidth',help=colours[5]+'Template bin width: 1 (Double)'+colours[0],default=[1.],type='str',action='callback',callback=optsplitfloat,dest='dX')
	mg2.add_option('--lumi',help=colours[5]+'Luminosity: 35900. (double)'+colours[0],default=[35900.],type='str',action='callback',callback=optsplitfloat,dest='LUMI')
	mg2.add_option('--TF',help=colours[5]+'Transfer function label: POL1 (Double)'+colours[0],default=['POL1'],type='str',action='callback',callback=optsplit)
	mg2.add_option('--BRN',help=colours[5]+'Bernstein order: 8 (Double)'+colours[0],default=[2],type='str',action='callback',callback=optsplitint)
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
	gROOT.ProcessLineSync("gErrorIgnoreLevel = kWarning;")
	RooMsgService.instance().setSilentMode(kTRUE)
	for i in range(2): RooMsgService.instance().setStreamStatus(i,kFALSE)

def hStyle(h,i,sumw2=True):
	if sumw2: h.Sumw2()
	h.SetMarkerStyle([20,20,23,21][i])
	h.SetMarkerSize(1.2)
	h.SetMarkerColor([kBlack,kBlue,kRed,kGreen+2][i])
	h.SetLineColor([kBlack,kBlue,kRed,kGreen+2][i])

def gpStyle(g,i):
	g.SetFillColor([kBlack,kBlue,kRed,kGreen+2][i])
	g.SetFillStyle(3004)

def gaStyle(g,i):
	g.SetFillColor(kGray)
	g.SetFillStyle(1001)

####################################################################################################
def RooDraw(opts,can,C,x,rh,model,qPDF,wPDF,tPDF,archive):
	can.cd()
#	can.SetLogy()
	print 'plotting'
	gPad.SetBottomMargin(0.3)
	frametop = x.frame()
	framebot = x.frame()
# part 1
	rh.plotOn(frametop)
	model.plotOn(frametop,RooFit.LineWidth(2))


#	rhres = frametop.residHist()
	rhres = frametop.pullHist()
	p1 = RooArgSet(qPDF)
	p2 = RooArgSet(wPDF)
	p3 = RooArgSet(tPDF)
	model.plotOn(frametop,RooFit.Components(p1),RooFit.LineWidth(2),RooFit.LineColor(kBlack),RooFit.LineStyle(kDashed))
#	model.plotOn(frametop,RooFit.Components(p2),RooFit.LineWidth(2),RooFit.LineColor(kBlue))
#	model.plotOn(frametop,RooFit.Components(p3),RooFit.LineWidth(2),RooFit.LineColor(kGreen+1))
	frametop.GetXaxis().SetTitleSize(0)
	frametop.GetXaxis().SetLabelSize(0)
	frametop.GetYaxis().SetLabelSize(0.035)
	frametop.GetYaxis().SetTitleSize(0.055)
	frametop.GetYaxis().SetTitleOffset(1.5);
	frametop.GetYaxis().SetTitle("Events / %.1f GeV"%opts.dX[0])
	frametop.Draw()
	gPad.Update()
	
	ndf = 110 - (opts.BRN[0]+1)
	chi2 = frametop.chiSquare()  
	prob = ROOT.TMath.Prob(chi2*ndf,ndf)

# part 2
	pad = TPad("pad","pad",0.,0.,1.,1.)
	pad.SetTopMargin(0.7)
	pad.SetFillStyle(-1)
	pad.SetBorderSize(0)
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
	line.SetLineColor(kGreen+2)
	framebot.addObject(line)

	L00 = TLegend(0.7,0.6,1.-gStyle.GetPadRightMargin()-gStyle.GetPadTopMargin()*0.3333,1.-gStyle.GetPadTopMargin()*1.3333)
	L00.SetHeader("Z-boson search")
	L00.SetFillColor(-1)
	L00.SetBorderSize(0)
	L00.SetTextFont(42)
	L00.SetTextSize(gStyle.GetPadTopMargin()*0.60)
	L00.SetY1(L00.GetY2()-L00.GetNRows()*gStyle.GetPadTopMargin()*0.95)
	L00.Draw()
	framebot.addObject(L00)

	pave = TPaveText(0.75,0.75,0.9,0.85,"NDC")
	pave.SetTextAlign(11)
	pave.SetFillStyle(-1)
	pave.SetBorderSize(0)
	pave.SetTextFont(42)
	pave.SetTextSize(gStyle.GetPadTopMargin()*0.60)
	pave.SetTextColor(kBlue+1)
	pave.AddText("#chi^{2}=%.2f"%chi2)
	pave.AddText("Prob=%.2f"%prob)
	pave.Draw()
	framebot.addObject(pave)


	framebot.Draw()
	gPad.Update()
	archive += [pad]
	archive += [rhres]

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
	longtag1 = "_B%d"%(opts.X[0])
	longtag2 = "_B%d_TF%s"%(opts.X[0],opts.TF[0])
	longtag3 = "_B%d_BRN%d_TF%s"%(opts.X[0],opts.BRN[0],opts.TF[0])
	
# Setup
	SC = opts.SC if not type(opts.SC)==str else SELsetup(opts.SC)
	fTF = json.loads(filecontent("%s/vbfZbb_transfer_2015_run2_linear.json"%basepath))
	LUMI = opts.LUMI
	NBINS = [int((opts.X[1]-opts.X[0])/x) for x in opts.dX]	
	print NBINS

# Files
	fBKG = TFile.Open("%s/root/shapes_workspace%s.root"%(opts.workdir,"" if not opts.long else longtag1),"read")
	fSIG = TFile.Open("%s/root/shapes_workspace%s.root"%(opts.workdir,"" if not opts.long else longtag1),"read")
	fTRF = TFile.Open("%s/root/TransferFunctions%s.root"%(opts.workdir,"" if not opts.long else longtag2),"read")
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
	for iS,S in enumerate(SC.selections):
## Load tree
		fin = TFile.Open("/afs/cern.ch/work/n/nchernya/VBFHbb_2016/Zbb_combine/inputs/root/FitZbb_BTagCSV_analysis_double_trignone_v25_trigWeightNewRegBtagCut_ZbbBdtInclFit.root","read")
		T = fin.Get("Zbb/events")
## Containers
		brn          = {}
		brn_params   = {}
		qcd_pdf      = {}
		qcd_pdf_aux  = {}
		qcd_pdf_aux2 = {}
		model        = {}
		nQCD         = {}
		transfer     = {}
		wPDF,tPDF,qPDF,sPDF = {},{},{},{}
		yW, yZ,yT,yQ     =  {} ,{},{},{}
		can = TCanvas("canD","canD",600,600)
#		can.Divide(2,2)
## Category loop
#		for C in range(S.ncat):
		if (1>0):
			C=0
			Cp = C + sum([x for x in SC.ncats[0:iS]])
####################################################################################################
#### Start of RooFit part 
  ## BRN parameters
			for ib in range(opts.BRN[iS]+1):
				if C>0: continue
				nb = "b%d_CAT%d"%(ib,Cp)
				brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
				nbp = "CAT%d"%(Cp)
				if not nbp in brn_params: brn_params[nbp] = RooArgList()
				brn_params[nbp].add(brn[nb])
  ## x
  			x = RooRealVar("mbbReg_CAT%d"%Cp,"mbbReg_CAT%d"%Cp,opts.X[0],opts.X[1])

  ### QCD part
	  ### Containers
			N = "CAT%d"%(Cp)
			h[N]   = TH1F("hMbb_%s"%N,"hMbb_%s"%N,NBINS[iS],opts.X[0],opts.X[1])
			hb[N]  = TH1F("hMbb_blind_%s"%N,"hMbb_blind_%s"%N,NBINS[iS],opts.X[0],opts.X[1])
	  ### Define cut
			cut    = TCut("bdt_Zbb>%1.4f && bdt_Zbb<=%1.4f"%(S.boundaries[C],S.boundaries[C+1]))
			cutB   = TCut("bdt_Zbb>%1.4f && bdt_Zbb<=%1.4f && (mbbRegFSR<80 || mbbRegFSR>100)"%(S.boundaries[C],S.boundaries[C+1]))
	  ### Draw
	  		c0.cd()
			T.Draw("mbbRegFSR>>%s"%(h[N].GetName()),cut)
			T.Draw("mbbRegFSR>>%s"%(hb[N].GetName()),cutB)
		#	for i in range(NBINS[iS]):
		#		if h[N].GetBinContent(i+1) ==0 : print i
			#	print h[N].GetBinContent(i+1)
			#	print i
	  ### Yields
			Y[N]   = RooRealVar("yield_data_CAT%d"%Cp,"yield_data_CAT%d"%Cp,h[N].Integral())
	  ### Histograms
			rh[N]  = RooDataHist("data_hist_CAT%d"%Cp,"data_hist_CAT%d"%Cp,RooArgList(x),h[N])
			rhb[N] = RooDataHist("data_hist_blind_CAT%d"%Cp,"data_hist_blind_CAT%d"%Cp,RooArgList(x),hb[N])

  ### Model
	  ### Bernstein 0 or 4
			if C==0:
				for ib in range(opts.BRN[iS]+1):
					nb = "b%d_CAT%d"%(ib,Cp)
					brn[nb].setConstant(kFALSE)
			#	qcd_pdf_aux[N] = RooBernstein("qcd_model_CAT%d"%(Cp),"qcd_model_CAT%d"%(Cp),x,brn_params[N])
			#	qcd_pdf[N]     = qcd_pdf_aux[N]#RooAbsPdf(qcd_pdf_aux[N])
				const30 = '30'
				formula = "TMath::Power((mbbReg_CAT0-%s),b0_CAT0)*TMath::Exp(-1*b1_CAT0*TMath::Power(mbbReg_CAT0,b2_CAT0))"%(const30)
				brn_params[N].add(x);
				qcd_pdf[N] = ROOT.RooGenericPdf( "qcd_model_CAT0","",formula,brn_params[N] )

	  ### PDFs
	  		wPDF[N] = wBKG.pdf("W_model_CAT%d"%Cp)
			tPDF[N] = wBKG.pdf("Top_model_CAT%d"%Cp)
			sPDF[N] = wSIG.pdf("Z_model_CAT%d"%Cp)

	  ### Yields
	  		yW[N]    = wBKG.var("yield_WJets_CAT%d"%Cp)
	  		yT[N]    = wBKG.var("yield_Top_CAT%d"%Cp)
			yQ[N]    = RooRealVar("yield_QCD_CAT%d"%Cp,"yield_QCD_CAT%d"%Cp,1000,0,1e+10)
			yW[N].setConstant(kTRUE)
			yT[N].setConstant(kTRUE)
#			ySVBF = wSIG.var("yield_signalVBF_mass%d_CAT%d"%(MASS,Cp))
#			ySGF  = wSIG.var("yield_signalGF_mass%d_CAT%d"%(MASS,Cp))
#			yS    = RooRealVar("yield_signal_mass%d_CAT%d"%(MASS,Cp),"yield_signal_mass%d_CAT%d"%(MASS,Cp),ySVBF.getValV()+ySGF.getValV(),0.,400.)

  ### Combined model
	  		model[N] = RooAddPdf("bkg_model_CAT%d"%(Cp),"bkg_model_CAT%d"%(Cp),RooArgList(wPDF[N],tPDF[N],qcd_pdf[N]),RooArgList(yW[N],yT[N],yQ[N]))
			
  ### Fit
	  		res   = model[N].fitTo(rh[N],RooFit.Save())
			if opts.verbosity>0 and not opts.quiet: res.Print()

  ### Draw
  			RooDraw(opts,can,C,x,rh[N],model[N],qcd_pdf[N],wPDF[N],tPDF[N],archive)
			can.cd()
			pCMS1.Draw()
			pCMS2.Draw()

  ### Reset constant settings
			for ib in range(opts.BRN[iS]+1):
				if C==0:
					nb = "b%d_CAT%d"%(ib,Cp)
					brn[nb].setConstant(kFALSE if opts.forfit else kTRUE)  ## false for final fit, true if fixed
  ### Saving
			for o in [rh[N],rhb[N],model[N],Y[N]]:
				getattr(w,'import')(o,RooFit.RenameConflictNodes("(1)"))
				if opts.verbosity>0 and not opts.quiet: o.Print()
###
###--- end of CAT loop
###
		if opts.verbosity>0 and not opts.quiet: w.Print()

		makeDirs("%s/plot/biasFunctions/"%opts.workdir)
		can.SaveAs("%s/plot/biasFunctions/%s.pdf"%(opts.workdir,can.GetName()))
		can.SaveAs("%s/plot/biasFunctions/%s.png"%(opts.workdir,can.GetName()))
	
#
#--- end of SEL loop
#

	makeDirs("%s/root/"%opts.workdir)
	w.writeToFile("%s/root/bias_shapes_workspace%s.root"%(opts.workdir,"" if not opts.long else longtag3))
	print "Done."

####################################################################################################
if __name__=='__main__':
	main()
