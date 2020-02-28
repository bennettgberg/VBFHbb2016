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
#   mp.add_option('','',help=colours[5]+''+colours[0],default='',type='',dest='')
#
    mp.add_option('--workdir',help=colours[5]+'Case workdir.'+colours[0],default='case0',type='str')
    mp.add_option('--long',help=colours[5]+'Long filenames.'+colours[0],default=False,action='store_true')
    mp.add_option('--bias',help=colours[5]+'Bias filenames.'+colours[0],default=False,action='store_true')
    mp.add_option('-v','--verbosity',help=colours[5]+'Verbosity.'+colours[0],default=0,action='count')
    mp.add_option('-q','--quiet',help=colours[5]+'Quiet.'+colours[0],default=True,action='store_false')

    mg1 = OptionGroup(mp,'Selection setup')
    mg1.add_option('--SELCATs',help=colours[5]+'Selection/Category setup.'+colours[0],default='double;DoubleB;-1#0.0875#0.5775#0.7875#1.,single;SingleB;-1#0.0775#0.5775#0.7775#0.8775#1.',type='str',action='callback',callback=SELsetup,dest='SC')# final
    mp.add_option_group(mg1)

    mg2 = OptionGroup(mp,'Datacard settings')
    mg2.add_option('--bounds',help=colours[5]+'Template boundaries: 80,200 (min,max)'+colours[0],default=[80.,200.],type='str',action='callback',callback=optsplitfloat,dest='X')
    mg2.add_option('--lumi',help=colours[5]+'Luminosity: 35900(single,double)'+colours[0],default=[35900,35900.],type='str',action='callback',callback=optsplitfloat,dest='LUMI')
    mg2.add_option('--TF',help=colours[5]+'Transfer function label: POL1,POL1 (Single,Double)'+colours[0],default=['POL1','POL1'],type='str',action='callback',callback=optsplit)
    mg2.add_option('--BRN',help=colours[5]+'Bernstein order: 5,5 '+colours[0],default=[5,5],type='str',action='callback',callback=optsplitint)
#   mg2.add_option('--MASS',help=colours[5]+'Signal masspoint: 125'+colours[0],default=[115,120,125,130,135],type='int')
    mg2.add_option('--MASS',help=colours[5]+'Signal masspoint: 125'+colours[0],default=[125],type='int')
    mg2.add_option('--CATveto',help=colours[5]+'CATveto: (blank)'+colours[0],default=[],type='str',action='callback',callback=optsplitint)
    mg2.add_option('--CATS',help=colours[5]+'CATMIN-CATMAX: 0,8'+colours[0],default=[0,8],type='str',action='callback',callback=optsplitint)
    mg2.add_option('--function',help=colours[5]+'Function for bias study :  expPow'+colours[0],default='expPow',type='str')
    mg2.add_option('--prefix',help=colours[5]+'Datacard prefix (folder): datacards/'+colours[0],default="datacards/",type='str')
    mp.add_option_group(mg2)
#   mg1.add_option('','',help=colours[5]+''+colours[0],default='',type='',dest='')
    mp.add_option_group(mg1)
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

####################################################################################################
def main():
    mp = parser()
    opts,args = mp.parse_args()
# Style
    style()

# Create directories if needed
    makeDirs('%s'%opts.workdir)
    makeDirs('%s/plot'%opts.workdir)
    makeDirs('%s/root'%opts.workdir)
    longtag1 = "_B%d-%d"%(opts.X[0],opts.X[1])
    longtag2 = "_B%d-%d_TF%s-%s"%(opts.X[0],opts.X[1],opts.TF[0],opts.TF[1])
    longtag3 = "_B%d-%d_BRN%d-%d_TF%s-%s"%(opts.X[0],opts.X[1],opts.BRN[0],opts.BRN[1],opts.TF[0],opts.TF[1])
    
# Settings
    #set seven_cats true if want to discard cats 0 and 4.
    seven_cats = True 

    SC = opts.SC if not type(opts.SC)==str else SELsetup(opts.SC)
    TF = opts.TF
#   currentMASS = [opts.MASS]
#   MASS = [115,120,125,130,135]
    MASS = [125]
#   MASS = [opts.MASS]
    ncats = sum(SC.ncats)
    #nmass = len(MASS)
    prefix = (opts.prefix).strip('/')
    makeDirs('%s/%s'%(opts.workdir,prefix))
    CATS = list(set(range(opts.CATS[0],opts.CATS[1]+1)) - set(opts.CATveto))
    #if seven_cats is set to true, then remove cats 0 and 4
    if seven_cats:
        del CATS[4]
        del CATS[0]
        print("WARNING: only using 7 categories. To use all 9 change seven_cats to false in the code.")
    print("Length of CATS: %d"% len(CATS))

# Read files
    uncertainties = json.loads(filecontent("%s/vbfHbb_uncertainties_2016.json"%basepath))
    fBKG = TFile.Open("%s/root/bkg_shapes_workspace%s.root"%(opts.workdir,"" if not opts.long else longtag1),"read")
    fSIG = TFile.Open("%s/root/sig_shapes_workspace%s.root"%(opts.workdir,"" if not opts.long else longtag1),"read")
    fDAT = TFile.Open("%s/root/data_shapes_workspace%s.root"%(opts.workdir,"" if not opts.long else longtag3),"read")
    wBKG = fBKG.Get("w")
    wSIG = fSIG.Get("w")
    wDAT = fDAT.Get("w")

# Load yields
    yD,yZ,yT,yVBF,yGF = [],[],[],[],[]
    for i in range(sum(SC.ncats)):
        yD += [wDAT.var("yield_data_CAT%d"%i).getValV()]
        yZ += [wDAT.var("yield_ZJets_CAT%d"%i).getValV()]
        yT += [wDAT.var("yield_Top_CAT%d"%i).getValV()]
    #   for j,jmass in enumerate(range(115,140,5)):
        for j,jmass in enumerate(range(125,130,10)):
            if j==0: 
                yVBF += [[]]
                yGF  += [[]]
            if not jmass in MASS: continue
            yVBF[i] += [wSIG.var("yield_signalVBF_mass%d_CAT%d"%(jmass,i)).getValV()]
            yGF[i]  += [wSIG.var("yield_signalGF_mass%d_CAT%d"%(jmass,i)).getValV()]
# Print card
#   for j,jmass in enumerate(range(115,140,5)):
    for j,jmass in enumerate(range(125,130,10)):
        print jmass
        if not jmass in MASS: continue
        if opts.bias :
    #       fout = open("%s/%s/datacard_vbfHbb_combo_bias%s_m%d%s%s.txt"%(opts.workdir,prefix,opts.function,jmass,"" if CATS==range(7) else "_CAT%d-CAT%d_CATveto%s"%(opts.CATS[0],opts.CATS[1],"".join(["%d"%x for x in opts.CATveto])),"" if not opts.long else longtag3),"w+")
            if not seven_cats:
                fout = open("%s/%s/datacard_vbfHbb_bias%s_m%d%s%s.txt"%(opts.workdir,prefix,opts.function,jmass,"" if CATS==range(7) else "_CAT%d-CAT%d_CATveto%s"%(opts.CATS[0],opts.CATS[1],"".join(["%d"%x for x in opts.CATveto])),"" if not opts.long else longtag3),"w+")
            else:
                fout = open("%s/%s/datacard_vbfHbb_bias%s_m%d%s%s_7cats.txt"%(opts.workdir,prefix,opts.function,jmass,"" if CATS==range(7) else "_CAT%d-CAT%d_CATveto%s"%(opts.CATS[0],opts.CATS[1],"".join(["%d"%x for x in opts.CATveto])),"" if not opts.long else longtag3),"w+")
            fout_table = open("%s/%s/table_bias%s.txt"%(opts.workdir,prefix,opts.function),"w+")
        if not opts.bias :
    #       fout = open("%s/%s/datacard_vbfHbb_combo_m%d%s%s.txt"%(opts.workdir,prefix,jmass,"" if CATS==range(7) else "_CAT%d-CAT%d_CATveto%s"%(opts.CATS[0],opts.CATS[1],"".join(["%d"%x for x in opts.CATveto])),"" if not opts.long else longtag3),"w+")
            fout = open("%s/%s/datacard_vbfHbb_names_m%d%s%s.txt"%(opts.workdir,prefix,jmass,"" if CATS==range(7) else "_CAT%d-CAT%d_CATveto%s"%(opts.CATS[0],opts.CATS[1],"".join(["%d"%x for x in opts.CATveto])),"" if not opts.long else longtag3),"w+")
            fout_table = open("%s/%s/table.txt"%(opts.workdir,prefix),"w+")
 ## numbers
        fout.write("imax %d\n"%len(CATS))
        fout.write("jmax *\nkmax *\n")
        fout.write("-"*100+"\n")
        fname = os.path.split(fDAT.GetName())[1]
        fname_bias = "bias_shapes_workspace_%s.root"%opts.function;
        sfname = os.path.split(fSIG.GetName())[1]
 ## files
        if opts.bias : 
            fout.write("%-15s %-15s *   %s/%-30s %-30s\n"%("shapes","data_obs","../root",fname,"w:data_hist_$CHANNEL"))
            fout.write("%-15s %-15s *   %s/%-30s %-30s\n"%("shapes","qcd     ","../root",fname_bias,"w:qcd_model_%s_$CHANNEL"%(''.join(opts.function))))
            fout.write("%-15s %-15s *   %s/%-30s %-30s\n"%("shapes","top     ","../root",fname,"w:Top_model_$CHANNEL"))
            fout.write("%-15s %-15s *   %s/%-30s %-30s\n"%("shapes","zjets   ","../root",fname,"w:Z_model_$CHANNEL"))
            fout.write("%-15s %-15s *   %s/%-30s %-30s\n"%("shapes","qqH_hbb     ","../root",sfname,"w:signal_model_m%s_$CHANNEL"%jmass))
            fout.write("%-15s %-15s *   %s/%-30s %-30s\n"%("shapes","ggH_hbb     ","../root",sfname,"w:signal_model_m%s_$CHANNEL"%jmass))

        if not opts.bias : 
            fout.write("%-15s %-15s *   %s/%-30s %-30s\n"%("shapes","qcd     ","../root",fname,"w:qcd_model_%s_$CHANNEL"%(''.join(opts.TF))))
            fout.write("%-15s %-15s *   %s/%-30s %-30s\n"%("shapes","data_obs","../root",fname,"w:data_hist_$CHANNEL"))
            fout.write("%-15s %-15s *   %s/%-30s %-30s\n"%("shapes","top     ","../root",fname,"w:Top_model_$CHANNEL"))
            fout.write("%-15s %-15s *   %s/%-30s %-30s\n"%("shapes","zjets   ","../root",fname,"w:Z_model_$CHANNEL"))
            fout.write("%-15s %-15s *   %s/%-30s %-30s\n"%("shapes","qqH_hbb     ","../root",sfname,"w:signal_model_m%s_$CHANNEL"%jmass))
            fout.write("%-15s %-15s *   %s/%-30s %-30s\n"%("shapes","ggH_hbb     ","../root",sfname,"w:signal_model_m%s_$CHANNEL"%jmass))
 ## top info
        fout.write("-"*100+"\n")
        fout.write("%-15s"%"bin")
        for i in CATS:
            fout.write("%-10s"%("CAT%d"%i))
        fout.write("\n")
        fout.write("%-15s"%"observation")
        for i in CATS: 
            fout.write("%-10d"%(-1))
        fout.write("\n"+"-"*100+"\n")
        fout.write("%-15s"%"bin")
        labels = ["qqH_hbb","ggH_hbb","qcd","top","zjets"]
        for i in CATS: 
            for l in labels:
                fout.write("%-8s"%("CAT%d"%i))
        fout.write("\n")
        fout.write("%-15s"%"process")
        for i in CATS: 
            for l in labels:
                fout.write("%-8s"%l)
        fout.write("\n")
        fout.write("%-15s"%"process")
        for i in CATS: 
            fout.write("%-8d"%(0))
            fout.write("%-8d"%(-1))
            fout.write("%-8d"%(1))
            fout.write("%-8d"%(1))
            fout.write("%-8d"%(1))
        fout.write("\n")
        fout.write("%-15s"%"rate")
        for i in CATS:
            fout.write("%-11s"%("%.2f"%yVBF[i][j]))
            fout.write("%-11s"%("%.2f"%yGF[i][j]))
            fout.write("%-11d"%(yD[i]))
            fout.write("%-11s"%("%.2f"%yT[i]))
            fout.write("%-11s"%("%.2f"%yZ[i]))
            fout_table.write("CAT%d"%i)
            fout_table.write("\n")
            fout_table.write("%s\t%.2f"%("VBF",yVBF[i][j]))
            fout_table.write("\n")
            fout_table.write("%s\t%.2f"%("GF",yGF[i][j]))
            fout_table.write("\n")
            fout_table.write("%s\t%.2f"%("QCD",yD[i]))
            fout_table.write("\n")
            fout_table.write("%s\t%.2f"%("Top",yT[i]))
            fout_table.write("\n")
            fout_table.write("%s\t%.2f"%("Z/W + jets",yZ[i]))
            fout_table.write("\n")
            fout_table.write("\n")
        fout_table.close()
        fout.write("\n"+"-"*100+"\n")
  ## fixed parts
        labels=[]
    #   if not opts.bias : labels = ["BR_hbb","QCDscale_qqh","QCDscale_ggh","pdf_Higgs_qqbar","pdf_Higgs_gg","lumi_13TeV","QGL_vbfbb","CMS_scale_j_ACCEPT_13TeV","CMS_res_j_ACCEPT_13TeV","CMS_qqH_hbb_trigger_13TeV","UE_vbfbb_13TeV","PS_vbfbb_13TeV","CMS_qqH_hbb_QCDscale_13TeV","CMS_ggH_hbb_QCDscale_13TeV","CMS_qqH_hbb_pdf_13TeV","CMS_ggH_hbb_pdf_13TeV","CMS_vbfbb_pu_unc_13TeV","CMS_vbfbb_bTagJES_13TeV","CMS_vbfbb_bTagLF_13TeV","CMS_vbfbb_bTagHF_13TeV","CMS_vbfbb_bTagLFStats1_13TeV","CMS_vbfbb_bTagLFStats2_13TeV","CMS_vbfbb_bTagHFStats1_13TeV","CMS_vbfbb_bTagHFStats2_13TeV","CMS_vbfbb_bTagCErr1_13TeV","CMS_vbfbb_bTagCErr2_13TeV"]
        if not opts.bias : labels = ["BR_hbb","QCDscale_qqH","QCDscale_ggH","pdf_Higgs_qqbar","pdf_Higgs_gg","lumi_13TeV_2016","QGL_vbfbb","CMS_scale_j_ACCEPT_13TeV","CMS_res_j_ACCEPT_13TeV","CMS_qqH_hbb_trigger_13TeV","UE","PS","QCDscale_qqH_ACCEPT","QCDscale_ggH_ACCEPT","pdf_qqH_ACCEPT","pdf_ggH_ACCEPT","CMS_vbfbb_pu_unc_13TeV","CMS_vbfbb_bTagJES_13TeV","CMS_vbfbb_bTagLF_13TeV","CMS_vbfbb_bTagHF_13TeV","CMS_vbfbb_bTagLFStats1_13TeV","CMS_vbfbb_bTagLFStats2_13TeV","CMS_vbfbb_bTagHFStats1_13TeV","CMS_vbfbb_bTagHFStats2_13TeV","CMS_vbfbb_bTagCErr1_13TeV","CMS_vbfbb_bTagCErr2_13TeV"]
    #   if opts.bias : labels = ["QCDscale_ggh","pdf_Higgs_qqbar","pdf_Higgs_gg","lumi_13TeV","CMS_scale_j_ACCEPT_13TeV","CMS_qqH_hbb_trigger_13TeV","UE_vbfbb_13TeV","PS_vbfbb_13TeV","CMS_qqH_hbb_QCDscale_13TeV","CMS_ggH_hbb_QCDscale_13TeV","CMS_vbfbb_bTagLF_13TeV"]
        if opts.bias : labels = ["QCDscale_ggH","pdf_Higgs_qqbar","pdf_Higgs_gg","lumi_13TeV_2016","CMS_scale_j_ACCEPT_13TeV","CMS_qqH_hbb_trigger_13TeV","UE","PS","QCDscale_qqH_ACCEPT","QCDscale_ggH_ACCEPT","CMS_vbfbb_bTagLF_13TeV"]
        for l in labels:
            fout.write("%-25s %-8s"%(l,"lnN"))
            for i in CATS: 
                k = j if len(uncertainties[l][0])==5 else i if len(uncertainties[l][0])==ncats else 0
                if len(uncertainties[l])==2: 
                    fout.write("%-10s %-10s %-10s %-10s %-10s"%("%g"%uncertainties[l][0][k],"%g"%uncertainties[l][1][k],"-","-","-"))
                elif len(uncertainties[l])==1 and "qq" in l:
                    fout.write("%-10s %-10s %-10s %-10s %-10s"%("%g"%uncertainties[l][0][k],"-","-","-","-"))
                elif len(uncertainties[l])==1 and "UE" in l:
                    fout.write("%-10s %-10s %-10s %-10s %-10s"%("%g"%uncertainties[l][0][k],"-","-","-","-"))
                elif len(uncertainties[l])==1 and "gg" in l:
                    fout.write("%-10s %-10s %-10s %-10s %-10s"%("-","%g"%uncertainties[l][0][k],"-","-","-"))
                elif len(uncertainties[l])==1:
                    fout.write("%-10s %-10s %-10s %-10s %-10s"%("%g"%uncertainties[l][0][k],"%g"%uncertainties[l][0][k],"-","-","-"))
            fout.write("\n")
        fout.write("-"*100+"\n")
  ## QCD part
        for i in CATS: 
            fout.write("%-25s %-8s"%("CMS_vbfbb_qcd_norm_CAT%d_13TeV"%i,"lnU"))
            for j in CATS: 
                if i==j: fout.write("%-10s %-10s %-10s %-10s %-10s"%("-","-","%g"%uncertainties['qcd_norm'],"-","-"))
                else:    fout.write("%-10s %-10s %-10s %-10s %-10s"%("-","-","-","-","-"))
            fout.write("\n")
        fout.write("-"*100+"\n")
  ## Z part
        for i in CATS: 
            fout.write("%-25s %-8s"%("CMS_vbfbb_zjets_norm_CAT%d_13TeV"%i,"lnN"))
            for j in CATS: 
                if i==j: fout.write("%-10s %-10s %-10s %-10s %-10s"%("-","-","-","-","%g"%uncertainties['zjets_norm']))
                else:    fout.write("%-10s %-10s %-10s %-10s %-10s"%("-","-","-","-","-"))
            fout.write("\n")
  ## T part
        for i in CATS: 
            fout.write("%-25s %-8s"%("CMS_vbfbb_top_norm_CAT%d_13TeV"%i,"lnN"))
            for j in CATS: 
                if i==j: fout.write("%-10s %-10s %-10s %-10s %-10s"%("-","-","-","%g"%uncertainties['top_norm'],"-"))
                else:    fout.write("%-10s %-10s %-10s %-10s %-10s"%("-","-","-","-","-"))
            fout.write("\n")
        fout.write("-"*100+"\n")
  ## scale/res part
        labels = ["CMS_vbfbb_scale_mbb_selsingle_13TeV_2016","CMS_vbfbb_scale_mbb_seldouble_13TeV_2016","CMS_vbfbb_res_mbb_selsingle_13TeV_2016","CMS_vbfbb_res_mbb_seldouble_13TeV_2016"]
        for l in labels:
            fout.write("%-30s %-10s %9s %9s\n"%(l,"param","%.1f"%uncertainties[l][0],"%g"%uncertainties[l][1]))
        fout.write("\n"+"-"*100+"\n")
  ## mean/sigma part
        labels = ["mean","sigma"]
        for l in labels:
            for i in CATS: 
                n = "%s_m%d_CAT%d_13TeV_2016"%(l,jmass,i)
                v = wSIG.var(n)
                fout.write("%-30s %-10s %9s %9s\n"%(n,"param","%g"%v.getValV(),"%g"%v.getError()))
        labels = ["Top_mean","Top_sigma","Z_mean","Z_sigma"]
        for l in labels:
            for i in CATS: 
                n = "%s_CAT%d_13TeV_2016"%(l,i)
                v = wDAT.var(n)
                fout.write("%-30s %-10s %9s %9s\n"%(n,"param","%g"%v.getValV(),"%g"%(v.getError() if not 'Top' in l else v.getValV()*(0.01 if 'mean' in l else 0.05))))
  ## transfer functions
        for iS,S in enumerate(SC.selections):
            start = sum(SC.ncats[:iS])
            end   = start + S.ncat
            if 'Fix' in opts.TF[iS]:
                for i in CATS: 
                    if not (i>start and i<end): continue
                    o = int(re.search('.*([0-9]{1})',opts.TF[iS]).group(1))
                    for j in range(o+1):
                        v = wDAT.var("trans_%s_CAT%d_p%d"%(opts.TF[iS],i,j))
                        fout.write("%-30s %-10s %9.2e %9.2e\n"%(v.GetName(),"param",v.getValV(),v.getError()))
        
# Close
        fout.close()
        fDAT.Close()
        fSIG.Close()
        fBKG.Close()

if __name__=='__main__':
    main()
