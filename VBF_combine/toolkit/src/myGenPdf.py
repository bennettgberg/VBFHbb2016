
import ROOT
from ROOT import RooFit,RooRealVar
import math
import sys
import os
from  pdf_param_cfi import *

#real_cat is the real category, 0-8 (for naming purposes only)
def generate_pdf(x=ROOT.RooRealVar(), pdf_name="Pol5",x_name='mbbReg_CAT0', selection="double",gcs=[], real_cat=0, val140=10000):
    pdf = None
    coeff = ROOT.RooArgList()
#    n_param = Nparam[pdf_name]
#    print n_param
    cat_num=0
#TEMPORARY COMMENTS, DON'T FORGET TO CHANGE BACK!!!!!!!!!!!!!!11
    if selection=="double" : cat_num=0
    if selection=="single" : cat_num=4
#    if not pdf_name in Parameters:
#        sys.exit("Error: pdf name %s is not recognized."%pdf_name)
    name = "qcd_model_%s_CAT%d"%(pdf_name, real_cat)
    #name = "qcd_model_%s_CAT%d"%(pdf_name, cat_num) #temporary change, for unified singleB and doubleB
    if pdf_name=="expPow":            
        coeff.removeAll()
        brn = {}
        brn_names = []
        for p in xrange(n_param):   
            nb = "b%d_CAT%d"%(p,real_cat)
        #   nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
            brn_names.append(nb)
        #   brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
            [pmin,pmax] = Parameters[pdf_name][nb] 
            brn[nb] = RooRealVar(nb,nb,pmin,pmax)
            coeff.add(brn[nb])
            gcs.append(brn[nb])
            const30 = '30'
        formula = "TMath::Power((%s-%s),%s)*TMath::Exp(-1*%s*TMath::Power(%s,%s))"%(x_name,const30,brn_names[0],brn_names[1],x_name,brn_names[2])
        coeff.add(x)
        pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="expPow2":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "TMath::Power((%s-%s),%s)*TMath::Exp(-1*%s*TMath::Power(%s,%s)-%s)"%(x_name,const30,brn_names[0],brn_names[1],x_name,brn_names[2],brn_names[3])
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif "x^Pol" in pdf_name:
                coeff.removeAll()
                brn = {}
                brn_names = []
                #x^Pols have one more param than their order (start indexing from 0)
                nparam = int(pdf_name[len(pdf_name)-1]) + 1
                formula = "TMath::Power(%s,(0+"%(x_name)
                for p in xrange(nparam):
                        nb = "b%d_%s_CAT%d"%(p,selection,real_cat)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                 #      brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                 #       [pmin,pmax] = Parameters[pdf_name][nb]
                        [pmin,pmax] = [-10,10]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        formula += "+%s*TMath::Power(%s, %d)"%(brn_names[p], x_name, p)
                #first parenthesis for closing the exponent, then for closing the power function
                formula += "))"
                print("formula:%s"%formula)
#                sys.exit(0)
                #formula = "TMath::Power(%s,(%s+%s*%s))"%(x_name,brn_names[0],brn_names[1],x_name) 
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="x^Pol1":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                 #      brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                 #       [pmin,pmax] = Parameters[pdf_name][nb]
                        [pmin,pmax] = [-10,10]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "TMath::Power(%s,(%s+%s*%s))"%(x_name,brn_names[0],brn_names[1],x_name)
        
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="x^Pol2":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                 #      brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "TMath::Power(%s,%s+%s*TMath::Power(%s,1)+%s*TMath::Power(%s,2))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="x^pol3":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                 #       nb = "b%d_CAT%d"%(p,real_cat)
                        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                 #      brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const = '100000000'
                formula = "TMath::Power(%s,%s+%s*TMath::Power(%s,1)+%s*TMath::Power(%s,2)+%s*TMath::Power(%s,3))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="x^Pol4":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                 #      brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "TMath::Power(%s,%s+%s*TMath::Power(%s,1)+%s*TMath::Power(%s,2)+%s*TMath::Power(%s,3)+%s*TMath::Power(%s,4))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="x^Pol5":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                 #      brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "TMath::Power(%s,%s+%s*TMath::Power(%s,1)+%s*TMath::Power(%s,2)+%s*TMath::Power(%s,3)+%s*TMath::Power(%s,4)+%s*TMath::Power(%s,5))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )    

    elif pdf_name=="Pol1exp":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "TMath::Exp(-30*%s)*((1+%s*TMath::Power(%s,1)))"%(x_name,brn_names[0],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="Pol2exp":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "TMath::Exp(-30*%s)*((1+%s*TMath::Power(%s,1)+%s*TMath::Power(%s,2)))"%(x_name,brn_names[0],x_name,brn_names[1],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="Pol1_dijet":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                                   
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                          
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "%s^(%s+%s*TMath::Log(%s))*(1+%s*TMath::Power(%s,1))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )


    elif pdf_name=="Pol2_dijet":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                                   
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                          
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "%s^(%s+%s*TMath::Log(%s))*(1+%s*TMath::Power(%s,1)+%s*TMath::Power(%s,2))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="pol3_dijet":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                #        nb = "b%d_CAT%d"%(p,real_cat)
                        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                            
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                   $
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "%s^(%s+%s*TMath::Log(%s))*(1+%s*TMath::Power(%s,1)+%s*TMath::Power(%s,2)+%s*TMath::Power(%s,3))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="Pol4_dijet":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                            $
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                   $
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "%s^(%s+%s*TMath::Log(%s))*(1+%s*TMath::Power(%s,1)+%s*TMath::Power(%s,2)+%s*TMath::Power(%s,3)+%s*TMath::Power(%s,4))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="Pol5_dijet":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                            $
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                   $
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "%s^(%s+%s*TMath::Log(%s))*(1+%s*TMath::Power(%s,1)+%s*TMath::Power(%s,2)+%s*TMath::Power(%s,3)+%s*TMath::Power(%s,4)+%s*TMath::Power(%s,5))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],x_name,brn_names[6],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="Pol6_dijet":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                            $
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                   $
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "%s^(%s+%s*TMath::Log(%s))*(1+%s*TMath::Power(%s,1)+%s*TMath::Power(%s,2)+%s*TMath::Power(%s,3)+%s*TMath::Power(%s,4)+%s*TMath::Power(%s,5)+%s*TMath::Power(%s,6))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],x_name,brn_names[6],x_name,brn_names[7],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="Pol7_dijet":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                            $
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                   $
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "%s^(%s+%s*TMath::Log(%s))*(1+%s*TMath::Power(%s,1)+%s*TMath::Power(%s,2)+%s*TMath::Power(%s,3)+%s*TMath::Power(%s,4)+%s*TMath::Power(%s,5)+%s*TMath::Power(%s,6)+%s*TMath::Power(%s,7))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],x_name,brn_names[6],x_name,brn_names[7],x_name,brn_names[8],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="Pol8_dijet":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                            $
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                   $
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "%s^(%s+%s*TMath::Log(%s))*(1+%s*TMath::Power(%s,1)+%s*TMath::Power(%s,2)+%s*TMath::Power(%s,3)+%s*TMath::Power(%s,4)+%s*TMath::Power(%s,5)+%s*TMath::Power(%s,6)+%s*TMath::Power(%s,7)+%s*TMath::Power(%s,8))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],x_name,brn_names[6],x_name,brn_names[7],x_name,brn_names[8],x_name,brn_names[9],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )


    elif pdf_name=="Pol1_invdijet":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                                   
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                          
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "(%s^(%s+%s*TMath::Log(%s)))*(1+%s*TMath::Power(1/%s,1))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="pol2_invdijet":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                #        nb = "b%d_CAT%d"%(p,real_cat)
                        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                                   
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                          
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "(%s^(%s+%s*TMath::Log(%s)))*(1+%s*TMath::Power(1/%s,1)+%s*TMath::Power(1/%s,2))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="Pol3_invdijet":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                                   
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                          
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "(%s^(%s+%s*TMath::Log(%s)))*(1+%s*TMath::Power(1/%s,1)+%s*TMath::Power(1/%s,2)+%s*TMath::Power(1/%s,3))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="Pol4_invdijet":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                                   
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                          
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "(%s^(%s+%s*TMath::Log(%s)))*(1+%s*TMath::Power(1/%s,1)+%s*TMath::Power(1/%s,2)+%s*TMath::Power(1/%s,3)+%s*TMath::Power(1/%s,4))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="Pol5_invdijet":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                                   
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                          
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "(%s^(%s+%s*TMath::Log(%s)))*(1+%s*TMath::Power(1/%s,1)+%s*TMath::Power(1/%s,2)+%s*TMath::Power(1/%s,3)+%s*TMath::Power(1/%s,4)+%s*TMath::Power(1/%s,5))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],x_name,brn_names[6],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="Pol6_invdijet":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                                   
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                          
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "(%s^(%s+%s*TMath::Log(%s)))*(1+%s*TMath::Power(1/%s,1)+%s*TMath::Power(1/%s,2)+%s*TMath::Power(1/%s,3)+%s*TMath::Power(1/%s,4)+%s*TMath::Power(1/%s,5)+%s*TMath::Power(1/%s,6))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],x_name,brn_names[6],x_name,brn_names[7],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="Pol7_invdijet":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,real_cat)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                                   
                        brn_names.append(nb)
                 #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                          
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "(%s^(%s+%s*TMath::Log(%s)))*(1+%s*TMath::Power(1/%s,1)+%s*TMath::Power(1/%s,2)+%s*TMath::Power(1/%s,3)+%s*TMath::Power(1/%s,4)+%s*TMath::Power(1/%s,5)+%s*TMath::Power(1/%s,6)+%s*TMath::Power(1/%s,7))"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],x_name,brn_names[6],x_name,brn_names[7],x_name,brn_names[8],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="modG":            
        coeff.removeAll()
        brn = {}
        brn_names = []
        for p in xrange(n_param):   
        #   nb = "b%d_%s_CAT0"%(p,pdf_name)
        #   nb = "b%d_CAT%d"%(p,real_cat)
            nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
            brn_names.append(nb)
        #   brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
            [pmin,pmax] = Parameters[pdf_name][nb] 
            brn[nb] = RooRealVar(nb,nb,pmin,pmax)
            coeff.add(brn[nb])
            gcs.append(brn[nb])
        formula = "TMath::Erfc(%s-%s*%s)*TMath::Exp(-1.*%s*%s)"%(brn_names[0],x_name,brn_names[1],brn_names[2],x_name)
        coeff.add(x)
        pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="modG2":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                #       nb = "b%d_%s_CAT0"%(p,pdf_name)
                        nb = "b%d_CAT%d"%(p,real_cat)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                formula = "TMath::Erfc(%s-%s*%s)*TMath::Exp(-1.*%s*%s-%s*TMath::Power(%s,2))"%(brn_names[0],x_name,brn_names[1],x_name,brn_names[2],brn_names[3],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="modG3":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                #       nb = "b%d_%s_CAT0"%(p,pdf_name)
                        nb = "b%d_CAT%d"%(p,real_cat)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                formula = "TMath::Erfc(%s-%s*%s)*TMath::Exp(-1.*%s*%s-%s*TMath::Power(%s,2)-%s*TMath::Power(%s,3))"%(brn_names[0],x_name,brn_names[1],x_name,brn_names[2],x_name,brn_names[3],brn_names[4],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="modG4":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                #       nb = "b%d_%s_CAT0"%(p,pdf_name)
                #       nb = "b%d_CAT%d"%(p,real_cat)
                        nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                formula = "TMath::Exp(-1.*%s*%s-%s*TMath::Power(%s,2)-%s*TMath::Power(%s,3)-%s*TMath::Power(%s,4))*TMath::Erfc(%s-%s*%s)"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],brn_names[2],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="modG5":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                #       nb = "b%d_%s_CAT0"%(p,pdf_name)
                        nb = "b%d_CAT%d"%(p,real_cat)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                formula = "TMath::Exp(-1.*%s*%s-%s*TMath::Power(%s,2)-%s*TMath::Power(%s,3)-%s*TMath::Power(%s,4)-%s*TMath::Power(%s,5))*TMath::Erfc(%s-%s*%s)"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],brn_names[6],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="modG6":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                #       nb = "b%d_%s_CAT0"%(p,pdf_name)
                        nb = "b%d_CAT%d"%(p,real_cat)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                formula = "TMath::Exp(-1.*%s*%s-%s*TMath::Power(%s,2)-%s*TMath::Power(%s,3)-%s*TMath::Power(%s,4)-%s*TMath::Power(%s,5)-%s*TMath::Power(%s,6))*TMath::Erfc(%s-%s*%s)"%(x_name,brn_names[0],brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],x_name,brn_names[6],brn_names[7],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="tanh2":            
        coeff.removeAll()
        brn = {}
        brn_names = []
        for p in xrange(n_param):   
        #   nb = "b%d_%s_CAT%d"%(p,pdf_name,real_cat)
        #   nb2 = "b%d_CAT%d"%(p,real_cat)
            nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
            brn_names.append(nb)
        #   brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
            [pmin,pmax] = Parameters[pdf_name][nb] 
            brn[nb] = RooRealVar(nb,nb,pmin,pmax)
            coeff.add(brn[nb])
            gcs.append(brn[nb])
        formula = "%s-%s*TMath::TanH(%s*%s-%s)-%s*TMath::Power(TMath::TanH(%s*%s-%s),2)"%(brn_names[0],brn_names[1],x_name,brn_names[2],brn_names[3],brn_names[4],x_name,brn_names[2],brn_names[3])
        coeff.add(x)
        pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="tanh":            
        coeff.removeAll()
        brn = {}
        brn_names = []
        for p in xrange(n_param):   
        #   nb = "b%d_%s_CAT0"%(p,pdf_name)
        #   nb2 = "b%d_CAT0"%(p)
            nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
        #   nb = "b%d_%s_CAT%d"%(p,pdf_name,real_cat)
            brn_names.append(nb)
        #   brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
            [pmin,pmax] = Parameters[pdf_name][nb] 
            brn[nb] = RooRealVar(nb,nb,pmin,pmax)
            coeff.add(brn[nb])
            gcs.append(brn[nb])
        formula = "%s-%s*TMath::TanH(%s*%s-%s)"%(brn_names[0],brn_names[1],x_name,brn_names[2],brn_names[3])
        coeff.add(x)
        pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="tanh3":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_%s_CAT%d"%(p,pdf_name,real_cat)
                #   nb2 = "b%d_CAT%d"%(p,real_cat)
                #   nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                #   brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                formula = "%s-%s*TMath::TanH(%s*%s-%s)-%s*TMath::Power(TMath::TanH(%s*%s-%s),2)-%s*TMath::Power(TMath::TanH(%s*%s-%s),3)"%(brn_names[0],brn_names[1],x_name,brn_names[2],brn_names[3],brn_names[4],x_name,brn_names[2],brn_names[3],brn_names[5],x_name,brn_names[2],brn_names[3])
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="tanh4":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_%s_CAT%d"%(p,pdf_name,real_cat)
                #       nb2 = "b%d_CAT%d"%(p,real_cat)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                formula = "%s-%s*TMath::TanH(%s*%s-%s)-%s*TMath::Power(TMath::TanH(%s*%s-%s),2)-%s*TMath::Power(TMath::TanH(%s*%s-%s),3)-%s*TMath::Power(TMath::TanH(%s*%s-%s),4)"%(brn_names[0],brn_names[1],x_name,brn_names[2],brn_names[3],brn_names[4],x_name,brn_names[2],brn_names[3],brn_names[5],x_name,brn_names[2],brn_names[3],brn_names[6],x_name,brn_names[2],brn_names[3])
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="tanh5":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_%s_CAT%d"%(p,pdf_name,real_cat)
                #       nb2 = "b%d_CAT%d"%(p,real_cat)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                formula = "%s-%s*TMath::TanH(%s*%s-%s)-%s*TMath::Power(TMath::TanH(%s*%s-%s),2)-%s*TMath::Power(TMath::TanH(%s*%s-%s),3)-%s*TMath::Power(TMath::TanH(%s*%s-%s),4)-%s*TMath::Power(TMath::TanH(%s*%s-%s),5)"%(brn_names[0],brn_names[1],x_name,brn_names[2],brn_names[3],brn_names[4],x_name,brn_names[2],brn_names[3],brn_names[5],x_name,brn_names[2],brn_names[3],brn_names[6],x_name,brn_names[2],brn_names[3],brn_names[7],x_name,brn_names[2],brn_names[3])
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="tanh6":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_%s_CAT%d"%(p,pdf_name,real_cat)
                #       nb2 = "b%d_CAT%d"%(p,real_cat)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                formula = "%s-%s*TMath::TanH(%s*%s-%s)-%s*TMath::Power(TMath::TanH(%s*%s-%s),2)-%s*TMath::Power(TMath::TanH(%s*%s-%s),3)-%s*TMath::Power(TMath::TanH(%s*%s-%s),4)-%s*TMath::Power(TMath::TanH(%s*%s-%s),5)-%s*TMath::Power(TMath::TanH(%s*%s-%s),6)"%(brn_names[0],brn_names[1],x_name,brn_names[2],brn_names[3],brn_names[4],x_name,brn_names[2],brn_names[3],brn_names[5],x_name,brn_names[2],brn_names[3],brn_names[6],x_name,brn_names[2],brn_names[3],brn_names[7],x_name,brn_names[2],brn_names[3],brn_names[8],x_name,brn_names[2],brn_names[3])
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="sine1":            
        coeff.removeAll()
        brn = {}
        brn_names = []
        for p in xrange(n_param):   
        #   nb = "b%d_%s_CAT0"%(p,pdf_name)
        #   nb = "b%d_CAT%d"%(p,real_cat)
            #nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
            nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
            brn_names.append(nb)
        #   brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
            [pmin,pmax] = Parameters[pdf_name][nb] 
            brn[nb] = RooRealVar(nb,nb,pmin,pmax)
            coeff.add(brn[nb])
            gcs.append(brn[nb])
        formula = "(1+%s*TMath::Sin(%s*%s+%s))"%(brn_names[0],x_name,brn_names[1],brn_names[2])
        coeff.add(x)
        pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="sine2":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                #       nb = "b%d_%s_CAT0"%(p,pdf_name)                                                       
                        #nb = "b%d_CAT%d"%(p,real_cat)
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)                                          
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)                                                 
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                formula = "(1+%s*TMath::Sin(%s*%s+%s)+%s*TMath::Power(TMath::Sin(%s*%s+%s),2))"%(brn_names[0],x_name,brn_names[1],brn_names[2],brn_names[3],x_name,brn_names[1],brn_names[2])
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )
    elif pdf_name=="sine3":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                #   nb = "b%d_%s_CAT0"%(p,pdf_name)
                        #nb = "b%d_CAT%d"%(p,real_cat)
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                formula = "(1+%s*TMath::Sin(%s*%s+%s)+%s*TMath::Power(TMath::Sin(%s*%s+%s),2)+%s*TMath::Power(TMath::Sin(%s*%s+%s),3))"%(brn_names[0],x_name,brn_names[1],brn_names[2],brn_names[3],x_name,brn_names[1],brn_names[2],brn_names[4],x_name,brn_names[1],brn_names[2])
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="sine4":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                #       nb = "b%d_%s_CAT0"%(p,pdf_name)
                        #nb = "b%d_CAT%d"%(p,real_cat)
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                formula = "(1+%s*TMath::Sin(%s*%s+%s)+%s*TMath::Power(TMath::Sin(%s*%s+%s),2)+%s*TMath::Power(TMath::Sin(%s*%s+%s),3)+%s*TMath::Power(TMath::Sin(%s*%s+%s),4))"%(brn_names[0],x_name,brn_names[1],brn_names[2],brn_names[3],x_name,brn_names[1],brn_names[2],brn_names[4],x_name,brn_names[1],brn_names[2],brn_names[5],x_name,brn_names[1],brn_names[2])
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="sine5":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                #        nb = "b%d_%s_CAT0"%(p,pdf_name)
                      #  nb = "b%d_CAT%d"%(p,real_cat)
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                formula = "(1+%s*TMath::Sin(%s*%s+%s)+%s*TMath::Power(TMath::Sin(%s*%s+%s),2)+%s*TMath::Power(TMath::Sin(%s*%s+%s),3)+%s*TMath::Power(TMath::Sin(%s*%s+%s),4)+%s*TMath::Power(TMath::Sin(%s*%s+%s),5))"%(brn_names[0],x_name,brn_names[1],brn_names[2],brn_names[3],x_name,brn_names[1],brn_names[2],brn_names[4],x_name,brn_names[1],brn_names[2],brn_names[5],x_name,brn_names[1],brn_names[2],brn_names[6],x_name,brn_names[1],brn_names[2])
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif pdf_name=="sineErf":            
        coeff.removeAll()
        brn = {}
        brn_names = []
        for p in xrange(n_param):   
        #   nb = "b%d_%s_CAT0"%(p,pdf_name)
            nb2 = "b%d_CAT%d"%(p,real_cat)
            nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
            brn_names.append(nb)
        #   brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
            [pmin,pmax] = Parameters[pdf_name][nb2] 
            brn[nb] = RooRealVar(nb,nb,pmin,pmax)
            coeff.add(brn[nb])
            gcs.append(brn[nb])
        formula = "(1+%s*TMath::Sin(%s*%s+%s))*TMath::Erfc(%s-%s*%s)"%(brn_names[0],x_name,brn_names[1],brn_names[2],brn_names[3],brn_names[4],x_name)
        coeff.add(x)
        pdf = ROOT.RooGenericPdf(name,"",formula,coeff )
    elif pdf_name=="erfPol2":            
        coeff.removeAll()
        brn = {}
        brn_names = []
        for p in xrange(n_param):   
        #   nb = "b%d_%s_CAT0"%(p,pdf_name)
            nb2 = "b%d_CAT0"%(p)
            nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
            brn_names.append(nb)
        #   brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
            [pmin,pmax] = Parameters[pdf_name][nb2] 
            brn[nb] = RooRealVar(nb,nb,pmin,pmax)
            coeff.add(brn[nb])
            gcs.append(brn[nb])
        formula = "TMath::Erf((%s-%s)/%s)*(%s*TMath::Power(%s,2)+%s*%s+1.)"%(x_name,brn_names[0],brn_names[1],brn_names[2],x_name,brn_names[3],x_name)
        coeff.add(x)
        pdf = ROOT.RooGenericPdf(name,"",formula,coeff )
    elif pdf_name=="erfPol3":            
        coeff.removeAll()
        brn = {}
        brn_names = []
        for p in xrange(n_param):   
        #   nb = "b%d_%s_CAT0"%(p,pdf_name)
            nb2 = "b%d_CAT0"%(p)
            nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
            brn_names.append(nb)
        #   brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
            [pmin,pmax] = Parameters[pdf_name][nb2] 
            brn[nb] = RooRealVar(nb,nb,pmin,pmax)
            coeff.add(brn[nb])
            gcs.append(brn[nb])
        formula = "TMath::Erf((%s-%s)/%s)*(%s*TMath::Power(%s,3)+%s*TMath::Power(%s,2)+%s*%s+1.)"%(x_name,brn_names[0],brn_names[1],brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name)
        coeff.add(x)
        pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

    elif "expPol" in pdf_name:
                coeff.removeAll()
                brn = {}
                brn_names = []
            #number of parameters (excluding the normalization) is the order of the expPol.
                nparam = int(pdf_name[len(pdf_name)-1])
                formula = "TMath::Exp(0"
            #first open the file for reading (we'll write to it later)
                filename = "expPol_params%d.txt"%real_cat
                read_from_file = (os.path.exists(filename) and nparam > 2)
                if read_from_file:
                    initfile = open(filename, "r")
            #starts at expPol2
                #we can read the init values of all but one parameter (expPol2 shouldn't read any).
                for p in range(nparam-2):
                    line = initfile.readline()
                sparams = []
                #only can read into sparams if the file to read from exists.
                if read_from_file:
                    sparams = line.split()
                if read_from_file and not len(sparams) == nparam-1:
                    print("Error! %d params read from file %s (expected %d params for %s"%(len(sparams), filename, nparam-1, pdf_name))
                    sys.exit(0)
                for p in xrange(nparam):
                    nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
                    brn_names.append(nb)
               #     init_val = (-0.01)**(nparam+1) #works ok!
                    np = nparam #min(nparam, 4)
                    init_val = 0.0
                    if p < nparam-1 and read_from_file:
                        init_val = float(sparams[p])
                    #elif p < nparam-1:
                    else:
                        init_val = (-0.01)**(nparam+1)
                    brn[nb] = RooRealVar(nb,nb,init_val,-1,1)
                    coeff.add(brn[nb])
                    gcs.append(brn[nb])
                    formula += " + " + brn_names[p] + "*TMath::Power(" + x_name + ", " + str(p+1) + ")"
                #formula = "TMath::Exp(-1*%s*TMath::Power(%s,1) + %s)"%(brn_names[0],x_name, brn_names[1])
                formula += ")"
                print(formula)
               # raw_input()
               #now add x (mass of bb)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )
                print("cat %d expPol%d fit: %s"%(real_cat, nparam, str(pdf)))
    elif "BernPol" in pdf_name:
        #Bernstein polynomial
        #Pol of order n has n+1 parameters.
        nparam = int(pdf_name[len(pdf_name)-1]) + 1
        coeff.removeAll()
        brn = {}
        brn_names = []
        for p in xrange(nparam):   
        #   nb = "b%d_%s_CAT0"%(p,pdf_name)
            nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
            brn_names.append(nb)
            brn[nb] = RooRealVar(nb,nb,-10,10)
        #need to square this before fitting. 
        #how to do that???
            coeff.add(brn[nb]**2)
            gcs.append(brn[nb]**2)
        pdf = ROOT.RooBernstein(name,"",x,coeff)

    elif "Pol" in pdf_name: 
# #   if pdf_name=="Pol2":            
        #Pol of order n has n+1 parameters.
        nparam = int(pdf_name[len(pdf_name)-1]) + 1
        coeff.removeAll()
        brn = {}
        brn_names = []
        #TEMPORARY CHANGE, MAKE SURE YOU FIX LATER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #    sel = "double"#for fitcat0 
        for p in xrange(nparam):   
        #   nb = "b%d_%s_CAT0"%(p,pdf_name)
            nb = "b%d_sel%s_CAT%d"%(p,selection,real_cat)
            #nb = "b%d_sel%s_CAT%d"%(p,sel,cat_num) #for fitcat0
            brn_names.append(nb)
            brn[nb] = RooRealVar(nb,nb,-1,1)
            coeff.add(brn[nb])
            gcs.append(brn[nb])
        pdf = ROOT.RooChebychev(name,"",x,coeff )
    #form: a0*exp(x^0)+a1*exp(x^1)+...+an*exp(x^n)
    elif "expx" in pdf_name:
        nparam = int(pdf_name[len(pdf_name)-1]) + 1
        coeff.removeAll()
        brn = {}
        brn_names = []
        formula = ""
        for p in xrange(nparam):
            nb = "b%d_sel%s_CAT%d"%(p, selection, real_cat)
            if p != 0:
                formula += " + "
            formula += "%s*TMath::Exp(TMath::Power(%s,%d))"%(nb, x_name, p)
            brn_names.append(nb)
            brn[nb] = RooRealVar(nb, nb, 0, -10, 10)
            coeff.add(brn[nb])
            gcs.append(brn[nb])
        print("formula: %s"%formula)
        coeff.add(x)
        pdf = ROOT.RooGenericPdf(name, "", formula, coeff)
    else:
        sys.exit("Error: Unknown function %s"%pdf_name)

    return [pdf, coeff]
