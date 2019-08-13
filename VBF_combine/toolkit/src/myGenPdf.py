
import ROOT
from ROOT import RooFit,RooRealVar
import math
import sys
from  pdf_param_cfi import *

#real_cat is the real category, 0-8 (for naming purposes only)
def generate_pdf(x=ROOT.RooRealVar(), pdf_name="Pol5",x_name='mbbReg_CAT0', selection="double",gcs=[], real_cat=0):
	pdf = None
	coeff = ROOT.RooArgList()
	n_param = Nparam[pdf_name]
	print n_param
	cat_num=0
	if selection=="double" : cat_num=0
	if selection=="single" : cat_num=4
        if not pdf_name in Parameters:
                sys.exit("Error: pdf name %s is not recognized."%pdf_name)
        name = "qcd_model_%s_CAT%d"%(pdf_name, real_cat)
	if pdf_name=="expPow":            
		coeff.removeAll()
		brn = {}
		brn_names = []
		for p in xrange(n_param):   
			nb = "b%d_CAT%d"%(p,cat_num)
		#	nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
			brn_names.append(nb)
		#	brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
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
                        nb = "b%d_CAT%d"%(p,cat_num)
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


	elif pdf_name=="expPol1":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                #        nb = "b%d_CAT%d"%(p,cat_num)
                	nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
                        brn_names.append(nb)
                #	brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
                formula = "TMath::Exp(-1*%s*TMath::Power(%s,1))"%(brn_names[0],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

#	elif pdf_name=="exp_pol2": #BG comment bc wtf
        elif pdf_name=="expPol2":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                #        nb = "b%d_CAT%d"%(p,cat_num)
                        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
            #    formula = "TMath::Exp(1-1*%s*TMath::Power(%s,1)-1*%s*TMath::Power(%s,2)-1*%s*TMath::Power(%s,3)-1*%s*TMath::Power(%s,4))"%(brn_names[0],x_name,brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name)
                formula = "TMath::Exp(-1*%s*TMath::Power(%s,1)-1*%s*TMath::Power(%s,2))"%(brn_names[0],x_name,brn_names[1],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

	elif pdf_name=="expPol3":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
            #    formula = "TMath::Exp(1-1*%s*TMath::Power(%s,1)-1*%s*TMath::Power(%s,2)-1*%s*TMath::Power(%s,3)-1*%s*TMath::Power(%s,4))"%(brn_names[0],x_name,brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name)
                formula = "TMath::Exp(-1*%s*TMath::Power(%s,1)-1*%s*TMath::Power(%s,2)-1*%s*TMath::Power(%s,3))"%(brn_names[0],x_name,brn_names[1],x_name,brn_names[2],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

	elif pdf_name=="expPol4":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                 #       nb = "b%d_CAT%d"%(p,cat_num)
                        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
            #    formula = "TMath::Exp(1-1*%s*TMath::Power(%s,1)-1*%s*TMath::Power(%s,2)-1*%s*TMath::Power(%s,3)-1*%s*TMath::Power(%s,4))"%(brn_names[0],x_name,brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name)
                formula = "TMath::Exp(-1*%s*TMath::Power(%s,1)-1*%s*TMath::Power(%s,2)-1*%s*TMath::Power(%s,3)-1*%s*TMath::Power(%s,4))"%(brn_names[0],x_name,brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )
	
	elif pdf_name=="expPol5":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
            #    formula = "TMath::Exp(1-1*%s*TMath::Power(%s,1)-1*%s*TMath::Power(%s,2)-1*%s*TMath::Power(%s,3)-1*%s*TMath::Power(%s,4))"%(brn_names[0],x_name,brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name)
                formula = "TMath::Exp(-1*%s*TMath::Power(%s,1)-1*%s*TMath::Power(%s,2)-1*%s*TMath::Power(%s,3)-1*%s*TMath::Power(%s,4)-1*%s*TMath::Power(%s,5))"%(brn_names[0],x_name,brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

	elif pdf_name=="expPol6":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
            #    formula = "TMath::Exp(1-1*%s*TMath::Power(%s,1)-1*%s*TMath::Power(%s,2)-1*%s*TMath::Power(%s,3)-1*%s*TMath::Power(%s,4))"%(brn_names[0],x_name,brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name)
                formula = "TMath::Exp(-1*%s*TMath::Power(%s,1)-1*%s*TMath::Power(%s,2)-1*%s*TMath::Power(%s,3)-1*%s*TMath::Power(%s,4)-1*%s*TMath::Power(%s,5)-1*%s*TMath::Power(%s,6))"%(brn_names[0],x_name,brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

	elif pdf_name=="expPol7":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
            #    formula = "TMath::Exp(1-1*%s*TMath::Power(%s,1)-1*%s*TMath::Power(%s,2)-1*%s*TMath::Power(%s,3)-1*%s*TMath::Power(%s,4))"%(brn_names[0],x_name,brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name)
                formula = "TMath::Exp(-1*%s*TMath::Power(%s,1)-1*%s*TMath::Power(%s,2)-1*%s*TMath::Power(%s,3)-1*%s*TMath::Power(%s,4)-1*%s*TMath::Power(%s,5)-1*%s*TMath::Power(%s,6)-1*%s*TMath::Power(%s,7))"%(brn_names[0],x_name,brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],x_name,brn_names[6],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

	elif pdf_name=="expPol8":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                        const30 = '30'
            #    formula = "TMath::Exp(1-1*%s*TMath::Power(%s,1)-1*%s*TMath::Power(%s,2)-1*%s*TMath::Power(%s,3)-1*%s*TMath::Power(%s,4))"%(brn_names[0],x_name,brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name)
                formula = "TMath::Exp(1-1*%s*TMath::Power(%s,1)-1*%s*TMath::Power(%s,2)-1*%s*TMath::Power(%s,3)-1*%s*TMath::Power(%s,4)-1*%s*TMath::Power(%s,5)-1*%s*TMath::Power(%s,6)-1*%s*TMath::Power(%s,7)-1*%s*TMath::Power(%s,8))"%(brn_names[0],x_name,brn_names[1],x_name,brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name,brn_names[5],x_name,brn_names[6],x_name,brn_names[7],x_name)
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

	elif pdf_name=="x^Pol1":
                coeff.removeAll()
                brn = {}
                brn_names = []
                for p in xrange(n_param):
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
                        brn_names.append(nb)
                 #      brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
                 #       nb = "b%d_CAT%d"%(p,cat_num)
                        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                                   
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                                   
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
                #        nb = "b%d_CAT%d"%(p,cat_num)
                        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                            
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                            $
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                            $
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                            $
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                            $
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                            $
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                                   
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
                #        nb = "b%d_CAT%d"%(p,cat_num)
                        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                                   
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                                   
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                                   
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                                   
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                                   
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                                   
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
		#	nb = "b%d_%s_CAT0"%(p,pdf_name)
		#	nb = "b%d_CAT%d"%(p,cat_num)
			nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
			brn_names.append(nb)
		#	brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
                #       nb = "b%d_CAT%d"%(p,cat_num)
                        nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
		#	nb = "b%d_%s_CAT%d"%(p,pdf_name,cat_num)
		#	nb2 = "b%d_CAT%d"%(p,cat_num)
			nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
			brn_names.append(nb)
		#	brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
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
		#	nb = "b%d_%s_CAT0"%(p,pdf_name)
		#	nb2 = "b%d_CAT0"%(p)
			nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
		#	nb = "b%d_%s_CAT%d"%(p,pdf_name,cat_num)
			brn_names.append(nb)
		#	brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
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
                        nb = "b%d_%s_CAT%d"%(p,pdf_name,cat_num)
                #	nb2 = "b%d_CAT%d"%(p,cat_num)
                #	nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
                        brn_names.append(nb)
                #	brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
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
                        nb = "b%d_%s_CAT%d"%(p,pdf_name,cat_num)
                #       nb2 = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
                        nb = "b%d_%s_CAT%d"%(p,pdf_name,cat_num)
                #       nb2 = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
                        nb = "b%d_%s_CAT%d"%(p,pdf_name,cat_num)
                #       nb2 = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
                        brn_names.append(nb)
                #       brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
                        [pmin,pmax] = Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,pmin,pmax)
                        coeff.add(brn[nb])
                        gcs.append(brn[nb])
                formula = "%s-%s*TMath::TanH(%s*%s-%s)-%s*TMath::Power(TMath::TanH(%s*%s-%s),2)-%s*TMath::Power(TMath::TanH(%s*%s-%s),3)-%s*TMath::Power(TMath::TanH(%s*%s-%s),4)-%s*TMath::Power(TMath::TanH(%s*%s-%s),5)-%s*TMath::Power(TMath::TanH(%s*%s-%s),6)"%(brn_names[0],brn_names[1],x_name,brn_names[2],brn_names[3],brn_names[4],x_name,brn_names[2],brn_names[3],brn_names[5],x_name,brn_names[2],brn_names[3],brn_names[6],x_name,brn_names[2],brn_names[3],brn_names[7],x_name,brn_names[2],brn_names[3],brn_names[8],x_name,brn_names[2],brn_names[3])
                coeff.add(x)
                pdf = ROOT.RooGenericPdf(name,"",formula,coeff )

	elif pdf_name=="sine":            
		coeff.removeAll()
		brn = {}
		brn_names = []
		for p in xrange(n_param):   
		#	nb = "b%d_%s_CAT0"%(p,pdf_name)
		#	nb = "b%d_CAT%d"%(p,cat_num)
			nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
			brn_names.append(nb)
		#	brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)                                          
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
                #	nb = "b%d_%s_CAT0"%(p,pdf_name)
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
                #       nb = "b%d_%s_CAT0"%(p,pdf_name)
                        nb = "b%d_CAT%d"%(p,cat_num)
                #       nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
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
		#	nb = "b%d_%s_CAT0"%(p,pdf_name)
			nb2 = "b%d_CAT%d"%(p,cat_num)
			nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
			brn_names.append(nb)
		#	brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
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
		#	nb = "b%d_%s_CAT0"%(p,pdf_name)
			nb2 = "b%d_CAT0"%(p)
			nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
			brn_names.append(nb)
		#	brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
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
		#	nb = "b%d_%s_CAT0"%(p,pdf_name)
			nb2 = "b%d_CAT0"%(p)
			nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
			brn_names.append(nb)
		#	brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
			[pmin,pmax] = Parameters[pdf_name][nb2] 
			brn[nb] = RooRealVar(nb,nb,pmin,pmax)
			coeff.add(brn[nb])
			gcs.append(brn[nb])
		formula = "TMath::Erf((%s-%s)/%s)*(%s*TMath::Power(%s,3)+%s*TMath::Power(%s,2)+%s*%s+1.)"%(x_name,brn_names[0],brn_names[1],brn_names[2],x_name,brn_names[3],x_name,brn_names[4],x_name)
		coeff.add(x)
		pdf = ROOT.RooGenericPdf(name,"",formula,coeff )
	if pdf_name=="Pol2":            
		coeff.removeAll()
		brn = {}
		brn_names = []
		for p in xrange(n_param):   
		#	nb = "b%d_%s_CAT0"%(p,pdf_name)
			nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
			brn_names.append(nb)
			#brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
			brn[nb] = RooRealVar(nb,nb,-1,1)
			coeff.add(brn[nb])
			gcs.append(brn[nb])
		pdf = ROOT.RooChebychev(name,"",x,coeff )
	if pdf_name=="Pol3":            
		coeff.removeAll()
		brn = {}
		brn_names = []
		for p in xrange(n_param):   
		#	nb = "b%d_%s_CAT%d"%(p,pdf_name,cat_num)
	#		nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
			nb = "b%d_Pol3_CAT%d"%(p,cat_num) #BG change
			brn_names.append(nb)
			pmin = Parameters[pdf_name][nb]
		#	brn[nb] = RooRealVar(nb,nb,0.48,0,10.)
			brn[nb] = RooRealVar(nb,nb,-1,1)
			coeff.add(brn[nb])
			gcs.append(brn[nb])
		pdf = ROOT.RooChebychev(name,"",x,coeff )
	if pdf_name=="Pol4":            
		coeff.removeAll()
		brn = {}
		brn_names = []
		for p in xrange(n_param):   
		#	nb = "b%d_%s_CAT%d"%(p,pdf_name,cat_num)
		#	nb = "b%d_%s_CAT0"%(p,pdf_name)
			nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
			brn_names.append(nb)
		#	brn[nb] = RooRealVar(nb,nb,0.38,0,10.)
			pmin = Parameters[pdf_name][nb]
			brn[nb] = RooRealVar(nb,nb,-1,1)
			coeff.add(brn[nb])
			gcs.append(brn[nb])
		pdf = ROOT.RooChebychev(name,"",x,coeff )
	if pdf_name=="Pol5":            
		coeff.removeAll()
		brn = {}
		brn_names = []
		for p in xrange(n_param):   
		#	nb = "b%d_%s_CAT0"%(p,pdf_name)
			nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
		#	nb = "b%d_%s_CAT%d"%(p,pdf_name,cat_num)
			brn_names.append(nb)
		#	brn[nb] = RooRealVar(nb,nb,0.5,0,10.)
			pmin = 0#Parameters[pdf_name][nb]
                        brn[nb] = RooRealVar(nb,nb,-1,1) #BG changed 10->1
			coeff.add(brn[nb])
			print "test",coeff
			gcs.append(brn[nb])
		#pdf = ROOT.RooBernstein(name,"",x,coeff )
		pdf = ROOT.RooChebychev(name,"",x,coeff )
		#pdf = ROOT.RooBernstein(name,"",(x-80.0)/120.0,coeff )
#(x-min)/(mass-min) so range b.t. 0 and 1
                print brn, brn_names
	if pdf_name=="Pol6":            
		coeff.removeAll()
		brn = {}
		brn_names = []
		for p in xrange(n_param):   
		#	nb = "b%d_%s_CAT0"%(p,pdf_name)
			#nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num) #BG comment
			nb = "b%d_Pol6_CAT%d"%(p,cat_num)
		#	nb = "b%d_%s_CAT%d"%(p,pdf_name,cat_num)
			brn_names.append(nb)
		#	brn[nb] = RooRealVar(nb,nb,0.32,0,10.)
			pmin = Parameters[pdf_name][nb]
                        #brn[nb] = RooRealVar(nb,nb,pmin,10)
                        brn[nb] = RooRealVar(nb,nb,-1,1)
			coeff.add(brn[nb])
			gcs.append(brn[nb])
		pdf = ROOT.RooChebychev(name,"",x,coeff )
	if pdf_name=="Pol7":            
		coeff.removeAll()
		brn = {}
		brn_names = []
		for p in xrange(n_param):   
		#	nb = "b%d_%s_CAT0"%(p,pdf_name)
			nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
		#	nb = "b%d_%s_CAT%d"%(p,pdf_name,cat_num)
			brn_names.append(nb)
		#	brn[nb] = RooRealVar(nb,nb,0.30,0,10.)
			pmin = Parameters[pdf_name][nb]
                        #brn[nb] = RooRealVar(nb,nb,pmin,10)
                        brn[nb] = RooRealVar(nb,nb,-1,1)
			coeff.add(brn[nb])
			gcs.append(brn[nb])
		pdf = ROOT.RooChebychev(name,"",x,coeff )
	if pdf_name=="Pol8":            
		coeff.removeAll()
		brn = {}
		brn_names = []
		for p in xrange(n_param):   
		#	nb = "b%d_%s_CAT0"%(p,pdf_name)
	#		nb = "b%d_sel%s_CAT%d"%(p,selection,cat_num)
			nb = "b%d_Pol8_CAT%d"%(p,cat_num)  #BG change
		#	nb = "b%d_%s_CAT%d"%(p,pdf_name,cat_num)
			brn_names.append(nb)
		#	brn[nb] = RooRealVar(nb,nb,3.54,0,10.)
			pmin = Parameters[pdf_name][nb]
                        #brn[nb] = RooRealVar(nb,nb,pmin,10)
                        brn[nb] = RooRealVar(nb,nb,-1,1)
			coeff.add(brn[nb])
			gcs.append(brn[nb])
		pdf = ROOT.RooChebychev(name,"",x,coeff )


	pdf.Print()
	coeff.Print()
	return [pdf, coeff]
