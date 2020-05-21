import os
import sys

#script to setup bias tests for multiple categories and/or orders of polynomial
#special name to call the new directories
spec_name = "/TFPOL1/"
#use systematics or nah
systematics = True 
#use --saveShapes in the combine call or nah
saveShapes = False
sig_strength = 1 #how much signal to expect

#which categories to setup for
for c in range(2, 7):#9):
#  if not os.path.exists("cat%d"%c):
#    os.system("mkdir cat%d"%c)
  #now run the 4 tests
#  if c in [2, 3, 6, 7]: continue
#  if c == 3: continue
#Polynomial orders to generate the toys with
  for j in range(3,4):#orders[c]-1, orders[c]+1):
    #Polynomial orders to fit the toys with
    for k in range(5,6):#orders[c]-1, orders[c]+1):
#      if not abs(j - k) <= 1:
#        continue
      if os.path.exists("cat%d%sPol%d_Pol%d"%(c, spec_name, j, k)):
        #get confirmation before overwriting a directory
        overwrite = raw_input("Path cat%d%sPol%d_Pol%d exists already. Overwrite all contents? [y/n] "%(c, spec_name, j, k))
        if not overwrite in ["y", "yes", "Y", "Yes", "YES"]:
            sys.exit("Terminating program.") 
      if not os.path.exists("cat%d%s"%(c, spec_name)):
        os.system("mkdir cat%d%s"%(c, spec_name))
#must do all in one os.system call or else resets back to current dir.
      os.system("mkdir cat%d%sPol%d_Pol%d"%(c, spec_name, j, k) + "; cd cat%d%sPol%d_Pol%d"%(c, spec_name, j, k) + "; ~/public/make_scripts.py --njobs 500 --ncalls 25 --ntoys 1 --gen_func Pol%d --fit_func Pol%d --seed_params %d,%d,%d,%d,%d --CATS %d,%d --dir %s %s -e %f %s" %(j,k,c*c*j+k,j*k+c,16*j*j+k*k+1739-37*c,678+c*k-j,43+k*j, c, c, spec_name, " --saveShapes 1" if saveShapes else "", sig_strength, "-S 0" if not systematics else ""))
