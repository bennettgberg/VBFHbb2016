#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "skim_trees_macro_vbfhbb.h"
#include "skim_trees_macro_vbfhbb.C"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>

using namespace std;

int main(int argc, char* argv[]){
		TString tree_number = std::string(argv[1]); 	
		TString output_dir = std::string(argv[2]);
		stringstream ss;
		ss.str(argv[3]); //list of files
		vector<string> list_of_files; //comma separated string and then parse before the loop for AddFile ;
		while( ss.good() )
		{
   		string substr;
    		getline( ss, substr, ',' );
    		list_of_files.push_back( substr );
		}
	
		TFile *f = TFile::Open(list_of_files[0].c_str());
		if (f!=0){	
			skim_trees_macro_vbfhbb *c = new skim_trees_macro_vbfhbb(list_of_files[0]); 
			for (int i=1;i<list_of_files.size();i++){
				c->AddFile(list_of_files[i]); // also adds content of histograms
			}
			c->Loop(output_dir+"/skimmed_tree_"+tree_number+".root");
		}
	return 0;
}
