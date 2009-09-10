/*
 * main.cpp
 *
 *  Created on: Aug 14, 2009
 *      Author: Stephen A. Smith
 */

#include <ctime>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <string.h>
#include <map>
#include <math.h>
using namespace std;

#include "RateMatrixUtils.h"
#include "BioGeoTreeTools.h"
#include "RateModel.h"
#include "BioGeoTree.h"
#include "OptimizeBioGeo.h"
#include "OptimizeBioGeoPowell.h"
#include "InputReader.h"
#include "Utils.h"

#include "expm.h"

#include <Phyl/TreeTemplate.h>
#include <Phyl/TreeTemplateTools.h>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/Newick.h>
#include <Phyl/Tree.h>
#include <Utils/BppVector.h>
using namespace bpp;

int main(int argc, char* argv[]){
	if(argc != 2){
		cout << "you need more arguments." << endl;
		cout << "usage: lagrange configfile" << endl;
		exit(0);
	}else{
		string treefile;
		string datafile;
		string ratematrixfile;
		string logfile;
		int maxareas=0;
		vector<double> periods;
		map<string,vector<string> > mrcas;
		map<string,vector<int> > fixnodewithmrca;
		vector<vector<int> > excludedists;
		vector<vector<int> > includedists;
		vector<string> areanames;
		map<string,int> areanamemap;
		map<int,string> areanamemaprev;
		vector<string> ancstates;//string should be the mrca or if it is just _all_
		//then everything will be computed
		bool treecolors;
		vector<string> areacolors;
		vector<string> fossilmrca;
		vector<string> fossiltype;
		vector<string> fossilarea;
		vector<double> fossilage;//0's for N type and # for B type

		bool marginal = true; // false means joint
		bool splits = false;
		bool states = false;
		int numthreads = 0;
		bool sparse = false;

		BioGeoTreeTools tt;

		//read file
		ifstream ifs(argv[1]);
		string line;
		while(getline(ifs,line)){
			if(line.size()>0){
				if(&line[0]!="#"){
					vector<string> tokens;
					string del("=");
					tokens.clear();
					Tokenize(line, tokens, del);
					for(unsigned int j=0;j<tokens.size();j++){
						TrimSpaces(tokens[j]);
					}
					if(!strcmp(tokens[0].c_str(), "treefile")){
						treefile = tokens[1];
					}else if(!strcmp(tokens[0].c_str(),  "datafile")){
						datafile = tokens[1];
					}else if(!strcmp(tokens[0].c_str(),  "ratematrix")){
						ratematrixfile = tokens[1];
						if(ratematrixfile == "d" || ratematrixfile == "D"){
							ratematrixfile = "";
						}
					}else if(!strcmp(tokens[0].c_str(), "areanames")){
						vector<string> searchtokens;
						Tokenize(tokens[1], searchtokens, ", 	");
						for(unsigned int j=0;j<searchtokens.size();j++){
							TrimSpaces(searchtokens[j]);
						}
						areanames = searchtokens;
					}else if(!strcmp(tokens[0].c_str(), "fixnode")){
						vector<string> searchtokens;
						Tokenize(tokens[1], searchtokens, ", 	");
						for(unsigned int j=0;j<searchtokens.size();j++){
							TrimSpaces(searchtokens[j]);
						}
						vector<int> dist;
						for(unsigned int j=0;j<searchtokens[1].size();j++){
							char c = (searchtokens[1].c_str())[j];
							dist.push_back(atoi(&c));
						}
						fixnodewithmrca[searchtokens[0]] = dist;
					}else if(!strcmp(tokens[0].c_str(), "excludedists")){
						vector<string> searchtokens;
						Tokenize(tokens[1], searchtokens, ", 	");
						for(unsigned int j=0;j<searchtokens.size();j++){
							TrimSpaces(searchtokens[j]);
						}
						for(unsigned int j=0;j<searchtokens.size();j++){
							vector<int> dist;
							for(unsigned int k=0;k<searchtokens[j].size();k++){
								char c = (searchtokens[j].c_str())[k];
								dist.push_back(atoi(&c));
							}
							excludedists.push_back(dist);
						}
					}else if(!strcmp(tokens[0].c_str(), "includedists")){
						vector<string> searchtokens;
						Tokenize(tokens[1], searchtokens, ", 	");
						for(unsigned int j=0;j<searchtokens.size();j++){
							TrimSpaces(searchtokens[j]);
						}
						if(searchtokens[0].size()==1){
							maxareas = atoi(searchtokens[0].c_str());
						}else{
							for(unsigned int j=0;j<searchtokens.size();j++){
								vector<int> dist;
								for(unsigned int k=0;k<searchtokens[j].size();k++){
									char c = (searchtokens[j].c_str())[k];
									dist.push_back(atoi(&c));
								}
								includedists.push_back(dist);
							}
						}
					}else if(!strcmp(tokens[0].c_str(), "areacolors")){
						vector<string> searchtokens;
						Tokenize(tokens[1], searchtokens, ", 	");
						for(unsigned int j=0;j<searchtokens.size();j++){
							TrimSpaces(searchtokens[j]);
						}
						areacolors = searchtokens;
					}else if(!strcmp(tokens[0].c_str(), "periods")){
						vector<string> searchtokens;
						Tokenize(tokens[1], searchtokens, ", 	");
						for(unsigned int j=0;j<searchtokens.size();j++){
							TrimSpaces(searchtokens[j]);
							periods.push_back(atof(searchtokens[j].c_str()));
						}
					}else if(!strcmp(tokens[0].c_str(),  "treecolors")){
						treecolors = true;
					}else if(!strcmp(tokens[0].c_str(), "mrca")){
						vector<string> searchtokens;
						Tokenize(tokens[1], searchtokens, ", 	");
						for(unsigned int j=0;j<searchtokens.size();j++){
							TrimSpaces(searchtokens[j]);
						}
						vector<string> mrc;
						for(unsigned int j=1;j<searchtokens.size();j++){
							mrc.push_back(searchtokens[j]);
						}
						mrcas[searchtokens[0]] = mrc;
					}else if(!strcmp(tokens[0].c_str(), "ancstate")){
						vector<string> searchtokens;
						Tokenize(tokens[1], searchtokens, ", 	");
						for(unsigned int j=0;j<searchtokens.size();j++){
							TrimSpaces(searchtokens[j]);
						}
						ancstates.push_back(searchtokens[0]);
					}else if(!strcmp(tokens[0].c_str(), "fossil")){
						vector<string> searchtokens;
						Tokenize(tokens[1], searchtokens, ", 	");
						for(unsigned int j=0;j<searchtokens.size();j++){
							TrimSpaces(searchtokens[j]);
						}
						fossiltype.push_back(searchtokens[0]);
						fossilmrca.push_back(searchtokens[1]);
						fossilarea.push_back(searchtokens[2]);
						if(searchtokens.size()>3){
							fossilage.push_back(atof(searchtokens[3].c_str()));
						}else{
							fossilage.push_back(0.0);
						}
					}else if(!strcmp(tokens[0].c_str(),  "calctype")){
						string calctype = tokens[1];
						if(calctype.compare("m") != 0 && calctype.compare("M") != 0){
							marginal = false;
						}
					}else if(!strcmp(tokens[0].c_str(),  "report")){
						if(tokens[1].compare("split") != 0){
							splits = false;
						}
					}else if(!strcmp(tokens[0].c_str(),  "sparse")){
						sparse = true;
					}else if(!strcmp(tokens[0].c_str(),  "splits")){
						splits = true;
					}else if(!strcmp(tokens[0].c_str(),  "states")){
						states = true;
					}else if(!strcmp(tokens[0].c_str(),  "numthreads")){
						numthreads = atoi(tokens[1].c_str());
					}
				}
			}
		}
		ifs.close();
		/*
		 * after reading the input file
		 */
		InputReader ir;
		vector<TreeTemplate<Node> *> intrees= ir.readMultipleTreeFile(treefile);
		map<string,vector<int> > data = ir.readStandardInputData(datafile);
		ir.checkData(data,intrees);

		/*
		 * read area names
		 */
		if(areanames.size() > 0){
			for(unsigned int i=0;i<areanames.size();i++){
				areanamemap[areanames[i]] = i;
				areanamemaprev[i] = areanames[i];
			}
		}else{
			for(int i=0;i<ir.nareas;i++){
				std::ostringstream osstream;
				osstream << i;
				std::string string_x = osstream.str();
				areanamemap[string_x] = i;
				areanamemaprev[i] = string_x;
			}
		}
		/*
		 * need to figure out how to work with multiple trees best
		 */
		if(periods.size() < 1){
			periods.push_back(10000);
		}
		RateModel rm(ir.nareas,true,periods,sparse);
		if(numthreads != 0){
			rm.set_nthreads(numthreads);
			cout << "Setting the number of threads: " << numthreads << endl;
		}
		rm.setup_Dmask();
		/*
		 * if there is a ratematrixfile then it will be processed
		 */
		if(ratematrixfile != "" && ratematrixfile.size() > 0){
			cout << "Reading rate matrix file" << endl;
			vector<vector<vector<double> > > dmconfig = processRateMatrixConfigFile(ratematrixfile,ir.nareas,periods.size());
			for(unsigned int i=0;i<dmconfig.size();i++){
				for(unsigned int j=0;j<dmconfig[i].size();j++){
					for(unsigned int k=0;k<dmconfig[i][j].size();k++){
						if(dmconfig[i][j][k] != 1){
							cout << dmconfig[i][j][k];
						}else{
							cout << " . ";
						}
						cout << " ";
					}
					cout << endl;
				}cout << endl;
			}
			for(unsigned int i=0;i<dmconfig.size();i++){
				for(unsigned int j=0;j<dmconfig[i].size();j++){
					for(unsigned int k=0;k<dmconfig[i][j].size();k++){
						rm.set_Dmask_cell(i,j,k,dmconfig[i][j][k],false);
					}
				}
			}
		}
		if(includedists.size() > 0 || excludedists.size() > 0 || maxareas >= 2){
			if(excludedists.size() > 0){
				rm.setup_dists(excludedists,false);
			}else{
				if(maxareas >=2)
					includedists = generate_dists_from_num_max_areas(ir.nareas,maxareas);
				rm.setup_dists(includedists,true);
			}
		}else{
			rm.setup_dists();
		}
		rm.setup_D(0.01);
		rm.setup_E(0.01);
		rm.setup_Q();

		/*
		 * outfile for tree reconstructed states
		 */
		ofstream outTreeFile;
		ofstream outTreeKeyFile;

		for(unsigned int i =0;i<intrees.size();i++){
			BioGeoTree bgt(intrees[i],periods);
			/*
			 * record the mrcas
			 */
			map<string,int> mrcanodeint;
			map<string,vector<string> >::iterator it;
			for(it=mrcas.begin();it != mrcas.end();it++){
				vector<int> nodeIds;
				for(unsigned int k=0;k<(*it).second.size();k++){
					nodeIds.push_back(intrees[i]->getNode((*it).second[k])->getId());
				}
				mrcanodeint[(*it).first] = tt.getLastCommonAncestor(*intrees[i],nodeIds);
				cout << "Reading mrca: " << (*it).first << " " << mrcanodeint[(*it).first] << endl;
			}

			/*
			 * set fixed nodes
			 */
			map<string,vector<int> >::iterator fnit;
			for(fnit = fixnodewithmrca.begin(); fnit != fixnodewithmrca.end(); fnit++){
				vector<int> dista = (*fnit).second;
				for(unsigned int k=0;k<rm.getDists()->size();k++){
					bool isnot = true;
					for(unsigned int j=0;j<dista.size();j++){
						if(dista[j] != rm.getDists()->at(k)[j])
							isnot = false;
					}
					if(isnot == false){
						bgt.set_excluded_dist(rm.getDists()->at(k),intrees[0]->getNode(mrcanodeint[(*fnit).first]));
					}
				}
				cout << "fixing " << (*fnit).first << " = ";print_vector_int((*fnit).second);
			}


			bgt.set_default_model(&rm);
			bgt.set_tip_conditionals(data);

			/*
			 * setting up fossils
			 */
			for(unsigned int k=0;k<fossiltype.size();k++){
				if(fossiltype[k] == "n" || fossiltype[k] == "N"){
					bgt.setFossilatNodeByMRCA_id(mrcanodeint[fossilmrca[k]],areanamemap[fossilarea[k]]);
					cout << "Setting node fossil at mrca: " << fossilmrca[k] << " at area: " << fossilarea[k] << endl;
				}else{
					bgt.setFossilatBranchByMRCA_id(mrcanodeint[fossilmrca[k]],areanamemap[fossilarea[k]],fossilage[k]);
					cout << "Setting branch fossil at mrca: " << fossilmrca[k] << " at area: " << fossilarea[k] << " at age: " << fossilage[k] << endl;
				}
			}

			/*
			 * initial likelihood calculation
			 */
			cout << "initial -ln likelihood: " << -log(bgt.eval_likelihood(marginal)) <<endl;
			/*
			 * optimize likelihood
			 */
			if(sparse == false){
				cout << "Optimizing (simplex) -ln likelihood." << endl;
				OptimizeBioGeo opt(&bgt,&rm,marginal);
				vector<double> disext  = opt.optimize_global_dispersal_extinction();
				cout << "dis: " << disext[0] << " ext: " << disext[1] << endl;
				rm.setup_D(disext[0]);
				rm.setup_E(disext[1]);
				rm.setup_Q();
				bgt.update_default_model(&rm);
				bgt.set_store_p_matrices(true);
				cout << "final -ln likelihood: "<< -log(bgt.eval_likelihood(marginal)) <<endl;
				bgt.set_store_p_matrices(false);
			}else{
				cout << "Optimizing (Powell) -ln likelihood." << endl;
				OptimizeBioGeoPowell opt(&bgt,&rm,marginal);
				vector<double> disext  = opt.optimize_global_dispersal_extinction();
				cout << "dis: " << disext[0] << " ext: " << disext[1] << endl;
				rm.setup_D(disext[0]);
				rm.setup_E(disext[1]);
				rm.setup_Q();
				bgt.update_default_model(&rm);
				cout << "final -ln likelihood: "<< -log(bgt.eval_likelihood(marginal)) <<endl;
			}
			/*
			 * ancestral splits calculation
			 */
			if(ancstates.size() > 0){
				bgt.set_use_stored_matrices(true);
				bgt.prepare_ancstate_reverse();
				if(ancstates[0] == "_all_" || ancstates[0] == "_ALL_"){
					for(unsigned int j=0;j<intrees[i]->getNumberOfNodes();j++){
						if(intrees[i]->getNode(j)->isLeaf()==false){
							if(splits){
								cout << "Ancestral splits for:\t" << intrees[i]->getNode(j)->getId() <<endl;
								map<vector<int>,vector<AncSplit> > ras = bgt.calculate_ancsplit_reverse(*intrees[i]->getNode(j),marginal);
								//bgt.ancstate_calculation_all_dists(*intrees[i]->getNode(j),marginal);
								tt.summarizeSplits(intrees[i]->getNode(j),ras,areanamemaprev,&rm);
								cout << endl;
							}
							if(states){
								cout << "Ancestral states for:\t" << intrees[i]->getNode(j)->getId() <<endl;
								vector<double> rast = bgt.calculate_ancstate_reverse(*intrees[i]->getNode(j),marginal);
								tt.summarizeAncState(intrees[i]->getNode(j),rast,areanamemaprev,&rm);
								cout << endl;
							}
						}
					}
					/*
					 * key file output
					 */
					outTreeKeyFile.open((treefile+".bgkey.tre").c_str(),ios::app );
					TreeTemplateTools ttt;
					outTreeKeyFile << ttt.nodeToParenthesis(*intrees[i]->getRootNode(),true) << ";"<< endl;
					outTreeKeyFile.close();
				}else{
					for(unsigned int j=0;j<ancstates.size();j++){
						if(splits){
							cout << "Ancestral splits for: " << ancstates[j] <<endl;
							map<vector<int>,vector<AncSplit> > ras = bgt.calculate_ancsplit_reverse(*intrees[i]->getNode(mrcanodeint[ancstates[j]]),marginal);
							tt.summarizeSplits(intrees[i]->getNode(mrcanodeint[ancstates[j]]),ras,areanamemaprev,&rm);
						}
						if(states){
							cout << "Ancestral splits for: " << ancstates[j] <<endl;
							vector<double> rast = bgt.calculate_ancstate_reverse(*intrees[i]->getNode(mrcanodeint[ancstates[j]]),marginal);
							tt.summarizeAncState(intrees[i]->getNode(mrcanodeint[ancstates[j]]),rast,areanamemaprev,&rm);
						}
					}
				}
				if(splits){
					outTreeFile.open((treefile+".bgsplits.tre").c_str(),ios::app );
					TreeTemplateTools ttt;
					outTreeFile << ttt.nodeToParenthesis(*intrees[i]->getRootNode(),false,"split") << ";"<< endl;
					outTreeFile.close();
				}
				if(states){
					TreeTemplateTools ttt;
					outTreeFile.open((treefile+".bgstates.tre").c_str(),ios::app );
					outTreeFile << ttt.nodeToParenthesis(*intrees[i]->getRootNode(),false,"state") << ";"<< endl;
					outTreeFile.close();
				}
				//cout << bgt.ti/CLOCKS_PER_SEC << " secs for anc" << endl;
			}

		}
		for(unsigned int i=0;i<intrees.size();i++){
			delete intrees[i];
		}
	}

	return 0;
}
