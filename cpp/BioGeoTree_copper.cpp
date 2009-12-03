/*
 * BioGeoTree.cpp
 *
 *  Created on: Aug 15, 2009
 *      Author: Stephen A. Smith
 */
#include <vector>
#include <string>
#include <algorithm>
#include <ctime>
#include <functional>
#include <numeric>
#include <iostream>
using namespace std;

#include "BioGeoTree_copper.h"
#include "BioGeoTreeTools_copper.h"
#include "BranchSegment_copper.h"
#include "RateMatrixUtils.h"
#include "RateModel.h"
#include "AncSplit.h"

#include "tree.h"
#include "node.h"
#include "vector_node_object.h"

namespace {
	inline double MAX(const double &a, const double &b)
	        {return b > a ? (b) : double(a);}
}

BioGeoTree_copper::BioGeoTree_copper(Tree * tr, vector<double> ps){
	seg = "segments";
	age = "age";
	dc = "dist_conditionals";
	en = "excluded_dists";
	andc = "anc_dist_conditionals";
	/*
        reverse bit
	 */
	revB  = "revB";
	/*
	end of the reverse bits
	*/
	store_p_matrices = false;
	use_stored_matrices = false;
	tree = tr;

	periods = ps;
	/*
	 * initialize each node with segments
	 */
	cout << "initializing nodes..." << endl;
	for(int i=0;i<tree->getExternalNodeCount();i++){
		VectorNodeObject<BranchSegment> * segs = new VectorNodeObject<BranchSegment>();
		tree->getExternalNode(i)->assocObject(seg,*segs);
		VectorNodeObject<vector<int> > * ens = new VectorNodeObject<vector<int> >();
		tree->getExternalNode(i)->assocObject(en,*ens);
	}
	for(int i=0;i<tree->getInternalNodeCount();i++){
		VectorNodeObject<BranchSegment> * segs = new VectorNodeObject<BranchSegment>();
		tree->getInternalNode(i)->assocObject(seg,*segs);
		VectorNodeObject<vector<int> > * ens = new VectorNodeObject<vector<int> >();
		tree->getInternalNode(i)->assocObject(en,*ens);
	}


	/*
	 * initialize the actual branch segments for each node
	 */
	tree->setHeightFromTipToNodes();
	cout << "initializing branch segments..." << endl;
	for(int i=0;i<tree->getExternalNodeCount();i++){
		if (tree->getExternalNode(i)->hasParent()){
			vector<double> pers(periods);
			double anc = tree->getExternalNode(i)->getParent()->getHeight();
			double des = tree->getExternalNode(i)->getHeight();
			//assert anc > des:q
			double t = des;
			if (pers.size() > 0){
				for(unsigned int j=0;j<pers.size();j++){
					double s = 0;
					if(pers.size() == 1)
						s = pers[0];
					for (unsigned int k=0;k<j+1;k++){
						s += pers[k];
					}
					if (t < s){
						double duration = min(s-t,anc-t);
						if (duration > 0){
							BranchSegment tseg = BranchSegment(duration,j);
							((VectorNodeObject<BranchSegment>*) tree->getExternalNode(i)->getObject(seg))->push_back(tseg);
						}
						//t += pers[j];
						t += duration; // TODO: make sure that this is all working
					}
					if (t > anc || pers[j] > t){
						break;
					}
				}
			}else{
				BranchSegment tseg = BranchSegment(tree->getExternalNode(i)->getBL(),0);
				((VectorNodeObject<BranchSegment>*) tree->getExternalNode(i)->getObject(seg))->push_back(tseg);
			}
		}
	}
	for(int i=0;i<tree->getInternalNodeCount();i++){
		if (tree->getInternalNode(i)->hasParent()){
			vector<double> pers(periods);
			double anc = tree->getInternalNode(i)->getParent()->getHeight();
			double des = tree->getInternalNode(i)->getHeight();
			//assert anc > des:q
			double t = des;
			if (pers.size() > 0){
				for(unsigned int j=0;j<pers.size();j++){
					double s = 0;
					if(pers.size() == 1)
						s = pers[0];
					for (unsigned int k=0;k<j+1;k++){
						s += pers[k];
					}
					if (t < s){
						double duration = min(s-t,anc-t);
						if (duration > 0){
							BranchSegment tseg = BranchSegment(duration,j);
							((VectorNodeObject<BranchSegment>*) tree->getInternalNode(i)->getObject(seg))->push_back(tseg);
						}
						//t += pers[j];
						t += duration; // TODO: make sure that this is all working
					}
					if (t > anc || pers[j] > t){
						break;
					}
				}
			}else{
				BranchSegment tseg = BranchSegment(tree->getInternalNode(i)->getBL(),0);
				((VectorNodeObject<BranchSegment>*) tree->getInternalNode(i)->getObject(seg))->push_back(tseg);
			}
		}
	}
}

void BioGeoTree_copper::set_store_p_matrices(bool i){
	store_p_matrices = i;
}

void BioGeoTree_copper::set_use_stored_matrices(bool i){
	use_stored_matrices = i;
}

void BioGeoTree_copper::set_default_model(RateModel * mod){
	rootratemodel = mod;
	for(int i=0;i<tree->getExternalNodeCount();i++){
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) tree->getExternalNode(i)->getObject(seg));
		for(unsigned int j=0;j<tsegs->size();j++){
			tsegs->at(j).setModel(mod);
			VectorNodeObject<double> * distconds = new VectorNodeObject<double> (rootratemodel->getDists()->size(), 0);
			tsegs->at(j).distconds = distconds;
			VectorNodeObject<double> * ancdistconds = new VectorNodeObject<double> (rootratemodel->getDists()->size(), 0);
			tsegs->at(j).ancdistconds = ancdistconds;
		}
	}
	for(int i=0;i<tree->getInternalNodeCount();i++){
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) tree->getInternalNode(i)->getObject(seg));
		for(unsigned int j=0;j<tsegs->size();j++){
			tsegs->at(j).setModel(mod);
			VectorNodeObject<double> * distconds = new VectorNodeObject<double> (rootratemodel->getDists()->size(), 0);
			tsegs->at(j).distconds = distconds;
			VectorNodeObject<double> * ancdistconds = new VectorNodeObject<double> (rootratemodel->getDists()->size(), 0);
			tsegs->at(j).ancdistconds = ancdistconds;
		}
	}
	VectorNodeObject<double> * distconds = new VectorNodeObject<double> (rootratemodel->getDists()->size(), 0);
	tree->getRoot()->assocObject(dc,*distconds);
	VectorNodeObject<double> * ancdistconds = new VectorNodeObject<double> (rootratemodel->getDists()->size(), 0);
	tree->getRoot()->assocObject(andc,*ancdistconds);
}

void BioGeoTree_copper::update_default_model(RateModel * mod){
	rootratemodel = mod;
	for(int i=0;i<tree->getExternalNodeCount();i++){
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) tree->getExternalNode(i)->getObject(seg));
		for(unsigned int j=0;j<tsegs->size();j++){
			tsegs->at(j).setModel(mod);
		}
	}
	for(int i=0;i<tree->getInternalNodeCount();i++){
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) tree->getInternalNode(i)->getObject(seg));
		for(unsigned int j=0;j<tsegs->size();j++){
			tsegs->at(j).setModel(mod);
		}
	}
}

void BioGeoTree_copper::set_tip_conditionals(map<string,vector<int> > distrib_data){
	int numofleaves = tree->getExternalNodeCount();
	for(int i=0;i<numofleaves;i++){
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) tree->getExternalNode(i)->getObject(seg));
		RateModel * mod = tsegs->at(0).getModel();
		//int ind1 = get_vector_int_index_from_multi_vector_int(
		//		&distrib_data[tree->getLeaves().at(i)->getName()],mod->getDists());
		int ind1 = get_vector_int_index_from_multi_vector_int(
						&distrib_data[tree->getExternalNode(i)->getName()],mod->getDists());
		tsegs->at(0).distconds->at(ind1) = 1.0;
	}
}

void BioGeoTree_copper::set_excluded_dist(vector<int> ind,Node * node){
	((VectorNodeObject<vector<int> >*) node->getObject(en))->push_back(ind);
}

/*
 * **************************************************
 *
 *
 *
 * **************************************************
 */

double BioGeoTree_copper::eval_likelihood(bool marginal){
	if( rootratemodel->sparse == true){
		columns = new vector<int>(rootratemodel->getDists()->size());
		whichcolumns = new vector<int>();
	}
	ancdist_conditional_lh(*tree->getRoot(),marginal);
	if( rootratemodel->sparse == true){
		delete columns;
		delete whichcolumns;
	}
	return calculate_vector_double_sum(*
			(VectorNodeObject<double>*) tree->getRoot()->getObject(dc));

}

VectorNodeObject<double> BioGeoTree_copper::conditionals(Node & node, bool marginal,bool sparse){
	VectorNodeObject<double> distconds;
	VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getObject(seg));
	distconds = *tsegs->at(0).distconds;
	for(unsigned int i=0;i<tsegs->size();i++){
		for(unsigned int j=0;j<distconds.size();j++){
				tsegs->at(i).distconds->at(j) = distconds.at(j);
		}
		RateModel * rm = tsegs->at(i).getModel();
		VectorNodeObject<double> * v = new VectorNodeObject<double> (rootratemodel->getDists()->size(), 0);
		vector<int> distrange;
		if(tsegs->at(i).get_start_dist_int() != -666){
			int ind1 = tsegs->at(i).get_start_dist_int();
			distrange.push_back(ind1);
		}else if(tsegs->at(i).getFossilAreas().size()>0){
			for(unsigned int j=0;j<rootratemodel->getDists()->size();j++){
				distrange.push_back(j);
			}
			for(unsigned int k=0;k<distrange.size();k++){
				bool flag = true;
				for(unsigned int x = 0;x<tsegs->at(i).getFossilAreas().size();x++){
					if (tsegs->at(i).getFossilAreas()[x] == 1 && distrange.at(x) == 0){
						flag = false;
					}
				}
				if(flag == true){
                    distrange.erase(distrange.begin()+k);
				}
			}
		}else{
			for(unsigned int j=0;j<rootratemodel->getDists()->size();j++){
				distrange.push_back(j);
			}
		}
		/*
		 * marginal
		 */
		if(marginal == true){
			if(sparse == false){
				vector<vector<double > > p;
				if(use_stored_matrices == false){
					p= rm->setup_fortran_P(tsegs->at(i).getPeriod(),tsegs->at(i).getDuration(),
																 store_p_matrices);
				}else{
					p = rm->stored_p_matrices[tsegs->at(i).getPeriod()][tsegs->at(i).getDuration()];
				}
				for(unsigned int j=0;j<distrange.size();j++){
					for(unsigned int k=0;k<distconds.size();k++){
						v->at(distrange[j]) += (distconds.at(k)*p[distrange[j]][k]);
					}
				}
			}else{//sparse
				/*
				 testing pthread version
				 */
				if(rm->get_nthreads() > 0){
					vector<vector<double > > p = rm->setup_pthread_sparse_P(tsegs->at(i).getPeriod(),tsegs->at(i).getDuration(),*whichcolumns);
					for(unsigned int j=0;j<distrange.size();j++){
						for(unsigned int k=0;k<distconds.size();k++){
							v->at(distrange[j]) += (distconds.at(k)*p[distrange[j]][k]);
						}
					}
				}else{
					for(unsigned int j=0;j<distrange.size();j++){
						bool inthere = false;
						if(columns->at(distrange[j]) == 1)
							inthere = true;
						vector<double > p;
						if(inthere == true){
							p = rm->setup_sparse_single_column_P(tsegs->at(i).getPeriod(),tsegs->at(i).getDuration(),distrange[j]);
						}else{
							p = vector<double>(distconds.size(),0);
						}
						for(unsigned int k=0;k<distconds.size();k++){
							v->at(distrange[j]) += (distconds.at(k)*p[k]);
						}
					}
				}
			}
		}
		/*
		 * joint reconstruction
		 * NOT FINISHED YET -- DONT USE
		 */
		else{
			if(sparse == false){
				vector<vector<double > > p = rm->setup_fortran_P(tsegs->at(i).getPeriod(),tsegs->at(i).getDuration(),store_p_matrices);
				for(unsigned int j=0;j<distrange.size();j++){
					double maxnum = 0;
					for(unsigned int k=0;k<distconds.size();k++){
						maxnum = MAX((distconds.at(k)*p[distrange[j]][k]),maxnum);
					}
					v->at(distrange[j]) = maxnum;
				}
			}else{//sparse

			}
		}
		for(unsigned int j=0;j<distconds.size();j++){
			distconds[j] = v->at(j);
		}
		delete v;
	}
	/*
	 * if store is true we want to store the conditionals for each node
	 * for possible use in ancestral state reconstruction
	 */
	if(store_p_matrices == true){
		tsegs->at(0).alphas = distconds;
	}
	return distconds;
}

void BioGeoTree_copper::ancdist_conditional_lh(Node & node, bool marginal){
	VectorNodeObject<double> distconds(rootratemodel->getDists()->size(), 0);
	if (node.isExternal()==false){//is not a tip
		Node * c1 = &node.getChild(0);
		Node * c2 = &node.getChild(1);
		RateModel * model;
		if(node.hasParent()==true){
			VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getObject(seg));
			model = tsegs->at(0).getModel();
		}else{
			model = rootratemodel;
		}
		ancdist_conditional_lh(*c1,marginal);
		ancdist_conditional_lh(*c2,marginal);
		bool sparse = rootratemodel->sparse;
		VectorNodeObject<double> v1;
		VectorNodeObject<double> v2;
		if(sparse == true){
			//getcolumns
			VectorNodeObject<BranchSegment>* c1tsegs = ((VectorNodeObject<BranchSegment>*) c1->getObject(seg));
			VectorNodeObject<BranchSegment>* c2tsegs = ((VectorNodeObject<BranchSegment>*) c2->getObject(seg));
			vector<int> lcols = get_columns_for_sparse(*c1tsegs->at(0).distconds,rootratemodel);
			vector<int> rcols = get_columns_for_sparse(*c2tsegs->at(0).distconds,rootratemodel);
			whichcolumns->clear();
			for(unsigned int i=0;i<lcols.size();i++){
				if(lcols[i]==1 || rcols[i] ==1){
					columns->at(i)=1;
					if(i!=0 && count(whichcolumns->begin(),whichcolumns->end(),i) == 0)
						whichcolumns->push_back(i);
				}else{
					columns->at(i)=0;
				}
			}
			if(calculate_vector_int_sum(columns)==0){
				for(unsigned int i=0;i<lcols.size();i++){
					columns->at(i)=1;
				}
			}
			columns->at(0) = 0;
		}
		v1 =conditionals(*c1,marginal,sparse);
		v2 =conditionals(*c2,marginal,sparse);

		vector<vector<int> > * dists = rootratemodel->getDists();
		vector<int> leftdists;
		vector<int> rightdists;
		double weight;
		//cl1 = clock();
		for (unsigned int i=0;i<dists->size();i++){
			if(accumulate(dists->at(i).begin(),dists->at(i).end(),0) > 0){
				double lh = 0.0;
				VectorNodeObject<vector<int> >* exdist = ((VectorNodeObject<vector<int> >*) node.getObject(en));
				int cou = count(exdist->begin(),exdist->end(),dists->at(i));
				if(cou == 0){
					iter_ancsplits_just_int(rootratemodel,dists->at(i),leftdists,rightdists,weight);
					for (unsigned int j=0;j<leftdists.size();j++){
						int ind1 = leftdists[j];
						int ind2 = rightdists[j];
						double lh_part = v1.at(ind1)*v2.at(ind2);
						lh += (lh_part * weight);
					}
				}
				distconds.at(i)= lh;
			}
		}
		///cl2 = clock();
		//ti += cl2-cl1;
	}else{
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getObject(seg));
		distconds = *tsegs->at(0).distconds;
	}
	if(node.hasParent() == true){
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getObject(seg));
		for(unsigned int i=0;i<distconds.size();i++){
			tsegs->at(0).distconds->at(i) = distconds.at(i);
		}
	}else{
		for(unsigned int i=0;i<distconds.size();i++){
			((VectorNodeObject<double>*)node.getObject(dc))->at(i) = distconds.at(i);
		}
	}
}

/*
 * ********************************************
 *
 * adds fossils either at the node or along a branch
 *
 * ********************************************
 */
void BioGeoTree_copper::setFossilatNodeByMRCA(vector<string> nodeNames, int fossilarea){
	Node * mrca = tree->getMRCA(nodeNames);
	vector<vector<int> > * dists = rootratemodel->getDists();
	for(unsigned int i=0;i<dists->size();i++){
		if(dists->at(i).at(fossilarea) == 0){
			VectorNodeObject<vector<int> > * exd = ((VectorNodeObject<vector<int> > *) mrca->getObject(en));
			exd->push_back(dists->at(i));
		}
	}
}
void BioGeoTree_copper::setFossilatNodeByMRCA_id(int id, int fossilarea){
	/*vector<vector<int> > * dists = rootratemodel->getDists();
	for(unsigned int i=0;i<dists->size();i++){
		if(dists->at(i).at(fossilarea) == 0){
			VectorNodeObject<vector<int> > * exd = ((VectorNodeObject<vector<int> > *) tree->getObject(id,en));
			exd->push_back(dists->at(i));
		}
	}*/
}
void BioGeoTree_copper::setFossilatBranchByMRCA(vector<string> nodeNames, int fossilarea, double age){
	Node * mrca = tree->getMRCA(nodeNames);
	VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) mrca->getObject(seg));
	double startage = mrca->getHeight();
	for(unsigned int i=0;i<tsegs->size();i++){
		if(age > startage && age < (startage+tsegs->at(i).getDuration())){
			tsegs->at(i).setFossilArea(fossilarea);
		}
		startage += tsegs->at(i).getDuration();
	}
}
void BioGeoTree_copper::setFossilatBranchByMRCA_id(int id, int fossilarea, double age){
	/*TreeTemplateTools * ttt = new TreeTemplateTools();
	//VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) tree->getNode(id)->getObject(seg));
	//double startage = ttt->getHeight(*tree->getNode(id));
	VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) tree_get_node_from_id[id]->getObject(seg));
	double startage = ttt->getHeight(*tree_get_node_from_id[id]);

	for(unsigned int i=0;i<tsegs->size();i++){
		if(age > startage && age < (startage+tsegs->at(i).getDuration())){
			tsegs->at(i).setFossilArea(fossilarea);
		}
		startage += tsegs->at(i).getDuration();
	}
	delete ttt;*/
}


/************************************************************
 forward and reverse stuff
 ************************************************************/
//add joint
void BioGeoTree_copper::prepare_ancstate_reverse(){
    reverse(*tree->getRoot());
}

/*
 * called from prepare_ancstate_reverse and that is all
 */
void BioGeoTree_copper::reverse(Node & node){
	VectorNodeObject<double> * revconds = new VectorNodeObject<double> (rootratemodel->getDists()->size(), 0);//need to delete this at some point
	if (&node == tree->getRoot()) {
		for(unsigned int i=0;i<rootratemodel->getDists()->size();i++){
			revconds->at(i) = 1.0;//prior
		}
		node.assocObject(revB,*revconds);
		for(int i = 0;i<node.getChildCount();i++){
			reverse(node.getChild(i));
		}
	}else if(node.isExternal() == false){
		//calculate A i 
		//sum over all alpha k of sister node of the parent times the priors of the speciations 
		//(weights) times B of parent j
		VectorNodeObject<double> * parrev = ((VectorNodeObject<double>*)node.getParent()->getObject(revB));
		VectorNodeObject<double> sisdistconds;
		if(&node.getParent()->getChild(0) != &node){
			VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getParent()->getChild(0).getObject(seg));
			sisdistconds = tsegs->at(0).alphas;
		}else{
			VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getParent()->getChild(1).getObject(seg));
			sisdistconds = tsegs->at(0).alphas;
		}
		vector<vector<int> > * dists = rootratemodel->getDists();
		vector<int> leftdists;
		vector<int> rightdists;
		double weight;
		//cl1 = clock();
		VectorNodeObject<double> tempA (rootratemodel->getDists()->size(),0);
		for (unsigned int i = 0; i < dists->size(); i++) {
			if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
				VectorNodeObject<vector<int> >* exdist =
						((VectorNodeObject<vector<int> >*) node.getObject(en));
				int cou = count(exdist->begin(), exdist->end(), dists->at(i));
				if (cou == 0) {
					iter_ancsplits_just_int(rootratemodel, dists->at(i),
							leftdists, rightdists, weight);
					//root has i, curnode has left, sister of cur has right
					for (unsigned int j = 0; j < leftdists.size(); j++) {
						int ind1 = leftdists[j];
						int ind2 = rightdists[j];
						tempA[ind1] += (sisdistconds.at(ind2)*weight*parrev->at(i));
					}
				}
			}
		}

		//now calculate node B
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getObject(seg));
		for(unsigned int k=0;k<tsegs->size();k++){
			RateModel * rm = tsegs->at(k).getModel();
			vector<vector<double > > * p = &rm->stored_p_matrices[tsegs->at(k).getPeriod()][tsegs->at(k).getDuration()];
			for(unsigned int j=0;j < dists->size();j++){
				if(accumulate(dists->at(j).begin(), dists->at(j).end(), 0) > 0){
					for (unsigned int i = 0; i < dists->size(); i++) {
						if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
							revconds->at(j) += tempA[i]*((*p)[i][j]);
						}
					}
				}
			}
		}
		node.assocObject(revB,*revconds);
		for(int i = 0;i<node.getChildCount();i++){
			reverse(node.getChild(i));
		}
	}//there should be no else
}

/*
 * calculates the most likely split (not state) -- the traditional result for lagrange
 */

map<vector<int>,vector<AncSplit> > BioGeoTree_copper::calculate_ancsplit_reverse(Node & node,bool marg){
	VectorNodeObject<double> * Bs = (VectorNodeObject<double> *) node.getObject(revB);
	map<vector<int>,vector<AncSplit> > ret;
	for(unsigned int j=0;j<rootratemodel->getDists()->size();j++){
		vector<int> dist = rootratemodel->getDists()->at(j);
		vector<AncSplit> ans = iter_ancsplits(rootratemodel,dist);
		if (node.isExternal()==false){//is not a tip
			Node * c1 = &node.getChild(0);
			Node * c2 = &node.getChild(1);
			VectorNodeObject<BranchSegment>* tsegs1 = ((VectorNodeObject<BranchSegment>*) c1->getObject(seg));
			VectorNodeObject<BranchSegment>* tsegs2 = ((VectorNodeObject<BranchSegment>*) c2->getObject(seg));
			for (unsigned int i=0;i<ans.size();i++){
				VectorNodeObject<double> v1  =tsegs1->at(0).alphas;
				VectorNodeObject<double> v2 = tsegs2->at(0).alphas;
				double lh = (v1[ans[i].ldescdistint]*v2[ans[i].rdescdistint]*Bs->at(j)*ans[i].getWeight());
				ans[i].setLikelihood(lh);
				//cout << lh << endl;
			}
		}
		ret[dist] = ans;
	}
	return ret;
}

/*
 * calculates the ancestral area over all the possible splits
 */

vector<double> BioGeoTree_copper::calculate_ancstate_reverse(Node & node,bool marg){
	if (node.isExternal()==false){//is not a tip
		VectorNodeObject<double> * Bs = (VectorNodeObject<double> *) node.getObject(revB);
		vector<vector<int> > * dists = rootratemodel->getDists();
		vector<int> leftdists;
		vector<int> rightdists;
		double weight;
		Node * c1 = &node.getChild(0);
		Node * c2 = &node.getChild(1);
		VectorNodeObject<BranchSegment>* tsegs1 = ((VectorNodeObject<BranchSegment>*) c1->getObject(seg));
		VectorNodeObject<BranchSegment>* tsegs2 = ((VectorNodeObject<BranchSegment>*) c2->getObject(seg));
		VectorNodeObject<double> v1  =tsegs1->at(0).alphas;
		VectorNodeObject<double> v2 = tsegs2->at(0).alphas;
		//cl1 = clock();
		VectorNodeObject<double> LHOODS (rootratemodel->getDists()->size(),0);
		for (unsigned int i = 0; i < dists->size(); i++) {
			if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
				VectorNodeObject<vector<int> >* exdist =
				((VectorNodeObject<vector<int> >*) node.getObject(en));
				int cou = count(exdist->begin(), exdist->end(), dists->at(i));
				if (cou == 0) {
					iter_ancsplits_just_int(rootratemodel, dists->at(i),
											leftdists, rightdists, weight);
					for (unsigned int j=0;j<leftdists.size();j++){
						int ind1 = leftdists[j];
						int ind2 = rightdists[j];
						LHOODS[i] += (v1.at(ind1)*v2.at(ind2)*weight);
					}
					LHOODS[i] *= Bs->at(i);
				}
			}
		}
		return LHOODS;
	}
	
}


