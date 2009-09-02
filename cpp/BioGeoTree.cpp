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

using namespace std;

#include "BioGeoTree.h"
#include "BioGeoTreeTools.h"
#include "BranchSegment.h"
#include "RateMatrixUtils.h"
#include "RateModel.h"
#include "AncSplit.h"

#include <Phyl/TreeTemplate.h>
#include <Phyl/TreeTemplateTools.h>
#include <Phyl/Newick.h>
#include <Phyl/Tree.h>
#include <Utils/BppVector.h>
using namespace bpp;

namespace {
	inline double MAX(const double &a, const double &b)
	        {return b > a ? (b) : double(a);}
}

BioGeoTree::BioGeoTree(TreeTemplate<Node> * tr, vector<double> ps){
	seg = "segments";
	age = "age";
	dc = "dist_conditionals";
	en = "excluded_dists";
	nasp = "ancsplit";
	ast = "anstate";//current node with ancstate
	andc = "anc_dist_conditionals";
	tvec = "temporary vector";
	store_p_matrices = false;
	use_stored_matrices = false;
	tree = tr;
	periods = ps;
	TreeTemplateTools * ttt = new TreeTemplateTools();
	/*
	 * initialize each node with segments
	 */
	for (unsigned int i=0;i<tree->getNumberOfNodes();i++){
		Vector<BranchSegment> * segs = new Vector<BranchSegment>();
		tree->setNodeProperty(i,seg,*segs);
		Vector<vector<int> > * ens = new Vector<vector<int> >;
		tree->setNodeProperty(i,en,*ens);
		Vector<AncSplit> * ancsplits = new Vector<AncSplit>();
		tree->setNodeProperty(i,nasp,*ancsplits);
	}
	/*
	 * initialize the actual branch segments for each node
	 */
	for (unsigned int i=0;i<tree->getNumberOfNodes();i++){
		if (tree->getNode(i)->hasFather()){
			vector<double> pers(periods);
			double anc = ttt->getHeight(*tree->getNode(i)->getFather());
			double des = ttt->getHeight(*tree->getNode(i));
			//assert anc > des
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
							((bpp::Vector<BranchSegment>*) tree->getNodeProperty(i,seg))->push_back(tseg);
						}
						t += pers[j];
					}
					if (t > anc){
						break;
					}
				}
			}else{
				BranchSegment tseg = BranchSegment(tree->getNode(i)->getDistanceToFather(),0);
				((bpp::Vector<BranchSegment>*) tree->getNodeProperty(i,seg))->push_back(tseg);
			}
		}
	}
	delete ttt;
}

void BioGeoTree::set_store_p_matrices(bool i){
	store_p_matrices = i;
}

void BioGeoTree::set_use_stored_matrices(bool i){
	use_stored_matrices = i;
}

void BioGeoTree::cleanNodesAndSegs(){
	for(unsigned int i=0;i<tree->getNumberOfNodes();i++){
		((bpp::Vector<AncSplit>*) tree->getNodeProperty(i,nasp))->clear();
	}
}

void BioGeoTree::set_default_model(RateModel * mod){
	rootratemodel = mod;
	for(unsigned int i=0;i<tree->getNumberOfNodes();i++){
		bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) tree->getNode(i)->getNodeProperty(seg));
		for(unsigned int j=0;j<tsegs->size();j++){
			tsegs->at(j).setModel(mod);
			Vector<double> * distconds = new Vector<double> (rootratemodel->getDists()->size(), 0);
			tsegs->at(j).distconds = distconds;
			Vector<double> * ancdistconds = new Vector<double> (rootratemodel->getDists()->size(), 0);
			tsegs->at(j).ancdistconds = ancdistconds;
        }
        //Vector<double> * tempvec = new Vector<double>(rootratemodel->getDists()->size(), 0);
        //tree->setNodeProperty(i,tvec,*tempvec);
	}
	Vector<double> * distconds = new Vector<double> (rootratemodel->getDists()->size(), 0);
	tree->getRootNode()->setNodeProperty(dc,*distconds);
	Vector<double> * ancdistconds = new Vector<double> (rootratemodel->getDists()->size(), 0);
	tree->getRootNode()->setNodeProperty(andc,*ancdistconds);
}

void BioGeoTree::update_default_model(RateModel * mod){
	rootratemodel = mod;
	for(unsigned int i=0;i<tree->getNumberOfNodes();i++){
		bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) tree->getNode(i)->getNodeProperty(seg));
		for(unsigned int j=0;j<tsegs->size();j++){
			tsegs->at(j).setModel(mod);
		}
	}
}

void BioGeoTree::set_tip_conditionals(map<string,vector<int> > distrib_data){
	for(unsigned int i=0;i<tree->getNumberOfLeaves();i++){
		bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) tree->getLeaves().at(i)->getNodeProperty(seg));
		RateModel * mod = tsegs->at(0).getModel();
		int ind1 = get_vector_int_index_from_multi_vector_int(
				&distrib_data[tree->getLeaves().at(i)->getName()],mod->getDists());
		tsegs->at(0).distconds->at(ind1) = 1.0;
	}
}

void BioGeoTree::set_excluded_dist(vector<int> ind,Node * node){
	((bpp::Vector<vector<int> >*) node->getNodeProperty(en))->push_back(ind);
}

/*
 * **************************************************
 *
 *
 *
 * **************************************************
 */

double BioGeoTree::eval_likelihood(bool marginal){
	cleanNodesAndSegs();
	if( rootratemodel->sparse == true){
		columns = new vector<int>(rootratemodel->getDists()->size());
		whichcolumns = new vector<int>();
	}
	ancdist_conditional_lh(*tree->getRootNode(),marginal);
	if( rootratemodel->sparse == true){
		delete columns;
		delete whichcolumns;
	}
	return calculate_vector_double_sum(*
			(bpp::Vector<double>*) tree->getRootNode()->getNodeProperty(dc));

}

Vector<double> BioGeoTree::conditionals(Node & node, bool marginal,
						bool curancstate, bool calcancstate, bool sparse){
	Vector<double> distconds;
	bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) node.getNodeProperty(seg));
	if(calcancstate == false){
		distconds = *tsegs->at(0).distconds;
	}else{
		if (curancstate == false){
			distconds = *tsegs->at(0).distconds;
		}else{//calc == true and cur == true
			distconds = *tsegs->at(0).ancdistconds;
		}
	}
	for(unsigned int i=0;i<tsegs->size();i++){
		if(calcancstate == false){
			for(unsigned int j=0;j<distconds.size();j++){
				tsegs->at(i).distconds->at(j) = distconds.at(j);
			}
		}else{
			for(unsigned int j=0;j<distconds.size();j++){
				tsegs->at(i).ancdistconds->at(j) = distconds.at(j);
			}
		}
		RateModel * rm = tsegs->at(i).getModel();
		Vector<double> * v = new Vector<double> (rootratemodel->getDists()->size(), 0);
		vector<int> distrange;
		if(tsegs->at(i).get_start_dist_int() != -666){
		//if (tsegs->at(i).getStartDist().size()>0){//!=NULL
			int ind1 = tsegs->at(i).get_start_dist_int();
			//(*rootratemodel->get_dists_int_map())[tsegs->at(i).getStartDist()];
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
	return distconds;
}

void BioGeoTree::ancdist_conditional_lh(Node & node, bool marginal){
	Vector<double> distconds(rootratemodel->getDists()->size(), 0);
	if (node.isLeaf()==false){//is not a tip
		Node * c1 = node.getSon(0);
		Node * c2 = node.getSon(1);
		RateModel * model;
		if(node.hasFather()==true){
			bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) node.getNodeProperty(seg));
			model = tsegs->at(0).getModel();
		}else{
			model = rootratemodel;
		}
		ancdist_conditional_lh(*c1,marginal);
		ancdist_conditional_lh(*c2,marginal);
		bool sparse = rootratemodel->sparse;
		Vector<double> v1;
		Vector<double> v2;
		if(sparse == true){
			//getcolumns
			bpp::Vector<BranchSegment>* c1tsegs = ((bpp::Vector<BranchSegment>*) c1->getNodeProperty(seg));
			bpp::Vector<BranchSegment>* c2tsegs = ((bpp::Vector<BranchSegment>*) c2->getNodeProperty(seg));
			vector<int> lcols = get_columns_for_sparse(*c1tsegs->at(0).distconds,rootratemodel);
			vector<int> rcols = get_columns_for_sparse(*c2tsegs->at(0).distconds,rootratemodel);
			for(unsigned int i=0;i<lcols.size();i++){
				if(lcols[i]==1 || rcols[i] ==1){
					columns->at(i)=1;
					if(i!=0)
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
		v1 =conditionals(*c1,marginal,false,false,sparse);
		v2 =conditionals(*c2,marginal,false,false,sparse);

		vector<vector<int> > * dists = rootratemodel->getDists();
//		map<vector<int>,int> * distmap = rootratemodel->get_dists_int_map(); 
		Vector<AncSplit> * ancsplits = (Vector<AncSplit> *) node.getNodeProperty(nasp);
		//cl1 = clock();
		for (unsigned int i=0;i<dists->size();i++){
			//if (calculate_vector_int_sum(&dists->at(i)) > 0){
			if(accumulate(dists->at(i).begin(),dists->at(i).end(),0) > 0){
				double lh = 0.0;
				bpp::Vector<vector<int> >* exdist = ((bpp::Vector<vector<int> >*) node.getNodeProperty(en));
				int cou = count(exdist->begin(),exdist->end(),dists->at(i));
				if(cou == 0){
					vector<AncSplit> ans = iter_ancsplits(rootratemodel,dists->at(i));
					for (unsigned int j=0;j<ans.size();j++){
						int ind1 = ans[j].ldescdistint;
						//(*distmap)[ans[j].getLDescDist()];
						int ind2 = ans[j].rdescdistint;
						//(*distmap)[ans[j].getRDescDist()];
						double lh_part = v1.at(ind1)*v2.at(ind2);
						lh += (lh_part * ans[j].getWeight());
						ans[j].setLikelihood(lh_part);
						ancsplits->push_back(ans[j]);
					}
				}
				distconds.at(i)= lh;
			}
		}
		///cl2 = clock();
		//ti += cl2-cl1;
	}else{
		bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) node.getNodeProperty(seg));
		distconds = *tsegs->at(0).distconds;
	}
	if(node.hasFather() == true){
		bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) node.getNodeProperty(seg));
		for(unsigned int i=0;i<distconds.size();i++){
			tsegs->at(0).distconds->at(i) = distconds.at(i);
		}
	}else{
		for(unsigned int i=0;i<distconds.size();i++){
			((Vector<double>*)node.getNodeProperty(dc))->at(i) = distconds.at(i);
		}
	}

}

/*
 * all ancestral state calculations
 */
double BioGeoTree::eval_likelihood_ancstate(bool marginal,bpp::Node &startnode){
	cleanNodesAndSegs();
	if( rootratemodel->sparse == true){
		columns = new vector<int>(rootratemodel->getDists()->size());
	}
	ancstate_ancdist_conditional_lh(NULL,&startnode,marginal);
	if( rootratemodel->sparse == true){
		delete columns;
	}
	return calculate_vector_double_sum(*
			(bpp::Vector<double>*) tree->getRootNode()->getNodeProperty(dc));
}

/*
    This calculates the conditionals for internal nodes when calculating ancestral states.
    The major difference between this and typical conditional calculation is that it calculates
    down the backbone and therefore requires the node calling this procedure.

 */
void BioGeoTree::ancstate_ancdist_conditional_lh(Node * fromnode, Node * node, bool marginal){
	Vector<double> distconds(rootratemodel->getDists()->size(), 0);
	if (node->isLeaf()==false){//is not a tip
		Node * c1 = node->getSon(0);
		Node * c2 = node->getSon(1);
		RateModel * model;
		if(node->hasFather()==true){
			bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) node->getNodeProperty(seg));
			model = tsegs->at(0).getModel();
		}else{
			model = rootratemodel;
		}
		bool sparse = rootratemodel->sparse;
		Vector<double> v1;
		Vector<double> v2;
		if(sparse == true){
			bpp::Vector<BranchSegment>* c1tsegs = ((bpp::Vector<BranchSegment>*) c1->getNodeProperty(seg));
			bpp::Vector<BranchSegment>* c2tsegs = ((bpp::Vector<BranchSegment>*) c2->getNodeProperty(seg));
			vector<int> lcols = get_columns_for_sparse(*c1tsegs->at(0).distconds,rootratemodel);
			vector<int> rcols = get_columns_for_sparse(*c2tsegs->at(0).distconds,rootratemodel);
			for(unsigned int i=0;i<lcols.size();i++){
				if(lcols[i]==1 || rcols[i] ==1)
					columns->at(i)=1;
				else
					columns->at(i)=0;
			}
			if(calculate_vector_int_sum(columns)==0){
				for(unsigned int i=0;i<lcols.size();i++){
					columns->at(i)=1;
				}
			}
			columns->at(0) = 0;
		}
		if(node->getId() == curancstatenodeid){
			v1 =conditionals(*c1,marginal,false,true,sparse);
			v2 =conditionals(*c2,marginal,false,true,sparse);
		}else{
			if(c1 == fromnode){
				v1 =conditionals(*c1,marginal,true,true,sparse);
				v2 =conditionals(*c2,marginal,false,true,sparse);
			}else if(c2 == fromnode){
				v1 = conditionals(*c1,marginal,false,true,sparse);
				v2 = conditionals(*c2,marginal,true,true,sparse);
			}
		}

		vector<vector<int> > * dists = rootratemodel->getDists();
//		map<vector<int>,int> * distmap = rootratemodel->get_dists_int_map(); 
		Vector<AncSplit> * ancsplits = (Vector<AncSplit> *) node->getNodeProperty(nasp);
		//cl1 = clock();
		for (unsigned int i=0;i<dists->size();i++){
			//if (calculate_vector_int_sum(&dists->at(i)) > 0){
			if(accumulate(dists->at(i).begin(),dists->at(i).end(),0) > 0){
				double lh = 0.0;
				bpp::Vector<vector<int> >* exdist = ((bpp::Vector<vector<int> >*) node->getNodeProperty(en));
				int cou = count(exdist->begin(),exdist->end(),dists->at(i));
				if(cou == 0){
					vector<AncSplit> ans = iter_ancsplits(rootratemodel,dists->at(i));
					for (unsigned int j=0;j<ans.size();j++){
						int ind1 = ans[j].ldescdistint;
						//(*distmap)[ans[j].getLDescDist()];
						int ind2 = ans[j].rdescdistint;
						//(*distmap)[ans[j].getRDescDist()];
						double lh_part = v1.at(ind1)*v2.at(ind2);
						lh += (lh_part * ans[j].getWeight());
						ans[j].setLikelihood(lh_part);
						ancsplits->push_back(ans[j]);
					}
				}
				distconds.at(i)= lh;
			}
		}
		//cl2 = clock();
		//ti += cl2-cl1;
		//cout << ti/CLOCKS_PER_SEC << endl;
	}else{
		bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) node->getNodeProperty(seg));
		distconds = *tsegs->at(0).distconds;
	}
	if(node->hasFather() == true){
		bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) node->getNodeProperty(seg));
		for(unsigned int i=0;i<distconds.size();i++){
			tsegs->at(0).ancdistconds->at(i) = distconds.at(i);
		}
	}else{
		for(unsigned int i=0;i<distconds.size();i++){
			((Vector<double>*)node->getNodeProperty(dc))->at(i) = distconds.at(i);
		}
	}
	if(node->hasFather()){
		ancstate_ancdist_conditional_lh(node,node->getFather(),marginal);
	}
}

vector<AncSplit> BioGeoTree::ancstate_calculation(bpp::Node & node,vector<int> & dist, bool marginal){
	vector<AncSplit> ans = iter_ancsplits(rootratemodel,dist);
	if (node.isLeaf()==false){//is not a tip
		//cout << node->getId()<< " " << tree->getRootId() << endl;
		Node * c1 = node.getSon(0);
		Node * c2 = node.getSon(1);
		bpp::Vector<BranchSegment>* tsegs1 = ((bpp::Vector<BranchSegment>*) c1->getNodeProperty(seg));
		bpp::Vector<BranchSegment>* tsegs2 = ((bpp::Vector<BranchSegment>*) c2->getNodeProperty(seg));
		for (unsigned int i=0;i<ans.size();i++){
			//tsegs1->at(tsegs1->size()-1).setStartDist(ans[i].getLDescDist());
			tsegs1->at(tsegs1->size()-1).set_start_dist_int(ans[i].ldescdistint);
			//tsegs2->at(tsegs2->size()-1).setStartDist(ans[i].getRDescDist());
			tsegs2->at(tsegs2->size()-1).set_start_dist_int(ans[i].rdescdistint);
			double lh = eval_likelihood_ancstate(marginal,node);
			ans[i].setLikelihood(lh);
		}
		tsegs1->at(tsegs1->size()-1).clearStartDist();
		tsegs2->at(tsegs2->size()-1).clearStartDist();
	}
	return ans;
}

map<vector<int>,vector<AncSplit> > BioGeoTree::ancstate_calculation_all_dists(bpp::Node & node, bool marginal){
	curancstatenodeid = node.getId();
	map<vector<int>,vector<AncSplit> > ret;
	for(unsigned int j=0;j<rootratemodel->getDists()->size();j++){
		vector<int> dist = rootratemodel->getDists()->at(j);
		vector<AncSplit> ans = iter_ancsplits(rootratemodel,dist);
		if (node.isLeaf()==false){//is not a tip
			//cout << node->getId()<< " " << tree->getRootId() << endl;
			Node * c1 = node.getSon(0);
			Node * c2 = node.getSon(1);
			bpp::Vector<BranchSegment>* tsegs1 = ((bpp::Vector<BranchSegment>*) c1->getNodeProperty(seg));
			bpp::Vector<BranchSegment>* tsegs2 = ((bpp::Vector<BranchSegment>*) c2->getNodeProperty(seg));
			for (unsigned int i=0;i<ans.size();i++){
				//tsegs1->at(tsegs1->size()-1).setStartDist(ans[i].getLDescDist());
				tsegs1->at(tsegs1->size()-1).set_start_dist_int(ans[i].ldescdistint);
				//tsegs2->at(tsegs2->size()-1).setStartDist(ans[i].getRDescDist());
				tsegs2->at(tsegs2->size()-1).set_start_dist_int(ans[i].rdescdistint);
				double lh = eval_likelihood_ancstate(marginal,node);
				//cout << lh << endl;
				ans[i].setLikelihood(lh);
			}
			tsegs1->at(tsegs1->size()-1).clearStartDist();
			tsegs2->at(tsegs2->size()-1).clearStartDist();
		}
		ret[dist] = ans;
	}
	curancstatenodeid = NULL;
	return ret;
}


/*
 * ********************************************
 *
 * adds fossils either at the node or along a branch
 *
 * ********************************************
 */
void BioGeoTree::setFossilatNodeByMRCA(vector<string> nodeNames, int fossilarea){
	BioGeoTreeTools tt;
	vector<int> nodeIds;
	for(unsigned int i=0;i<nodeNames.size();i++){
		nodeIds.push_back(tree->getNode(nodeNames[i])->getId());
	}
	int id = tt.getLastCommonAncestor(*tree,nodeIds);
	vector<vector<int> > * dists = rootratemodel->getDists();
	for(unsigned int i=0;i<dists->size();i++){
		if(dists->at(i).at(fossilarea) == 0){
			bpp::Vector<vector<int> > * exd = ((bpp::Vector<vector<int> > *) tree->getNodeProperty(id,en));
			exd->push_back(dists->at(i));
		}
	}
}
void BioGeoTree::setFossilatNodeByMRCA_id(int id, int fossilarea){
	vector<vector<int> > * dists = rootratemodel->getDists();
	for(unsigned int i=0;i<dists->size();i++){
		if(dists->at(i).at(fossilarea) == 0){
			bpp::Vector<vector<int> > * exd = ((bpp::Vector<vector<int> > *) tree->getNodeProperty(id,en));
			exd->push_back(dists->at(i));
		}
	}
}
void BioGeoTree::setFossilatBranchByMRCA(vector<string> nodeNames, int fossilarea, double age){
	BioGeoTreeTools tt;
	TreeTemplateTools * ttt = new TreeTemplateTools();
	vector<int> nodeIds;
	for(unsigned int i=0;i<nodeNames.size();i++){
		nodeIds.push_back(tree->getNode(nodeNames[i])->getId());
	}
	int id = tt.getLastCommonAncestor(*tree,nodeIds);
	bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) tree->getNode(id)->getNodeProperty(seg));
	double startage = ttt->getHeight(*tree->getNode(id));
	for(unsigned int i=0;i<tsegs->size();i++){
		if(age > startage && age < (startage+tsegs->at(i).getDuration())){
			tsegs->at(i).setFossilArea(fossilarea);
		}
		startage += tsegs->at(i).getDuration();
	}
	delete ttt;
}
void BioGeoTree::setFossilatBranchByMRCA_id(int id, int fossilarea, double age){
	TreeTemplateTools * ttt = new TreeTemplateTools();
	bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) tree->getNode(id)->getNodeProperty(seg));
	double startage = ttt->getHeight(*tree->getNode(id));
	for(unsigned int i=0;i<tsegs->size();i++){
		if(age > startage && age < (startage+tsegs->at(i).getDuration())){
			tsegs->at(i).setFossilArea(fossilarea);
		}
		startage += tsegs->at(i).getDuration();
	}
	delete ttt;
}

