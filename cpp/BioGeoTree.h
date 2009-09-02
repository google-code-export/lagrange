/*
 * BioGeoTree.h
 *
 *  Created on: Aug 15, 2009
 *      Author: smitty
 */

#ifndef BIOGEOTREE_H_
#define BIOGEOTREE_H_

#include <vector>
#include <string>
#include <map>
using namespace std;

#include "RateModel.h"
#include "AncSplit.h"

#include <Phyl/TreeTemplateTools.h>
#include <Phyl/TreeTemplate.h>
#include <Phyl/Newick.h>
#include <Phyl/Tree.h>
#include <Utils/BppVector.h>

class BioGeoTree{
private:
	bpp::TreeTemplate<bpp::Node> * tree;
	vector<double> periods;
	string seg;
	string age;
	string dc;
	string en;
	string nasp;
	string ast;
	string andc;
	string tvec;
	clock_t cl1;
	clock_t cl2;
	clock_t c3;
	clock_t c4;
	clock_t c5;
	clock_t c6;
	int curancstatenodeid;
	vector<int> * columns;
	vector<int> * whichcolumns;
	RateModel * rootratemodel;
	bool store_p_matrices;
	bool use_stored_matrices;

public:
	BioGeoTree(bpp::TreeTemplate<bpp::Node> * tr, vector<double> ps);
	void cleanNodesAndSegs();
	void set_default_model(RateModel * mod);
	void update_default_model(RateModel * mod);
	double eval_likelihood(bool marg);
	void set_excluded_dist(vector<int> ind,bpp::Node * node);
	void set_tip_conditionals(map<string,vector<int> > distrib_data);
	bpp::Vector<double> conditionals(bpp::Node & node, bool marg, bool , bool, bool);
	//, bpp::Vector<double>&);
	void ancdist_conditional_lh(bpp::Node & node, bool marg);
	double eval_likelihood_ancstate(bool marginal,bpp::Node & startnode);
	void ancstate_ancdist_conditional_lh(bpp::Node * fromnode,bpp::Node * node, bool marginal);
	vector<AncSplit> ancstate_calculation(bpp::Node & node,vector<int> & dist,bool marg);
	map<vector<int>,vector<AncSplit> > ancstate_calculation_all_dists(bpp::Node & node, bool marginal);
	void setFossilatNodeByMRCA(vector<string> nodeNames, int fossilarea);
	void setFossilatNodeByMRCA_id(int id, int fossilarea);
	void setFossilatBranchByMRCA(vector<string> nodeNames, int fossilarea, double age);
	void setFossilatBranchByMRCA_id(int id, int fossilarea, double age);
	void set_store_p_matrices(bool);
	void set_use_stored_matrices(bool);

	double ti;
	double ti2;
	double ti3;
};

#endif /* BIOGEOTREE_H_ */
