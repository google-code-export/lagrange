/*
 * PhyloTree.h
 *
 *  Created on: Aug 15, 2009
 *      Author: smitty
 */

#ifndef PHYLOTREE_H_
#define PHYLOTREE_H_

#include <string>
#include <map>
#include <vector>
#include "AncSplit.h"
using namespace std;

/*
#include <Phyl/Node.h>
#include <Phyl/TreeTemplate.h>
using namespace bpp;
*/

#include "tree.h"
#include "node.h"

class BioGeoTreeTools {
public :
	//TreeTemplate<Node> * getTreeFromString(string treestring) throw (Exception);
	//vector<int> getAncestors(TreeTemplate<bpp::Node> & tree, int nodeId);
	//int getLastCommonAncestor(TreeTemplate<bpp::Node> & tree, const vector<int>& nodeIds);
	Tree * getTreeFromString(string treestring);
	vector<int> getAncestors(Tree & tree, int nodeId);
	int getLastCommonAncestor(Tree & tree, const vector<int>& nodeIds);

	void summarizeSplits(Node * node,map<vector<int>,vector<AncSplit> > & ans,map<int,string> &areanamemaprev, RateModel * rm);
	void summarizeAncState(Node * node,vector<double> & ans,map<int,string> &areanamemaprev, RateModel * rm);
};

#endif /* PHYLOTREE_H_ */
