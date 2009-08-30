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
#include <Phyl/Node.h>
#include <Phyl/TreeTemplate.h>
#include "AncSplit.h"
using namespace bpp;
using namespace std;

class BioGeoTreeTools {
	public :
		TreeTemplate<Node> * getTreeFromString(string treestring) throw (Exception);
		vector<int> getAncestors(TreeTemplate<bpp::Node> & tree, int nodeId);
		int getLastCommonAncestor(TreeTemplate<bpp::Node> & tree, const vector<int>& nodeIds);
		void summarizeSplits(Node *,map<vector<int>,vector<AncSplit> >,map<int,string>);
		void summarizeAncState(Node,map<vector<int>,vector<AncSplit> >);
};

#endif /* PHYLOTREE_H_ */
