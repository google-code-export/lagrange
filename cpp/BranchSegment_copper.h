/*
 * BranchSegment.h
 *
 *  Created on: Aug 16, 2009
 *      Author: smitty
 */

#ifndef BRANCHSEGMENT_COPPER_H_
#define BRANCHSEGMENT_COPPER_H_

#include <vector>
using namespace std;

#include "RateModel.h"

#include "vector_node_object.h"

class BranchSegment{
	private:
		double duration;
		int period;
		RateModel * model;
		vector<int> fossilareaindices;
		int startdistint;

	public:
		BranchSegment(double dur,int per);
		void setModel(RateModel * mod);
		//void setStartDist(vector<int> sd);
		void clearStartDist();
		double getDuration();
		int getPeriod();
		//vector<int> getStartDist();
		void set_start_dist_int(int d);
		int get_start_dist_int();
		RateModel * getModel();
		vector<int> getFossilAreas();
		void setFossilArea(int area);
		VectorNodeObject<double> * distconds;
		VectorNodeObject<double> alphas;
		VectorNodeObject<double> * ancdistconds;//for ancestral state reconstructions
};

#endif /* BRANCHSEGMENT_H_ */
