/*
 * BranchSegment.h
 *
 *  Created on: Aug 16, 2009
 *      Author: smitty
 */

#ifndef BRANCHSEGMENT_H_
#define BRANCHSEGMENT_H_

#include <vector>
using namespace std;

#include "RateModel.h"

#include <Utils/BppVector.h>
using namespace bpp;

class BranchSegment{
	private:
		double duration;
		int period;
		RateModel * model;
		vector<int> startdist;
		vector<int> fossilareaindices;

	public:
		BranchSegment(double dur,int per);
		void setModel(RateModel * mod);
		void setStartDist(vector<int> sd);
		void clearStartDist();
		double getDuration();
		int getPeriod();
		vector<int> getStartDist();
		RateModel * getModel();
		vector<int> getFossilAreas();
		void setFossilArea(int area);
		Vector<double> * distconds;
		Vector<double> * ancdistconds;//for ancestral state reconstructions
};

#endif /* BRANCHSEGMENT_H_ */
