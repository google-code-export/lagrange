/*
 * BranchSegment.cpp
 *
 *  Created on: Aug 16, 2009
 *      Author: smitty
 */


#include "BranchSegment.h"
#include "RateModel.h"

#include <vector>
using namespace std;

BranchSegment::BranchSegment(double dur,int per){
	duration = dur;
	period = per;
}

void BranchSegment::setModel(RateModel * mod){
	model = mod;
}

void BranchSegment::setStartDist(vector<int> sd){
	startdist = sd;
}

void BranchSegment::clearStartDist(){
	startdist.clear();
}

double BranchSegment::getDuration(){
	return duration;
}

int BranchSegment::getPeriod(){
	return period;
}

vector<int> BranchSegment::getStartDist(){
	return startdist;
}

RateModel * BranchSegment::getModel(){
	return model;
}

vector<int> BranchSegment::getFossilAreas(){
	return fossilareaindices;
}

void BranchSegment::setFossilArea(int area){
	fossilareaindices.push_back(area);
}
