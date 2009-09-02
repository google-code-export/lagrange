/*
 * AncSplit.cpp
 *
 *  Created on: Aug 15, 2009
 *      Author: Stephen A. Smith
 */

#include "AncSplit.h"
#include "RateModel.h"

#include <vector>
using namespace std;

/*AncSplit::AncSplit(RateModel * mod, vector<int> dist, vector<int > ldesc, vector<int> rdesc){
	model = mod;
	ancdist = dist;
	ldescdist = ldesc;
	rdescdist = rdesc;
	likelihood = 0;
	weight = 0;
}

AncSplit::AncSplit(RateModel * mod, vector<int> dist, vector<int > ldesc, vector<int> rdesc, double we){
	model = mod;
	ancdist = dist;
	ldescdist = ldesc;
	rdescdist = rdesc;
	weight = we;
	likelihood = 0;
}

AncSplit::AncSplit(RateModel * mod, vector<int> dist, vector<int > ldesc, vector<int> rdesc, double we, double like){
	model = mod;
	ancdist = dist;
	ldescdist = ldesc;
	rdescdist = rdesc;
	weight = we;
	likelihood = like;
}*/

AncSplit::AncSplit(RateModel * mod,int dist,int ldesc,int rdesc,double we){
	model = mod;
	ancdistint = dist;
	ldescdistint = ldesc;
	rdescdistint = rdesc;
	weight = we;
}

RateModel * AncSplit::getModel(){
	return model;
}

/*vector<int> AncSplit::getAncDist(){
	return ancdist;
}

vector<int> AncSplit::getLDescDist(){
	return ldescdist;
}

vector<int> AncSplit::getRDescDist(){
	return rdescdist;
}*/

double AncSplit::getWeight(){
	return weight;
}

double AncSplit::getLikelihood(){
	return likelihood;
}

void AncSplit::setLikelihood(double li){
	likelihood = li;
}
