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

double AncSplit::getWeight(){
	return weight;
}

double AncSplit::getLikelihood(){
	return likelihood;
}

void AncSplit::setLikelihood(double li){
	likelihood = li;
}
