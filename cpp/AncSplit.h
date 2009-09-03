/*
 * AncSplit.h
 *
 *  Created on: Aug 15, 2009
 *      Author: smitty
 */

#ifndef ANCSPLIT_H_
#define ANCSPLIT_H_

#include <vector>
//using namespace std;

#include "RateModel.h"

class AncSplit{
	private:
		RateModel * model;
		double weight;
		double likelihood;
	public:
		AncSplit(RateModel * mod,int,int,int,double);
		RateModel * getModel();
		double getWeight();
		double getLikelihood();
		void setLikelihood(double li);
		int ancdistint;
		int ldescdistint;
		int rdescdistint;
};


#endif /* ANCSPLIT_H_ */
