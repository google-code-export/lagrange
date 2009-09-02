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
		/*
		 * make these point to the dist number in the rootratemodel dists
		 */
		vector<int> ancdist;
		vector<int> ldescdist;
		vector<int> rdescdist;

		double weight;
		double likelihood;
	public:
		/*AncSplit(RateModel * mod, vector<int> dist, vector<int > ldesc, vector<int> rdesc);
		AncSplit(RateModel * mod, vector<int> dist, vector<int > ldesc, vector<int> rdesc, double we);
		AncSplit(RateModel * mod, vector<int> dist, vector<int > ldesc, vector<int> rdesc, double we, double like);*/
		AncSplit(RateModel * mod,int,int,int,double);
		RateModel * getModel();
/*		vector<int> getAncDist();
		vector<int> getLDescDist();
		vector<int> getRDescDist();*/
		double getWeight();
		double getLikelihood();
		void setLikelihood(double li);
		int ancdistint;
		int ldescdistint;
		int rdescdistint;
};


#endif /* ANCSPLIT_H_ */
