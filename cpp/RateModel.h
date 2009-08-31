/*
 * RateMatrix.h
 *
 *  Created on: Aug 14, 2009
 *      Author: smitty
 */

#ifndef RATEMODEL_H_
#define RATEMODEL_H_


//#include "AncSplit.h"

#include <vector>
#include <map>
#include <string>
using namespace std;

class RateModel{
	private:
		int nareas;
		bool globalext;
		vector<string> labels;
		vector<double> periods;
		vector<vector<int> > dists;
		map<vector<int>,vector<vector<vector<int> > > > iter_dists;
		map<vector<int>,string> distsmap;
		vector< vector< vector<double> > >D;
		vector< vector< vector<double> > >Dmask;
		vector< vector<double> > E;
		vector< vector< vector<double> > >Q;
		vector< vector< vector<double> > >QT;//transposed for sparse
		vector< vector< vector<double> > >P;
		vector<int> nzs;
		vector<vector<int> > ia_s;
		vector<vector<int> > ja_s;
		vector<vector<double> > a_s;

	public:
		RateModel(int na, bool ge, vector<double> pers,bool);
		void setup_dists();
		void setup_dists(vector<vector<int> >, bool);
		void setup_Dmask();
		void set_Dmask_cell(int period, int area, int area2, double prob, bool sym);
		void setup_D(double d);
		void setup_E(double e);
		void set_Qdiag(int period);
		void setup_Q();
		vector<vector<double > > setup_P(int period, double t);
		vector<vector<double > > setup_fortran_P(int period, double t);
		vector<vector<double > > setup_sparse_full_P(int period, double t);
		vector<double > setup_sparse_single_column_P(int period, double t, int column);
		string Q_repr(int period);
		string P_repr(int period);
		vector<vector<int> > enumerate_dists();
		vector<vector<vector<int> > > iter_dist_splits(vector<int> & dist);
		//vector<AncSplit> iter_ancsplits(vector<int> dist);
		vector<vector<int> > * getDists();
		vector<vector<vector<int> > > * get_iter_dist_splits(vector<int> & dist);
		void remove_dist(vector<int> dist);
		void iter_all_dist_splits();
		bool sparse;
};

#endif /* RATEMATRIX_H_ */
