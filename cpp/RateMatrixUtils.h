/*
 * RateMatrixUtils.h
 *
 *  Created on: Aug 13, 2009
 *      Author: Stephen A. Smith
 */

#ifndef RATEMATRIXUTILS_H_
#define RATEMATRIXUTILS_H_
#include "RateModel.h"
#include "AncSplit.h"
#include <vector>
using namespace std;
/*
 * this constructs a P matrix using a Q matrix
 */
vector< vector<double> > QMatrixToPmatrix(vector< vector<double> > & Q, double t);
void calcMatExp(int * ia,int * ja, double * a, int n);
double calculate_vector_double_sum(vector<double> & in);
int calculate_vector_int_sum(vector<int> * in);
int calculate_vector_int_sum_xor(vector<int> & in,vector<int> & in2);
int locate_vector_int_single_xor(vector<int> & in, vector<int> & in2);
//bool is_vector_int_equal_to_vector_int(vector<int> &,vector<int> &);
vector<AncSplit> iter_ancsplits(RateModel *rm, vector<int> & dist);
void iter_ancsplits_just_int(RateModel *rm, vector<int> & dist,vector<int> & leftdists, vector<int> & rightdists, double & weight);
void print_vector_int(vector<int> & in);
void print_vector_double(vector<double> & in);
int get_vector_int_index_from_multi_vector_int(vector<int> * in, vector<vector<int> > * in2);

vector<vector<int> > generate_dists_from_num_max_areas(int totalnum,int numareas);

int get_size_for_coo(vector<vector<double> > & inm, double t);
void convert_matrix_to_coo_for_fortran(vector<vector<double> > & inmatrix, double t, int * ia, int * ja, double * a);
void convert_matrix_to_coo_for_fortran_vector(vector<vector<double> > & inmatrix, vector<int> & ia, vector<int> & ja, vector<double> & a);
void convert_matrix_to_single_row_for_fortran(vector<vector<double> > & inmatrix, double t, double * H);

vector<vector<vector<double> > > processRateMatrixConfigFile(string filename, int numareas, int nperiods);

vector<int> get_columns_for_sparse(vector<double> &,RateModel *);

/*
	this is for pthread sparse columns
 */
struct sparse_thread_data{
	int thread_id;
	vector<int> columns;
	vector<vector<double> > presults;
	RateModel * rm;
	double t;
	int period;
};

void * sparse_column_pmatrix_pthread_go(void *threadarg);

#endif /* RATEMATRIXUTILS_H_ */
