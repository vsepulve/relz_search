#ifndef _FACTORS_INDEX_H
#define _FACTORS_INDEX_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <algorithm>
#include <vector>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>
#include <sdsl/inv_perm_support.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "FactorsIterator.h"
#include "FactorsIteratorComparator.h"

using namespace sdsl;
using namespace std;

class FactorsIndex {

private: 
	unsigned int len_text;
	const char *ref;
	unsigned int len_ref;
	unsigned int n_factors;
	
	bit_vector arr_s;
	csa_wt<> fm_index;
	rmq_succinct_sct<false, bp_support_sada<256,32,rank_support_v5<> > > rmq;
	
	rrr_vector<127> rrr_s;
	
	rrr_vector<127>::select_1_type select1_s;
	rrr_vector<127>::select_0_type select0_s;
	rrr_vector<127>::select_1_type select1_b;
	rrr_vector<127>::select_0_type select0_b;
	
	inv_perm_support<> perm_inv;
	inv_perm_support<> perm;
	
	bit_vector arr_b;
	rrr_vector<127> rrr_b;
	
	inv_perm_support<> perm_x;
	inv_perm_support<> perm_y;
	inv_perm_support<> perm_y_inv;
	
	wt_int<rrr_vector<63>> wt;
	
	// Solo para pruebas
	int_vector<> ez;
	int_vector<> pi;
	int_vector<> pi_inv;
	int_vector<> pre_x_inv;
	int_vector<> pre_y;
	int_vector<> pre_y_inv;
	int_vector<> values_wt;
	
	// Cache de iteradores
	unordered_map<unsigned int, FactorsIterator> mapa_iterators;
	unordered_map<unsigned int, FactorsIteratorReverse> mapa_iterators_rev;
	
	void recursive_rmq(unsigned int ini, unsigned int fin, unsigned int crit);
	
	char getChar(unsigned int factor, unsigned int pos);
	
	char getCharRev(unsigned int factor, unsigned int pos);
			
	pair<unsigned int, unsigned int> getRangeY(const char *pattern);
			
	pair<unsigned int, unsigned int> getRangeX(const char *pattern);
			
	
public: 
	FactorsIndex();
	FactorsIndex(vector<pair<unsigned int, unsigned int> > &factors, unsigned int _len_text, const char *_ref, unsigned int _len_ref);
	~FactorsIndex();
	
	void find(const string &pattern);
	
};

#endif //_FACTORS_INDEX_H
