#ifndef _FACTORS_INDEX_V3_H
#define _FACTORS_INDEX_V3_H

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
#include "FactorsFastIteratorComparator.h"
#include "NanoTimer.h"

#include "BitsUtils.h"
#include "KarpRabin.h"

using namespace sdsl;
using namespace std;

class FactorsIndexV3 {

private: 
	unsigned int len_text;
	const char *ref_text;
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
	
	KarpRabin *karp_rabin;
	vector<unsigned long long> arr_kr_ref;
	vector<unsigned long long> arr_kr_s;
	
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
	
	// Falta la estructura para agregar efectivamente los resultados, quiza un vector de posiciones
	void recursive_rmq(unsigned int ini, unsigned int fin, unsigned int min_pos, unsigned int occ_ref, vector<unsigned int> &results);
	
	char getChar(unsigned int factor, unsigned int pos);
	
	char getCharRev(unsigned int factor, unsigned int pos);
			
	pair<unsigned int, unsigned int> getRangeY(const char *pattern);
			
	pair<unsigned int, unsigned int> getRangeX(const char *pattern);
			
	
public: 
	FactorsIndexV3();
	FactorsIndexV3(vector<pair<unsigned int, unsigned int> > &factors, char *full_text, unsigned int _len_text, const char *_ref_text, unsigned int _len_ref, KarpRabin *_karp_rabin, const char *kr_frases_file, bool load_kr_frases = false);
	~FactorsIndexV3();
	
	void find(const string &pattern, vector<unsigned int> &results);
	
};

#endif //_FACTORS_INDEX_H