#ifndef _FACTORS_INDEX_H
#define _FACTORS_INDEX_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>

#include <algorithm>
#include <vector>
#include <unordered_map>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>
#include <sdsl/inv_perm_support.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "FactorsConfig.h"
#include "FactorsIterator.h"
#include "FactorsFastIteratorComparator.h"
#include "NanoTimer.h"

using namespace sdsl;
using namespace std;

class FactorsIndex {

private: 
	
	unsigned int len_text;
	char *ref_text;
	unsigned int len_ref;
	unsigned int n_factors;
	bool omit_text;
	
	fm_index_type fm_index;
	
	rmq_type rmq;
	
	bits_s_type bits_s;
	bits_s_type::select_1_type select1_s;
	bits_s_type::select_0_type select0_s;
	
	bits_b_type bits_b;
	bits_b_type::select_1_type select1_b;
	bits_b_type::select_0_type select0_b;
	
	int_vector<> pi;
	int_vector<> pi_inv;
	
	int_vector<> arr_x;
	int_vector<> arr_y;
	
	wt_type wt;
	
	// Prueba de aceleracion de recursive_rmq almacenando los datos de los factores descomprimidos
	static const bool precompute_rmq = false;
	vector<unsigned int> arr_tu;
	vector<unsigned int> arr_pu;
	vector<unsigned int> arr_lu;
	
	void recursive_rmq(unsigned int ini, unsigned int fin, unsigned int min_pos, unsigned int occ_ref, vector<unsigned int> &results);
			
	pair<unsigned int, unsigned int> getRangeY(const char *pattern);
	pair<unsigned int, unsigned int> getRangeX(const char *pattern);
	
	template <typename ItereatorType>
	bool factorLess(unsigned int factor, const char *pattern, unsigned int len, bool equal = false);
	
	// Testing methods, just to debbug
	char getChar(unsigned int factor, unsigned int pos, unsigned int max_len = 100);
	char getCharRev(unsigned int factor, unsigned int pos, unsigned int max_len = 100);
	// Iterators for getChar, just to debbug
	unordered_map<unsigned int, FactorsIterator> mapa_iterators;
	unordered_map<unsigned int, FactorsIteratorReverse> mapa_iterators_rev;
	
public: 
	FactorsIndex();
	FactorsIndex(vector<pair<unsigned int, unsigned int> > &factors, char *full_text, unsigned int _len_text, const char *_ref_text, unsigned int _len_ref, bool _omit_text = false);
	~FactorsIndex();
	
	void findTimes(const string &pattern, vector<unsigned int> &results);
	
	void printSize();
	
	void save(const string &file_base);
	
	void load(const string &file_base);
	
	unsigned long long querytime_p1;
	unsigned long long querytime_p2;
	unsigned long long querytime_p3;
	unsigned long long querytime_p4;
	
	unsigned int occs_a;
	unsigned int occs_b;
	unsigned int occs_c;
	
};








#endif //_FACTORS_INDEX_H
