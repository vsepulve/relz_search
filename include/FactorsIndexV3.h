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
#include "HashTrie.h"
#include "KarpRabinFactorsSuffixes.h"

using namespace sdsl;
using namespace std;

class FactorsIndexV3 {

private: 
	unsigned int len_text;
	const char *ref_text;
	unsigned int len_ref;
	unsigned int n_factors;
	
	csa_wt<> fm_index;
	rmq_succinct_sct<false, bp_support_sada<256,32,rank_support_v5<> > > rmq;
	
	rrr_vector<127> rrr_s;
	rrr_vector<127>::select_1_type select1_s;
	rrr_vector<127>::select_0_type select0_s;
	
//	sd_vector<> rrr_b;
//	sd_vector<>::select_1_type select1_b;
//	sd_vector<>::select_0_type select0_b;
	
	rrr_vector<127> rrr_b;
	rrr_vector<127>::select_1_type select1_b;
	rrr_vector<127>::select_0_type select0_b;
	
	inv_perm_support<> perm_inv;
	inv_perm_support<> perm;
	
	inv_perm_support<> perm_x;
	inv_perm_support<> perm_y;
	inv_perm_support<> perm_y_inv;
	
	wt_int<rrr_vector<63>> wt;
	
	// Prueba de aceleracion de recursive_rmq almacenando los datos de los factores descomprimidos
	bool acelerar_rmq;
	vector<unsigned int> arr_tu;
	vector<unsigned int> arr_pu;
	vector<unsigned int> arr_lu;
	
	KarpRabin *karp_rabin;
	vector<unsigned long long> arr_kr_ref;
	vector<unsigned long long> arr_kr_s;
	HashTrie tree_y;
	HashTrieRev tree_x;
	KarpRabinFactorsSuffixes *kr_factors;
	
	// Solo para pruebas
	bit_vector arr_s;
	bit_vector arr_b;
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
//	FactorsIndexV3(vector<pair<unsigned int, unsigned int> > &factors, char *full_text, unsigned int _len_text, const char *_ref_text, unsigned int _len_ref, KarpRabin *_karp_rabin, const char *kr_frases_file, bool load_kr_frases = false);
	FactorsIndexV3(vector<pair<unsigned int, unsigned int> > &factors, char *full_text, unsigned int _len_text, const char *_ref_text, unsigned int _len_ref, KarpRabin *_karp_rabin, const char *index_base_file);
	~FactorsIndexV3();
	
	void find(const string &pattern, vector<unsigned int> &results);
	void findTimes(const string &pattern, vector<unsigned int> &results);
	
	unsigned long long querytime_p1;
	unsigned long long querytime_p2;
	unsigned long long querytime_p3x;
	unsigned long long querytime_p3y;
	unsigned long long querytime_p4;
	
	KarpRabinFactorsSuffixes *getKRFactors(){
		return kr_factors;
	}
	
};

#endif //_FACTORS_INDEX_H
