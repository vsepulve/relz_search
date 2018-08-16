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

#include "FactorsIterator.h"
/*#include "FactorsIteratorComparator.h"*/
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
	
	
	/*
	template<class t_wt              = wt_huff<>,              // Wavelet tree type
		uint32_t t_dens         = 32,                     // Sample density for suffix array (SA) values
		uint32_t t_inv_dens     = 64,                     // Sample density for inverse suffix array (ISA) values
		class t_sa_sample_strat = sa_order_sa_sampling<>, // Policy class for the SA sampling.
		class t_isa_sample_strat= isa_sampling<>,         // Policy class for ISA sampling.
		class t_alphabet_strat  =                         // Policy class for the representation of the alphabet.
		typename wt_alphabet_trait<t_wt>::type
		>
	*/
   
/*	csa_wt<> fm_index;*/
	csa_wt<wt_huff<>, 16, 32, sa_order_sa_sampling<>, isa_sampling<>> fm_index;

	rmq_succinct_sct<false, bp_support_sada<256,32,rank_support_v5<> > > rmq;
	
/*	rrr_vector<127> rrr_s;*/
/*	rrr_vector<127>::select_1_type select1_s;*/
/*	rrr_vector<127>::select_0_type select0_s;*/
	
	bit_vector rrr_s;
	bit_vector::select_1_type select1_s;
	bit_vector::select_0_type select0_s;
	
/*	rrr_vector<127> rrr_b;*/
/*	rrr_vector<127>::select_1_type select1_b;*/
/*	rrr_vector<127>::select_0_type select0_b;*/
	
	sd_vector<> rrr_b;
	sd_vector<>::select_1_type select1_b;
	sd_vector<>::select_0_type select0_b;
	
	inv_perm_support<> perm;
	inv_perm_support<> perm_inv;
	
	inv_perm_support<> perm_x;
	inv_perm_support<> perm_y;
	inv_perm_support<> perm_y_inv;
	
	wt_int<rrr_vector<63>> wt;
	
	// Prueba de aceleracion de recursive_rmq almacenando los datos de los factores descomprimidos
	bool acelerar_rmq;
	vector<unsigned int> arr_tu;
	vector<unsigned int> arr_pu;
	vector<unsigned int> arr_lu;
	
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
	
	char getChar(unsigned int factor, unsigned int pos, unsigned int max_len = 100);
	
	char getCharRev(unsigned int factor, unsigned int pos, unsigned int max_len = 100);
			
	pair<unsigned int, unsigned int> getRangeY(const char *pattern);
			
	pair<unsigned int, unsigned int> getRangeX(const char *pattern);
	
	
public: 
	FactorsIndex();
	FactorsIndex(vector<pair<unsigned int, unsigned int> > &factors, char *full_text, unsigned int _len_text, const char *_ref_text, unsigned int _len_ref, bool _omit_text = false);
	~FactorsIndex();
	
	void find(const string &pattern, vector<unsigned int> &results);
	
	void findTimes(const string &pattern, vector<unsigned int> &results);
	
	void printSize();
	
	void save(const string &file_base);
	
	void load(const string &file_base);
	
	unsigned long long querytime_p1;
	unsigned long long querytime_p2;
	unsigned long long querytime_p3;
	unsigned long long querytime_p4;
	
	
	
};








#endif //_FACTORS_INDEX_H
