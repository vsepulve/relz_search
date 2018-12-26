#ifndef _KARP_RABIN_FACTORS_SUFFIXES_V2_H
#define _KARP_RABIN_FACTORS_SUFFIXES_V2_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <iostream>
#include <sstream>
#include <fstream>

#include <vector>

#include "NanoTimer.h"
#include "KarpRabin.h"
#include "KarpRabinFactorsSuffixes.h"
#include "FactorsConfig.h"
#include "BitsUtils.h"

using namespace std;
using namespace sdsl;

class KarpRabinFactorsSuffixesv2 : public KarpRabinFactorsSuffixes {

protected: 
	unsigned int n_factors;
	// Estructura con los hash de cada prefijo, ordenado posicionalmente por factor
	// Asi, el primer hash es 0, y el ultimmo (en la posicion *n_factors*) tiene el hash de la coleccion completa
	vector<unsigned long long> *arr_kr_s;
	bool delete_krs;
	
	KarpRabin *karp_rabin;
	
	vector<unsigned int> *factors_start;
	bool delete_start;
	
public: 
	
	
	// Public while testing, to give easy acces to HashTrie
	
	// Reference text
	// This should be replaced by the KarpRabinReference object
	const char *ref_text;
	
	// Select structures for S and B to process factors data (pos, len)
	bits_s_type::select_1_type *select1_s;
	bits_b_type::select_1_type *select1_b;
	bits_b_type::select_0_type *select0_b;
	
	// Permutation to turn positional factors to its id
	int_vector<> *pi_inv;
	
	
	
	KarpRabinFactorsSuffixesv2();
	
	KarpRabinFactorsSuffixesv2(unsigned int _n_factors, 
			vector<unsigned long long> *_arr_kr_s, 
			KarpRabin *_karp_rabin, 
			const char *_ref_text, 
			bits_s_type::select_1_type *_select1_s, 
			bits_b_type::select_1_type *_select1_b, 
			bits_b_type::select_0_type *_select0_b, 
			int_vector<> *_pi_inv, 
			vector<unsigned int> *_factors_start);
	
	KarpRabinFactorsSuffixesv2(const string &file, 
			KarpRabin *_karp_rabin, 
			const char *_ref_text, 
			bits_s_type::select_1_type *_select1_s, 
			bits_b_type::select_1_type *_select1_b, 
			bits_b_type::select_0_type *_select0_b, 
			int_vector<> *_pi_inv);
			
	virtual ~KarpRabinFactorsSuffixesv2();
	
	// Returns the hash from a suffix from a factor (positional), offset and length
	unsigned long long hash(unsigned int factor_ini, unsigned int offset, unsigned int length);
//	unsigned long long hashBin(unsigned int factor_ini, unsigned int offset, unsigned int length);
	
	// REturns the hash from the start of factor_ini (positional) of length length
	unsigned long long hash(unsigned int factor_ini, unsigned int length);
	
	// Returns the hash from a suffix from a factor (positional), offset and length
	// Fast version, only valid for offset and length < karp_rabin->getTableSize()
//	unsigned long long hashFast(unsigned int factor_ini, unsigned int offset, unsigned int length);
	
	// REturns the hash from the start of factor_ini (positional) of length length
	// Fast version, only valid for offset and length < karp_rabin->getTableSize()
//	unsigned long long hashFast(unsigned int factor_ini, unsigned int length);
	
	void save(const string &file);
	
	void load(const string &file);
	
	unsigned long long nano1;
	unsigned long long nano2;
	unsigned long long nano3;
	unsigned long long nano4;
	
};







#endif //_TEXT_FILTER_H





