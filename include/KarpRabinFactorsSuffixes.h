#ifndef _KARP_RABIN_FACTORS_SUFFIXES_H
#define _KARP_RABIN_FACTORS_SUFFIXES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <iostream>
#include <sstream>
#include <fstream>

#include <vector>

#include <sdsl/inv_perm_support.hpp>

#include "KarpRabin.h"

using namespace std;
using namespace sdsl;

class KarpRabinFactorsSuffixes{

protected: 
	unsigned int n_factors;
	// Estructura con los hash de cada prefijo, ordenado posicionalmente por factor
	// Asi, el primer hash es 0, y el ultimmo (en la posicion *n_factors*) tiene el hash de la coleccion completa
	vector<unsigned long long> *arr_kr_s;
	KarpRabin *karp_rabin;
	
	// Reference text
	// This should be replaced by the KarpRabinReference object
	const char *ref_text;
	
	// Permutation to turn posicional factors to its id in the other data structures
	inv_perm_support<> *perm_inv;
	
	// Si usamos la aceleracion de factores
	vector<unsigned int> *arr_tu;
	vector<unsigned int> *arr_pu;
	vector<unsigned int> *arr_lu;
	
public: 
	
	KarpRabinFactorsSuffixes();
	KarpRabinFactorsSuffixes(unsigned int _n_factors, vector<unsigned long long> *_arr_kr_s, KarpRabin *_karp_rabin, const char *_ref_text, inv_perm_support<> *_perm_inv, vector<unsigned int> *_arr_tu, vector<unsigned int> *_arr_pu, vector<unsigned int> *_arr_lu);
	virtual ~KarpRabinFactorsSuffixes();
	
	// Returns the hash from a suffix from a factor (positional), offset and length
	unsigned long long hash(unsigned int factor_ini, unsigned int offset, unsigned int length);
	
	// REturns the hash from the start of factor_ini (positional) of length length
	unsigned long long hash(unsigned int factor_ini, unsigned int length);
	
	// Returns the hash from a suffix from a factor (positional), offset and length
	// Fast version, only valid for offset and length < karp_rabin->getTableSize()
	unsigned long long hashFast(unsigned int factor_ini, unsigned int offset, unsigned int length);
	
	// REturns the hash from the start of factor_ini (positional) of length length
	// Fast version, only valid for offset and length < karp_rabin->getTableSize()
	unsigned long long hashFast(unsigned int factor_ini, unsigned int length);
	
};







#endif //_TEXT_FILTER_H




