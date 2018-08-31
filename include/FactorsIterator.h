#ifndef _FACTORS_ITERATOR_H
#define _FACTORS_ITERATOR_H

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

#include "FactorsConfig.h"

using namespace sdsl;
using namespace std;

class FactorsIterator {

private: 
	
	// Valores y estructuras globales, compartidas entre iteradores
	bits_s_type::select_1_type *select1_s;
	
	bits_b_type::select_1_type *select1_b;
	bits_b_type::select_0_type *select0_b;

	inv_perm_support<> *perm;
	inv_perm_support<> *perm_inv;
	
	// Permutaciones unc
	int_vector<> *pi;
	int_vector<> *pi_inv;
	
	// Texto de la referencia descomprimido
	const char *ref_text;
	
	// FM-Index en lugar del texto
	fm_index_type *fm_index;
	
	// Largo de la coleccion completa comprimida
	unsigned int full_size;
	
	// Posicion actual en la referencia para el factor actual
	unsigned int cur_pos;
	
	// Largo y posicion actual en el texto comprimido total 
	// Es decir, su largo y posicion como si fuese un string independiente
	// Estos datos pueden ser usados para simplificar comparaciones de strings
	unsigned int text_length;
	unsigned int text_pos;
	string str_f;
	unsigned int max_length;
	
	// Factor actual, con posiciones para la referencia
	unsigned int start_f;
	unsigned int n_factors;
	unsigned int cur_f;
	unsigned int cur_f_ini;
	unsigned int cur_f_fin;
	
	void loadFactor(unsigned int f, bool reset = false);
	
public: 
	
	FactorsIterator();
	
	FactorsIterator( unsigned int _start_f, unsigned int _n_factors, 
			bits_s_type::select_1_type *_select1_s, 
			bits_b_type::select_1_type *_select1_b, 
			bits_b_type::select_0_type *_select0_b, 
			inv_perm_support<> *_perm, 
			inv_perm_support<> *_perm_inv, 
			int_vector<> *_pi,
			int_vector<> *_pi_inv,
			const char *_ref_text,
			fm_index_type *_fm_index,
			unsigned int _full_size );
	
	void reset();
	
	char next();
	
	bool hasNext();
	
	unsigned int length();
	
	unsigned int position();
	
	void setMaxLength(unsigned int _max_length){
		max_length = _max_length;
	}
	
};

class FactorsIteratorReverse {

private: 
	
	// Valores y estructuras globales, compartidas entre iteradores
	bits_s_type::select_1_type *select1_s;
	
	bits_b_type::select_1_type *select1_b;
	bits_b_type::select_0_type *select0_b;

	inv_perm_support<> *perm;
	inv_perm_support<> *perm_inv;
	
	// Permutaciones unc
	int_vector<> *pi;
	int_vector<> *pi_inv;
	
	const char *ref_text;
	
	// FM-Index en lugar del texto
	fm_index_type *fm_index;
	
	// Largo de la coleccion completa comprimida
	unsigned int full_size;
	
	// Posicion actual en la referencia para el factor actual
	unsigned int cur_pos;
	
	// Largo y posicion actual en el texto comprimido total 
	// Es decir, su largo y posicion como si fuese un string independiente
	// Estos datos pueden ser usados para simplificar comparaciones de strings
	unsigned int text_length;
	unsigned int text_pos;
	string str_f;
	unsigned int max_length;
	
	// Factor actual, con posiciones para la referencia
	unsigned int start_f;
	unsigned int n_factors;
	unsigned int cur_f;
	unsigned int cur_f_ini;
	unsigned int cur_f_fin;
	
	void loadFactor(unsigned int f, bool reset = false);
	
public: 
	
	FactorsIteratorReverse();
	
	FactorsIteratorReverse( unsigned int _start_f, unsigned int _n_factors, 
			bits_s_type::select_1_type *_select1_s, 
			bits_b_type::select_1_type *_select1_b, 
			bits_b_type::select_0_type *_select0_b, 
			inv_perm_support<> *_perm, 
			inv_perm_support<> *_perm_inv, 
			int_vector<> *_pi,
			int_vector<> *_pi_inv,
			const char *_ref_text,
			fm_index_type *_fm_index,
			unsigned int _full_size );
	
	void reset();
	
	char next();
	
	bool hasNext();
	
	unsigned int length();
	
	unsigned int position();
	
	void setMaxLength(unsigned int _max_length){
		max_length = _max_length;
	}
	
};


#endif //_FACTORS_ITERATOR_H
