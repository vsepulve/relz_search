#ifndef _FACTORS_ITERATOR_H
#define _FACTORS_ITERATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <algorithm>
#include <vector>

#include "FactorsConfig.h"
#include "CompactedText.h"

using namespace sdsl;
using namespace std;

class FactorsIteratorCompacted {

private: 
	
	bits_s_type::select_1_type *select1_s;
	bits_b_type::select_1_type *select1_b;
	bits_b_type::select_0_type *select0_b;
	int_vector<> *pi_inv;
	
	CompactedText *ref_text;
	
	// Largo de la coleccion completa comprimida
	unsigned int full_size;
	
	// Posicion actual en la referencia para el factor actual
	unsigned int cur_pos;
	
	// Largo y posicion actual en el texto comprimido total 
	// Es decir, su largo y posicion como si fuese un string independiente
	// Estos datos pueden ser usados para simplificar comparaciones de strings
	unsigned int text_length;
	unsigned int text_pos;
	unsigned int max_length;
	
	// Factor actual, con posiciones para la referencia
	unsigned int start_f;
	unsigned int n_factors;
	unsigned int cur_f;
	unsigned int cur_f_ini;
	unsigned int cur_f_fin;
	
	void loadFactor(unsigned int f, bool reset = false);
	
public: 
	
	FactorsIteratorCompacted();
	
	FactorsIteratorCompacted( unsigned int _start_f, unsigned int _n_factors, 
			bits_s_type::select_1_type *_select1_s, 
			bits_b_type::select_1_type *_select1_b, 
			bits_b_type::select_0_type *_select0_b, 
			int_vector<> *_pi_inv,
			CompactedText *_ref_text,
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

class FactorsIteratorCompactedReverse {

private: 
	
	bits_s_type::select_1_type *select1_s;
	bits_b_type::select_1_type *select1_b;
	bits_b_type::select_0_type *select0_b;
	int_vector<> *pi_inv;
	
	CompactedText *ref_text;
	
	// Largo de la coleccion completa comprimida
	unsigned int full_size;
	
	// Posicion actual en la referencia para el factor actual
	unsigned int cur_pos;
	
	// Largo y posicion actual en el texto comprimido total 
	// Es decir, su largo y posicion como si fuese un string independiente
	// Estos datos pueden ser usados para simplificar comparaciones de strings
	unsigned int text_length;
	unsigned int text_pos;
	unsigned int max_length;
	
	// Factor actual, con posiciones para la referencia
	unsigned int start_f;
	unsigned int n_factors;
	unsigned int cur_f;
	unsigned int cur_f_ini;
	unsigned int cur_f_fin;
	
	void loadFactor(unsigned int f, bool reset = false);
	
public: 
	
	FactorsIteratorCompactedReverse();
	
	FactorsIteratorCompactedReverse( unsigned int _start_f, unsigned int _n_factors, 
			bits_s_type::select_1_type *_select1_s, 
			bits_b_type::select_1_type *_select1_b, 
			bits_b_type::select_0_type *_select0_b, 
			int_vector<> *_pi_inv,
			CompactedText *_ref_text,
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
