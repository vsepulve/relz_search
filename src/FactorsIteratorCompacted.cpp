#include "FactorsIteratorCompacted.h"

// Metodos de FactorsIteratorCompacted Normal

void FactorsIteratorCompacted::loadFactor(unsigned int f, bool reset){
//	cout << "FactorsIteratorCompacted::loadFactor - Cargando factor " << f << "\n";
	cur_f = f;
	if( cur_f == (unsigned int)(-1) || cur_f >= n_factors){
		return;
	}
	// Convertir el factor posicional creciente a la posicion EN la permutacion
	
	unsigned int cur_pi = (*pi_inv)[cur_f];
	unsigned int tu = select1_s->operator()(cur_pi + 1) - cur_pi;
	unsigned int pu = select1_b->operator()(cur_f + 1);
	unsigned int lu = select1_b->operator()(cur_f + 2) - pu;
	
//	cout << "FactorsIteratorCompacted::loadFactor - tu: " << tu << ", pu: " << pu << ", lu: " << lu << "\n";
	cur_f_ini = tu;
	cur_f_fin = tu + lu - 1;
	cur_pos = cur_f_ini;
	// valores como string independiente
	if(reset){
		text_length = full_size - pu;
		text_pos = 0;
	}
}

FactorsIteratorCompacted::FactorsIteratorCompacted(){
	start_f = 0;
	n_factors = 0;
	select1_s = NULL;
	select1_b = NULL;
	select0_b = NULL;
	pi_inv = NULL;
	ref_text = NULL;
	full_size = 0;
	cur_pos = 0;
	cur_f = 0;
	cur_f_ini = 0;
	cur_f_fin = 0;
	text_length = 0;
	text_pos = 0;
	max_length = 100;
	reset();
}

FactorsIteratorCompacted::FactorsIteratorCompacted( unsigned int _start_f, unsigned int _n_factors, 
		bits_s_type::select_1_type *_select1_s, 
		bits_b_type::select_1_type *_select1_b, 
		bits_b_type::select_0_type *_select0_b,
		int_vector<> *_pi_inv,
		CompactedText *_ref_text,
		unsigned int _full_size ){
	start_f = _start_f;
	n_factors = _n_factors;
	select1_s = _select1_s;
	select1_b = _select1_b;
	select0_b = _select0_b;
	pi_inv = _pi_inv;
	ref_text = _ref_text;
	full_size = _full_size;
	cur_pos = 0;
	cur_f = 0;
	cur_f_ini = 0;
	cur_f_fin = 0;
	max_length = 100;
	reset();
}

void FactorsIteratorCompacted::reset(){
//	cout << "FactorsIteratorCompacted::reset\n";
	loadFactor( start_f, true );
}

char FactorsIteratorCompacted::next(){
//	cout << "FactorsIteratorCompacted::next - cur_pos: " << cur_pos << ", cur_f_fin: " << cur_f_fin << ", cur_f: " << cur_f << " / " << n_factors << "\n";
	char ret = 0;
	ret = (*ref_text)[cur_pos];
	++cur_pos;
	++text_pos;
	if( (cur_pos > cur_f_fin) && (cur_f < n_factors-1) ){
		loadFactor(++cur_f);
	}
	return ret;
}

bool FactorsIteratorCompacted::hasNext(){
//	cout << "FactorsIteratorCompacted::hasNext - " << cur_pos << " <= " << cur_f_fin << "?\n";

	if( text_pos < text_length ){
		return true;
	}
	return false;
}

unsigned int FactorsIteratorCompacted::length(){
	return text_length;
}

unsigned int FactorsIteratorCompacted::position(){
	return text_pos;
}

// Metodos de FactorsIteratorCompactedReverse

void FactorsIteratorCompactedReverse::loadFactor(unsigned int f, bool reset){
//	cout << "FactorsIteratorCompactedReverse::loadFactor - Cargando factor " << f << "\n";
	cur_f = f;
	if( cur_f == (unsigned int)(-1) || cur_f >= n_factors){
//		cout << "FactorsIteratorCompactedReverse::loadFactor - Factor invalido\n";
		return;
	}
	// Convertir el factor posicional creciente a la posicion EN la permutacion 
	
	unsigned int cur_pi = (*pi_inv)[cur_f];
	unsigned int tu = select1_s->operator()(cur_pi + 1) - cur_pi;
	unsigned int pu = select1_b->operator()(cur_f + 1);
	unsigned int lu = select1_b->operator()(cur_f + 2) - pu;
	
//	cout << "FactorsIteratorCompactedReverse::loadFactor - cur_pi: " << cur_pi << ", tu: " << tu << ", pu: " << pu << ", lu: " << lu << "\n";
	cur_f_ini = tu;
	cur_f_fin = tu + lu - 1;
	cur_pos = cur_f_fin;
	// valores como string independiente
	if(reset){
		// Notar que este texto INCLUYE al factor actual
		// El llamador obviamente pide el prefijo de f inverso como (f-1)
		text_length = lu;
		text_pos = 0;
	}
}

FactorsIteratorCompactedReverse::FactorsIteratorCompactedReverse(){
	start_f = 0;
	n_factors = 0;
	select1_s = NULL;
	select1_b = NULL;
	select0_b = NULL;
	pi_inv = NULL;
	ref_text = NULL;
	full_size = 0;
	cur_pos = 0;
	cur_f = 0;
	cur_f_ini = 0;
	cur_f_fin = 0;
	text_length = 0;
	text_pos = 0;
	max_length = 100;
	reset();
}

FactorsIteratorCompactedReverse::FactorsIteratorCompactedReverse( unsigned int _start_f, unsigned int _n_factors, 
		bits_s_type::select_1_type *_select1_s, 
		bits_b_type::select_1_type *_select1_b, 
		bits_b_type::select_0_type *_select0_b, 
		int_vector<> *_pi_inv,
		CompactedText *_ref_text,
		unsigned int _full_size ){
	start_f = _start_f;
	n_factors = _n_factors;
	select1_s = _select1_s;
	select1_b = _select1_b;
	select0_b = _select0_b;
	pi_inv = _pi_inv;
	ref_text = _ref_text;
	full_size = _full_size;
	cur_pos = 0;
	cur_f = 0;
	cur_f_ini = 0;
	cur_f_fin = 0;
	max_length = 100;
	reset();
}

void FactorsIteratorCompactedReverse::reset(){
//	cout << "FactorsIteratorCompactedReverse::reset\n";
	loadFactor( start_f, true );
}

char FactorsIteratorCompactedReverse::next(){
//	cout << "FactorsIteratorCompactedReverse::next - cur_pos: " << cur_pos << ", cur_f_ini: " << cur_f_ini << ", cur_f: " << cur_f << " / " << n_factors << "\n";
	char ret = 0;
	ret = (*ref_text)[cur_pos];
	--cur_pos;
	++text_pos;
	
	return ret;
}

bool FactorsIteratorCompactedReverse::hasNext(){
//	cout << "FactorsIteratorCompactedReverse::hasNext - " << text_pos << " < " << text_length << " ? \n";
	if( text_pos < text_length ){
		return true;
	}
	return false;
}

unsigned int FactorsIteratorCompactedReverse::length(){
	return text_length;
}

unsigned int FactorsIteratorCompactedReverse::position(){
	return text_pos;
}

