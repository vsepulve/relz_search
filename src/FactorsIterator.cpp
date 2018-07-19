#include "FactorsIterator.h"

// Metodos de FactorsIterator Normal

void FactorsIterator::loadFactor(unsigned int f, bool reset){
//	cout << "FactorsIterator::loadFactor - Cargando factor " << f << "\n";
	cur_f = f;
	if( cur_f == (unsigned int)(-1) || cur_f >= n_factors){
		return;
	}
	// Convertir el factor posicional creciente a la posicion EN la permutacion (con perm_inv)
	unsigned int cur_perm = (*perm_inv)[cur_f];
	unsigned int tu = select1_s->operator()(cur_perm + 1) - cur_perm;
	unsigned int pu = select1_b->operator()( (*perm)[cur_perm] + 1);
	unsigned int lu = select1_b->operator()( (*perm)[cur_perm] + 2) - pu;
//	cout << "FactorsIterator::loadFactor - tu: " << tu << ", pu: " << pu << ", lu: " << lu << "\n";
	cur_f_ini = tu;
	cur_f_fin = tu + lu - 1;
	cur_pos = cur_f_ini;
	// valores como string independiente
	if(reset){
		text_length = full_size - pu;
		text_pos = 0;
	}
}

FactorsIterator::FactorsIterator(){
	start_f = 0;
	n_factors = 0;
	select1_s = NULL;
	select1_b = NULL;
	select0_b = NULL;
	perm = NULL;
	perm_inv = NULL;
	ref_text = NULL;
	full_size = 0;
	cur_pos = 0;
	cur_f = 0;
	cur_f_ini = 0;
	cur_f_fin = 0;
	reset();
}

FactorsIterator::FactorsIterator( unsigned int _start_f, unsigned int _n_factors, 
		rrr_vector<127>::select_1_type *_select1_s, 
		rrr_vector<127>::select_1_type *_select1_b, 
		rrr_vector<127>::select_0_type *_select0_b, 
		inv_perm_support<> *_perm, 
		inv_perm_support<> *_perm_inv, 
		const char *_ref_text,
		unsigned int _full_size ){
	start_f = _start_f;
	n_factors = _n_factors;
	select1_s = _select1_s;
	select1_b = _select1_b;
	select0_b = _select0_b;
	perm = _perm;
	perm_inv = _perm_inv;
	ref_text = _ref_text;
	full_size = _full_size;
	cur_pos = 0;
	cur_f = 0;
	cur_f_ini = 0;
	cur_f_fin = 0;
	reset();
}

void FactorsIterator::reset(){
//	cout << "FactorsIterator::reset\n";
	loadFactor( start_f, true );
}

char FactorsIterator::next(){
//	cout << "FactorsIterator::next - cur_pos: " << cur_pos << ", cur_f_fin: " << cur_f_fin << ", cur_f: " << cur_f << " / " << n_factors << "\n";
	char ret = ref_text[cur_pos];
	++cur_pos;
	++text_pos;
	if( (cur_pos > cur_f_fin) && (cur_f < n_factors-1) ){
		loadFactor(++cur_f);
	}
	return ret;
}

bool FactorsIterator::hasNext(){
//	cout << "FactorsIterator::hasNext - " << cur_pos << " <= " << cur_f_fin << "?\n";
//	if( cur_pos <= cur_f_fin ){
	if( cur_pos <= cur_f_fin && cur_f >= 0 && cur_f < n_factors ){
		return true;
	}
	return false;
}

unsigned int FactorsIterator::length(){
	return text_length;
}

unsigned int FactorsIterator::position(){
	return text_pos;
}

// Metodos de FactorsIteratorReverse

void FactorsIteratorReverse::loadFactor(unsigned int f, bool reset){
//	cout << "FactorsIteratorReverse::loadFactor - Cargando factor " << f << "\n";
	cur_f = f;
	if( cur_f == (unsigned int)(-1) || cur_f >= n_factors){
		return;
	}
	// Convertir el factor posicional creciente a la posicion EN la permutacion (con perm_inv)
	unsigned int cur_perm = (*perm_inv)[cur_f];
	unsigned int tu = select1_s->operator()(cur_perm + 1) - cur_perm;
	unsigned int pu = select1_b->operator()( (*perm)[cur_perm] + 1);
	unsigned int lu = select1_b->operator()( (*perm)[cur_perm] + 2) - pu;
//	cout << "FactorsIteratorReverse::loadFactor - tu: " << tu << ", pu: " << pu << ", lu: " << lu << "\n";
	cur_f_ini = tu;
	cur_f_fin = tu + lu - 1;
	cur_pos = cur_f_fin;
	// valores como string independiente
	if(reset){
		// Notar que este texto INCLUYE al factor actual
		// El llamador obviamente pide el prefijo de f inverso como (f-1)
		text_length = pu + lu;
		text_pos = 0;
	}
}

FactorsIteratorReverse::FactorsIteratorReverse(){
	start_f = 0;
	n_factors = 0;
	select1_s = NULL;
	select1_b = NULL;
	select0_b = NULL;
	perm = NULL;
	perm_inv = NULL;
	ref_text = NULL;
	full_size = 0;
	cur_pos = 0;
	cur_f = 0;
	cur_f_ini = 0;
	cur_f_fin = 0;
	reset();
}

FactorsIteratorReverse::FactorsIteratorReverse( unsigned int _start_f, unsigned int _n_factors, 
		rrr_vector<127>::select_1_type *_select1_s, 
		rrr_vector<127>::select_1_type *_select1_b, 
		rrr_vector<127>::select_0_type *_select0_b, 
		inv_perm_support<> *_perm, 
		inv_perm_support<> *_perm_inv, 
		const char *_ref_text,
		unsigned int _full_size ){
	start_f = _start_f;
	n_factors = _n_factors;
	select1_s = _select1_s;
	select1_b = _select1_b;
	select0_b = _select0_b;
	perm = _perm;
	perm_inv = _perm_inv;
	ref_text = _ref_text;
	full_size = _full_size;
	cur_pos = 0;
	cur_f = 0;
	cur_f_ini = 0;
	cur_f_fin = 0;
	reset();
}

void FactorsIteratorReverse::reset(){
//	cout << "FactorsIteratorReverse::reset\n";
	loadFactor( start_f, true );
}

char FactorsIteratorReverse::next(){
//	cout << "FactorsIteratorReverse::next - cur_pos: " << cur_pos << ", cur_f_ini: " << cur_f_ini << ", cur_f: " << cur_f << " / " << n_factors << "\n";
	char ret = ref_text[cur_pos];
	--cur_pos;
	++text_pos;
	if( (cur_pos < cur_f_ini || cur_pos == (unsigned int)(-1)) && (--cur_f != (unsigned int)(-1)) ){
		loadFactor(cur_f);
	}
	return ret;
}

bool FactorsIteratorReverse::hasNext(){
//	cout << "FactorsIteratorReverse::hasNext - " << cur_f << " != " << (unsigned int)(-1) << " (cur_pos: " << cur_pos << ")\n";
//	if( cur_f != (unsigned int)(-1) ){
	if( cur_f >= 0 && cur_f < n_factors ){
		return true;
	}
	return false;
}

unsigned int FactorsIteratorReverse::length(){
	return text_length;
}

unsigned int FactorsIteratorReverse::position(){
	return text_pos;
}

