#ifndef _FACTORS_ITERATOR_COMPARATOR_H
#define _FACTORS_ITERATOR_COMPARATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <algorithm>
#include <vector>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/inv_perm_support.hpp>

#include "FactorsIterator.h"

using namespace sdsl;
using namespace std;

class FactorsIteratorComparator : public std::binary_function<unsigned int, unsigned int, bool> {
private:
	unsigned int n_factors;
	rrr_vector<127>::select_1_type *select1_s;
	rrr_vector<127>::select_1_type *select1_b;
	rrr_vector<127>::select_0_type *select0_b;
	inv_perm_support<> *perm;
	inv_perm_support<> *perm_inv;
	const char *ref_text;
	unsigned int full_size;

public:
	FactorsIteratorComparator();
	
	FactorsIteratorComparator(unsigned int _n_factors, 
			rrr_vector<127>::select_1_type *_select1_s, 
			rrr_vector<127>::select_1_type *_select1_b, 
			rrr_vector<127>::select_0_type *_select0_b, 
			inv_perm_support<> *_perm, 
			inv_perm_support<> *_perm_inv, 
			const char *_ref_text, 
			unsigned int _full_size );
	
	inline bool operator()(const unsigned int a, const unsigned int b){
		FactorsIterator it_a( a, n_factors, select1_s, select1_b, select0_b, perm, perm_inv, ref_text, full_size );
		FactorsIterator it_b( b, n_factors, select1_s, select1_b, select0_b, perm, perm_inv, ref_text, full_size );
		char ch_a, ch_b;
		while( true ){
			if( ! it_a.hasNext() ){
				return true;
			}
			if( ! it_b.hasNext() ){
				return false;
			}
			ch_a = it_a.next();
			ch_b = it_b.next();
			if( ch_a < ch_b ){
				return true;
			}
			if( ch_b < ch_a){
				return false;
			}
		}
		return true;
	}
};

class FactorsIteratorReverseComparator : public std::binary_function<unsigned int, unsigned int, bool> {
private:
	unsigned int n_factors;
	rrr_vector<127>::select_1_type *select1_s;
	rrr_vector<127>::select_1_type *select1_b;
	rrr_vector<127>::select_0_type *select0_b;
	inv_perm_support<> *perm;
	inv_perm_support<> *perm_inv;
	const char *ref_text;
	unsigned int full_size;

public:
	FactorsIteratorReverseComparator();
	
	FactorsIteratorReverseComparator(unsigned int _n_factors, 
			rrr_vector<127>::select_1_type *_select1_s, 
			rrr_vector<127>::select_1_type *_select1_b, 
			rrr_vector<127>::select_0_type *_select0_b, 
			inv_perm_support<> *_perm, 
			inv_perm_support<> *_perm_inv, 
			const char *_ref_text, 
			unsigned int _full_size );
	
	inline bool operator()(const unsigned int a, const unsigned int b){
	FactorsIteratorReverse it_a( a - 1, n_factors, select1_s, select1_b, select0_b, perm, perm_inv, ref_text, full_size );
	FactorsIteratorReverse it_b( b - 1, n_factors, select1_s, select1_b, select0_b, perm, perm_inv, ref_text, full_size );
	char ch_a, ch_b;
	while( true ){
		if( ! it_a.hasNext() ){
			return true;
		}
		if( ! it_b.hasNext() ){
			return false;
		}
		ch_a = it_a.next();
		ch_b = it_b.next();
		if( ch_a < ch_b ){
			return true;
		}
		if( ch_b < ch_a){
			return false;
		}
	}
	return true;
}
};

#endif //_FACTORS_ITERATOR_COMPARATOR_H
