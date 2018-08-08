#include "FactorsIteratorComparator.h"

FactorsIteratorComparator::FactorsIteratorComparator(){
	n_factors = 0;
	select1_s = NULL;
	select1_b = NULL;
	select0_b = NULL;
	perm = NULL;
	perm_inv = NULL;
	ref_text = NULL;
	fm_index = NULL;
	full_size = 0;
}

FactorsIteratorComparator::FactorsIteratorComparator(unsigned int _n_factors, 
		rrr_vector<127>::select_1_type *_select1_s, 
		sd_vector<>::select_1_type *_select1_b, 
		sd_vector<>::select_0_type *_select0_b, 
		inv_perm_support<> *_perm, 
		inv_perm_support<> *_perm_inv, 
		const char *_ref_text, 
		csa_wt<> *_fm_index, 
		unsigned int _full_size ) {
	n_factors = _n_factors;
	select1_s = _select1_s;
	select1_b = _select1_b;
	select0_b = _select0_b;
	perm = _perm;
	perm_inv = _perm_inv;
	ref_text = _ref_text;
	fm_index = _fm_index;
	full_size = _full_size;
}

FactorsIteratorReverseComparator::FactorsIteratorReverseComparator(){
	n_factors = 0;
	select1_s = NULL;
	select1_b = NULL;
	select0_b = NULL;
	perm = NULL;
	perm_inv = NULL;
	ref_text = NULL;
	fm_index = NULL;
	full_size = 0;
}

FactorsIteratorReverseComparator::FactorsIteratorReverseComparator(unsigned int _n_factors, 
		rrr_vector<127>::select_1_type *_select1_s, 
		sd_vector<>::select_1_type *_select1_b, 
		sd_vector<>::select_0_type *_select0_b, 
		inv_perm_support<> *_perm, 
		inv_perm_support<> *_perm_inv, 
		const char *_ref_text, 
		csa_wt<> *_fm_index, 
		unsigned int _full_size ) {
	n_factors = _n_factors;
	select1_s = _select1_s;
	select1_b = _select1_b;
	select0_b = _select0_b;
	perm = _perm;
	perm_inv = _perm_inv;
	ref_text = _ref_text;
	fm_index = _fm_index;
	full_size = _full_size;
}


