#include "FactorsFastIteratorComparator.h"

FactorsFastIteratorComparator::FactorsFastIteratorComparator(){
	factors_starts = NULL;
	full_text = NULL;
	full_size = 0;
}

FactorsFastIteratorComparator::FactorsFastIteratorComparator(vector<unsigned int> *_factors_starts, char *_full_text, unsigned int _full_size) {
	factors_starts = _factors_starts;
	full_text = _full_text;
	full_size = _full_size;
}

FactorsFastIteratorReverseComparator::FactorsFastIteratorReverseComparator(){
	factors_starts = NULL;
	full_text = NULL;
	full_size = 0;
}

FactorsFastIteratorReverseComparator::FactorsFastIteratorReverseComparator(vector<unsigned int> *_factors_starts, char *_full_text, unsigned int _full_size) {
	factors_starts = _factors_starts;
	full_text = _full_text;
	full_size = _full_size;
}


