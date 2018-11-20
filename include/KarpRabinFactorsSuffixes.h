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

#include "NanoTimer.h"
#include "KarpRabin.h"

using namespace std;
using namespace sdsl;

// For now this class is just an interface

class KarpRabinFactorsSuffixes {

protected: 
	
public: 
	
	KarpRabinFactorsSuffixes(){}
	virtual ~KarpRabinFactorsSuffixes(){}
	
	// Returns the hash from a suffix from a factor (positional), offset and length
	virtual unsigned long long hash(unsigned int factor_ini, unsigned int offset, unsigned int length){
		cerr << "KarpRabinFactorsSuffixes::hash - Not Implemented.\n";
		return 0;
	}
	
	// REturns the hash from the start of factor_ini (positional) of length length
	virtual unsigned long long hash(unsigned int factor_ini, unsigned int length){
		cerr << "KarpRabinFactorsSuffixes::hash - Not Implemented.\n";
		return 0;
	}
	
	// Returns the hash from a suffix from a factor (positional), offset and length
	// Fast version, only valid for offset and length < karp_rabin->getTableSize()
	virtual unsigned long long hashFast(unsigned int factor_ini, unsigned int offset, unsigned int length){
		cerr << "KarpRabinFactorsSuffixes::hash - Not Implemented.\n";
		return 0;
	}
	
	// REturns the hash from the start of factor_ini (positional) of length length
	// Fast version, only valid for offset and length < karp_rabin->getTableSize()
	virtual unsigned long long hashFast(unsigned int factor_ini, unsigned int length){
		cerr << "KarpRabinFactorsSuffixes::hash - Not Implemented.\n";
		return 0;
	}
	
	unsigned long long kr_nano;
	unsigned int max_offset;
	unsigned int max_length;
	
};




#endif 





