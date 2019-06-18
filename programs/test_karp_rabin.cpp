#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

#include <algorithm>
#include <vector>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>
#include <sdsl/inv_perm_support.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bp_support_g.hpp>
#include <sdsl/bp_support_gg.hpp>

using namespace std;
using namespace sdsl;

#include "KarpRabin.h"
#include "FactorsConfig.h"
#include "KarpRabinFactorsSuffixes.h"

int main(int argc, char* argv[]){
	
	cout << "Start\n";
	
	unsigned int bits = 8;
//	unsigned int mod = 787;
	unsigned int mod = 15485863;

	if(argc == 2){
		mod = atoi(argv[1]);
	}
	
	cout << "Preparing KarpRabin(bits: " << bits << ", mod: " << mod << ")\n";
	KarpRabin karp_rabin(bits, mod, 10000000);
	
	while(true){
		cout << "Enter a word (void word to terminate)\n";
		string word;
//		cin >> word;
		getline(cin, word);
		if( word.length() < 1 ){
			break;
		}
		unsigned long long hash = karp_rabin.hash(word);
		cout << "Hash \"" << word << "\" : " << hash << "\n";
	}

	cout << "End\n";
}




















