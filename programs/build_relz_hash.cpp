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

#include "ReferenceIndexBasic.h"
#include "CompressorSingleBuffer.h"
#include "CoderBlocksRelz.h"
#include "DecoderBlocksRelz.h"
#include "TextFilterFull.h"

#include "RelzIndexHash.h"

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[]){

	if(argc != 4){
		cout<<"\nUsage: ./bin/build_relz_hash serialized_ref sequence_file index_directory\n";
		return 0;
	}
	
	const char *ref_file = argv[1];
	const char *input = argv[2];
	const char *index_directory = argv[3];
	
	ReferenceIndex *reference = new ReferenceIndexBasic();
	reference->load(ref_file);
	TextFilter *filter = new TextFilterFull();
	
	// Preparar Compresor
	Compressor compressor(
		NULL, 
		new CoderBlocksRelz(reference), 
		new DecoderBlocksRelz(reference->getText()), 
		filter
	);
	vector<pair<unsigned int, unsigned int> > factors;
	unsigned long long len_text = 0;
	char *text = compressor.compressFactors(input, 1000000, len_text, &factors);
	cout << "Full text loaded of " << len_text << " chars\n";
	
	const char *ref = reference->getText();
	unsigned int len_ref = reference->getLength();
	
	cout << "----- Building index -----\n";
	NanoTimer timer;
	unsigned int bits = 8;
//	unsigned int mod = 787;
//	unsigned int mod = 15485863;
	unsigned int mod = 44939137;
//	unsigned int mod = 2050009847;
//	KarpRabin karp_rabin(bits, mod, 1100000000);
	KarpRabin karp_rabin(bits, mod, 100000000);
	RelzIndexHash index(factors, text, len_text, ref, len_ref, &karp_rabin);
	cout << "----- index finished in " << timer.getMilisec() << " ms -----\n";
	index.printSize();
	
	cout << "----- Testing save -----\n";
	string output_path(index_directory);
	if( output_path.back() != '/' ){
		output_path += "/";
	}
	output_path += "index_hash";
	
	index.save(output_path);
	
	delete reference;

}




















