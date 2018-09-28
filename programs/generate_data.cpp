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

#include "FactorsIndex.h"
#include "FactorsIterator.h"

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[]){

	if(argc != 6){
		cout<<"\nUsage: generate_data serialized_ref input_text tmp_relz total_bytes output_text\n";
		return 0;
	}
	
	const char *ref_file = argv[1];
	const char *input_text = argv[2];
	const char *tmp_relz = argv[3];
	unsigned long long total_bytes = stoll(string(argv[4]));
	const char *output_text = argv[5];
	
	ReferenceIndex *reference = new ReferenceIndexBasic();
	reference->load(ref_file);
	
	// Preparar Compresor
	TextFilter *filter = new TextFilterFull();
	CompressorSingleBuffer compressor(
		tmp_relz, 
		new CoderBlocksRelz(reference), 
		new DecoderBlocksRelz(reference->getText()), 
		filter
	);
	
	vector<pair<unsigned int, unsigned int> > factors;
	compressor.compress(input_text, 1, 1000000, 0, &factors);
	
	const char *ref = reference->getText();
	unsigned int len_ref = reference->getLength();
	
	cout << "Generating " << total_bytes << " of text from " << len_ref << " chars of reference and " << factors.size() << " factors\n";
	
	random_device seed;
	mt19937 generator(seed());
	uniform_int_distribution<unsigned int> factors_dist(0, factors.size() - 1);
	uniform_int_distribution<unsigned int> ref_dist(0, len_ref - 1);
	
	fstream writer(output_text, fstream::trunc | fstream::out);
	unsigned long long written_bytes = 0;
	while( written_bytes < total_bytes ){
		unsigned int f = factors_dist(generator);
		unsigned long long len = factors[f].second;
		unsigned int ref_ini = ref_dist(generator);
		if( total_bytes - written_bytes < len ){
			len = total_bytes - written_bytes;
		}
		if( len_ref - ref_ini < len ){
			len = len_ref - ref_ini;
		}
//		cout << "Adding factor (" << ref_ini << ", " << len << ") from original (" << factors[f].first << ", " << factors[f].second << ")\n";
		writer.write(ref + ref_ini, len);
		written_bytes += len;
	}
	writer << "\n";
	writer.close();
	
	delete reference;
	
}




















