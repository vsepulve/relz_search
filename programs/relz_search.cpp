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
#include "FactorsIteratorComparator.h"

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[]){

	if(argc != 4){
		cout<<"\nModo de Uso: relz_search serialized_ref sequence_file output_relz\n";
		return 0;
	}
	
	const char *ref_file = argv[1];
	const char *input = argv[2];
	const char *output = argv[3];
	
	ReferenceIndex *reference = new ReferenceIndexBasic();
	reference->load(ref_file);
	
	// Preparar Compresor
	TextFilter *filter = new TextFilterFull();
	CompressorSingleBuffer compressor(
		output, 
		new CoderBlocksRelz(reference), 
		new DecoderBlocksRelz(reference->getText()), 
		filter
	);
	
	vector<pair<unsigned int, unsigned int> > factors;
	compressor.compress(input, 1, 1000000, 0, &factors);
	
	unsigned int len_ref = reference->getLength();
	unsigned int len_text = compressor.getTextSize();
	
	const char *ref = reference->getText();
	
	FactorsIndex index(factors, len_text, ref, len_ref);
	cout << "-----     -----\n";
	
	index.find("AB");
	cout << "-----     -----\n";
	
	index.find("ABA");
	cout << "-----     -----\n";
	
	index.find("ALA");
	cout << "-----     -----\n";
	
	
	
	
	
	delete reference;




}




















