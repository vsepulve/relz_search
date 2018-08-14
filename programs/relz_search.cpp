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

	if(argc != 5){
		cout<<"\nUsage: relz_search serialized_ref sequence_file output_relz queries_file\n";
		return 0;
	}
	
	const char *ref_file = argv[1];
	const char *input = argv[2];
	const char *output = argv[3];
	const char *queries_file = argv[4];
	
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
	
	// Recargo el texto de la entrada, solo para acelerar la construccion
	unsigned long long real_len_text = 0;
	char *text = filter->readText(input, real_len_text);
	unsigned int len_text = compressor.getTextSize();
	cout << "Full text loaded of " << len_text << " / " << real_len_text << " chars\n";
	
	const char *ref = reference->getText();
	unsigned int len_ref = reference->getLength();
	
	cout << "----- Building index -----\n";
	NanoTimer timer;
	vector<unsigned int> results;
	bool omitir_texto = true;
	FactorsIndex index(factors, text, len_text, ref, len_ref, omitir_texto);
	cout << "----- index finished in " << timer.getMilisec() << " ms -----\n";
	index.printSize();
	
//	cout << "----- Testing save -----\n";
//	index.save("test/index");
//	
//	cout << "----- Testing load -----\n";
//	FactorsIndex index2;
//	index2.load("test/index");
//	index2.printSize();
//	
//	cout << "----- Test Query -----\n";
//	index.find("CATC", results);
//	cout << "-----     -----\n";
//	results.clear();
	
//	cout << "----- Test Query -----\n";
//	index.find("BA", results);
//	cout << "-----     -----\n";
//	results.clear();
	
	cout << "----- Loading Queries from \"" << queries_file << "\" -----\n";
	vector<string> queries;
	unsigned int max_line = 1000000;
	char buff[max_line + 1];
	fstream reader(queries_file, fstream::in);
	while( reader.good() ){
		reader.getline(buff, max_line);
		unsigned int n_read = reader.gcount();
		if( n_read < 1 ){
			continue;
		}
		buff[n_read] = 0;
		string query(buff);
//		cout << "Query[" << queries.size() << "]: \"" << query << "\"\n";
		queries.push_back(query);
	}
	
	cout << "----- Searching Queres -----\n";
	timer.reset();
	unsigned long long total_occ = 0;
	index.querytime_p1 = 0;
	index.querytime_p2 = 0;
	index.querytime_p3 = 0;
	index.querytime_p4 = 0;
	for( string query : queries ){
		cout << "----- Query \"" << query << "\" -----\n";
//		index.find(query, results);
		index.findTimes(query, results);
		cout << "-----     -----\n";
		total_occ += results.size();
		results.clear();
	}
	cout << "----- Queries finished in " << timer.getMilisec() << " ms (" << (timer.getMilisec() / total_occ) << " ms/occ, " << queries.size() << " queries, " << total_occ << " occs) -----\n";
	
	unsigned long long total_nano = index.querytime_p1 + index.querytime_p2 + index.querytime_p3 + index.querytime_p4;
	cout << "Nanosec fm_index: " << index.querytime_p1 << " (" << ((long double)index.querytime_p1/total_nano)*100 << " \%)\n";
	cout << "Nanosec recursive_rmq: " << index.querytime_p2 << " (" << ((long double)index.querytime_p2/total_nano)*100 << " \%)\n";
	cout << "Nanosec getRange: " << index.querytime_p3 << " (" << ((long double)index.querytime_p3/total_nano)*100 << " \%)\n";
	cout << "Nanosec wt: " << index.querytime_p4 << " (" << ((long double)index.querytime_p4/total_nano)*100 << " \%)\n";
	cout << "Milisec total: " << (total_nano)/(1000000) << "\n";
	
	
	delete reference;




}




















