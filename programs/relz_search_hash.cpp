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
#include "FactorsIterator.h"

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[]){

	if(argc != 3){
		cout<<"\nUsage: ./bin/relz_search_hash index_directory queries_file\n";
		return 0;
	}
	
	const char *index_directory = argv[1];
	const char *queries_file = argv[2];
	
	string output_path(index_directory);
	if( output_path.back() != '/' ){
		output_path += "/";
	}
	output_path += "index_hash";
	
	cout << "----- Testing load -----\n";
	NanoTimer timer;
	vector<unsigned int> results;
	unsigned int bits = 8;
//	unsigned int mod = 787;
	unsigned int mod = 15485863;
//	KarpRabin karp_rabin(bits, mod, 1100000000);
	KarpRabin karp_rabin(bits, mod, 10000000);
	RelzIndexHash index;
	index.load(output_path, &karp_rabin);
	index.printSize();
	
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
	index.querytime_p3x = 0;
	index.querytime_p3y = 0;
	index.querytime_p4 = 0;
	index.occs_a = 0;
	index.occs_b = 0;
	index.occs_c = 0;
	for( string query : queries ){
//		cout << "----- Query \"" << query << "\" -----\n";
//		index.find(query, results);
		index.findTimes(query, results);
//		cout << "-----     -----\n";
		total_occ += results.size();
		results.clear();
	}
	double total_milisec = timer.getMilisec();
	cout << "----- Queries finished in " << total_milisec << " ms (" << (total_milisec / total_occ) << " ms/occ, " << queries.size() << " queries, " << total_occ << " occs) -----\n";
	
	cout << "Occs A: " << index.occs_a << "\n";
	cout << "Occs B: " << index.occs_b << "\n";
	cout << "Occs C: " << index.occs_c << "\n";
//	double total_occ_ref = index.occs_a + index.occs_b + index.occs_c;
	unsigned long long total_nano = index.querytime_p1 + index.querytime_p2 + index.querytime_p3x + index.querytime_p3y + index.querytime_p4;
	
	cout << "Fm_index: " << ((long double)(index.querytime_p1))/(total_occ * 1000) << " microsec/occ (" << ((long double)index.querytime_p1/total_nano)*100 << " \%)\n";
	cout << "Recursive RMQ: " << ((long double)(index.querytime_p2))/(total_occ * 1000) << " microsec/occ (" << ((long double)index.querytime_p2/total_nano)*100 << " \%)\n";
	cout << "getRange: " << ((long double)(index.querytime_p3x + index.querytime_p3y))/(total_occ * 1000) << " microsec/occ (" << ((long double)(index.querytime_p3x + index.querytime_p3y)/total_nano)*100 << " \%)\n";
	cout << "WT: " << ((long double)(index.querytime_p4))/(total_occ * 1000) << " microsec/occ (" << ((long double)index.querytime_p4/total_nano)*100 << " \%)\n";
	cout << "Milisec total: " << (total_nano)/(1000000) << "\n";
	
	cout << "kr_factors - max_offset: " <<  index.kr_factors->max_offset << ", max_length: " <<  index.kr_factors->max_length << "\n";
	cout << "KarpRabin max used length: " <<  karp_rabin.max_len << "\n";
	

}




















