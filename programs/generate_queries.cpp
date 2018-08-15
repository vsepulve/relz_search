#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

#include <map>
#include <vector>

#include <random>

using namespace std;

int main(int argc, char* argv[]){

	if(argc != 4){
		cout<<"\nUsage: generate_queries input n_queries query_length\n";
		cout<<"This program writes its output to cout\n";
		return 0;
	}
	
	const char *input = argv[1];
	unsigned int n_queries = atoi(argv[2]);
	unsigned int length = atoi(argv[3]);
	
//	cout << "Start - Reading \"" << input << "\" to extract " << n_queries << " of length " << length << "\n";
	
	unsigned int max_read = 10000000;
	char *buff = new char[max_read + 1];
	char *buff_output = new char[length + 1];
	
	// Tomo el largo para cargarlo completo
	fstream reader(input, fstream::in);
	
	reader.seekg (0, reader.end);
	unsigned long long file_size = reader.tellg();
	reader.seekg (0, reader.beg);
	
	char *text = new char[file_size + 1];
	
	// Total leido
	unsigned long long total = 0;
	// Total escrito (ya filtrado)
	unsigned long long text_size = 0;
	while( total < file_size && reader.good() ){
		reader.read(buff, max_read);
		unsigned int n_read = reader.gcount();
		total += n_read;
		for(unsigned int i = 0; i < n_read; ++i){
			char c = toupper(buff[i]);
//			if(c == 'A' || c == 'T' || c == 'C' || c == 'G'){
			if( (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') ){
				text[text_size++] = c;
			}
		}
	}
	text[text_size] = 0;
	reader.close();
	
	
	random_device seed;
	mt19937 generator(seed());
	uniform_int_distribution<unsigned long long> dist(0ULL, text_size - length - 1);
	
	for( unsigned int i = 0; i < n_queries; ++i ){
		unsigned long long pos = dist(generator);
		memcpy(buff_output, text + pos, length);
		buff_output[length] = 0;
//		cout << "query[" << i << "]: \"" << buff_output << "\"\n";
		cout << buff_output << "\n";
	}
	
	delete [] buff;
	delete [] buff_output;
	delete [] text;
	
//	cout << "End\n";
	
}

















