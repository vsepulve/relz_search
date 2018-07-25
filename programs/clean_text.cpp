#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

#include <map>
#include <vector>

using namespace std;

int main(int argc, char* argv[]){

	if(argc != 2){
		cout<<"\nUsage: clean_text input\n";
		cout<<"This program writes the clean text to cout\n";
		return 0;
	}
	
	const char *input = argv[1];
	
//	cout << "Start";
	
	unsigned int max_read = 10000000;
	char *buff = new char[max_read + 1];
	char *buff_output = new char[max_read + 1];
	memset(buff_output, 0, max_read);
	
	fstream reader(input, fstream::in);
	unsigned int cur_pos = 0;
	while( reader.good() ){
		reader.read(buff, max_read);
		unsigned int n_read = reader.gcount();
		if( n_read < 1 ){
			continue;
		}
		buff[n_read] = 0;
		for(unsigned int i = 0; i < n_read; ++i){
			char c = toupper(buff[i]);
			if(c == 'A' || c == 'T' || c == 'C' || c == 'G'){
				buff_output[cur_pos++] = buff[i];
				if( cur_pos == max_read ){
					buff_output[cur_pos] = 0;
					cout << buff_output;
					cur_pos = 0;
					memset(buff_output, 0, max_read);
				}
			}
		}
	}
	if( cur_pos > 0 ){
		buff_output[cur_pos] = 0;
		cout << buff_output;
	}
	// cout << "\n";
	reader.close();
	
	delete [] buff;
	delete [] buff_output;
	
//	cout << "End\n";
	
}

















