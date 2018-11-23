#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

#include <algorithm>
#include <vector>

using namespace std;

#include "CompactedText.h"

int main() {
	
	cout << "Start\n";
	
	string s = "";
	
	vector<char> table;
	table.push_back('A');
	table.push_back('C');
	table.push_back('G');
	table.push_back('T');
	unsigned int text_size = 1000000;
//	unsigned int text_size = 1023756213;
	
	random_device seed;
	mt19937 generator(seed());
	uniform_int_distribution<unsigned int> dist_table(0, table.size() - 1);
	for(unsigned int i = 0; i < text_size; ++i){
		s += table[ dist_table(generator) ];
	}
	
	CompactedText text(s.c_str(), s.length());
	
//	cout << "text[3]: " << text[3] << "\n";
	
	for(unsigned int i = 0; i < 10 && i < s.length(); ++i){
		cout << "Text[" << i << "]: " << text[i] << " / " << s[i] << "\n";
	}
	
	cout << "Checking...\n";
	if( text.length() != s.length() ){
		cerr << "Error, different sizes (" << text.length() << " / " << s.length() << ")\n";
	}
	
	for(unsigned int i = 0; i < text.length(); ++i){
		if( text[i] != s[i] ){
			cerr << "Error (" << text[i] << " / " << s[i] << ")\n";
		}
	}
	
	cout << "Testing save / load\n";
	text.save("./test_tmp.txt");
	
	CompactedText text2;
	text2.load("./test_tmp.txt");
	
	cout << "Checking...\n";
	if( text2.length() != s.length() ){
		cerr << "Error, different sizes (" << text2.length() << " / " << s.length() << ")\n";
	}
	
	for(unsigned int i = 0; i < text2.length(); ++i){
		if( text2[i] != s[i] ){
			cerr << "Error (" << text2[i] << " / " << s[i] << ")\n";
		}
	}
	
	cout << "End\n";
}




















