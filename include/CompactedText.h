#ifndef _COMPACTED_TEXT_H
#define _COMPACTED_TEXT_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <algorithm>
#include <vector>

#include "FactorsConfig.h"

using namespace sdsl;
using namespace std;

// This version only stores ACGT as 0123, 2 bits per char
class CompactedText {

private: 
	
	// in chars
	unsigned int len;
	// 1 + len / 4 
	unsigned char *bytes;
	
	static constexpr char table[] = "ACGT";
	
public: 
	
	CompactedText();
	
	CompactedText(const char *text, unsigned long long _len);
	
	~CompactedText();
	
	char operator[](unsigned long long pos);
	
	char at(unsigned long long pos);
	
	unsigned long long length();
	
	void save(const string &file);
	
	void load(const string &file);
	
	void save(fstream *writer);
	
	void load(fstream *reader);
	
};



#endif //_COMPACTED_TEXT_H
