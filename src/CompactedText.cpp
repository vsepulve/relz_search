#include "CompactedText.h"

constexpr char CompactedText::table[];
	
CompactedText::CompactedText(){
	len = 0;
	bytes = NULL;
}

CompactedText::CompactedText(const char *text, unsigned long long _len){
	len = 0;
	bytes = NULL;
	if( text == NULL || _len == 0 ){
		return;
	}
	// En construccion pido el peor caso por simplicidad
	bytes = new unsigned char[1 + _len/4];
	memset((char*)bytes, 0, 1 + _len/4);
	unsigned int in_pos = 0;
	unsigned long long out_pos = 0;
	unsigned char value = 0;
	for(unsigned long long i = 0; i < _len; ++i){
		if(text[i] == 'A'){
			value = 0;
			++len;
		}
		else if(text[i] == 'C'){
			value = 1;
			++len;
		}
		else if(text[i] == 'G'){
			value = 2;
			++len;
		}
		else if(text[i] == 'T'){
			value = 3;
			++len;
		}
		else{
			cerr << "CompactedText - Warning, omiting invalid char (" << text[i] << ")\n";
			continue;
		}
		value <<= (in_pos * 2);
//		cout << "CompactedText - Adding " << (unsigned int)value << " (to " << (unsigned int)bytes[out_pos] << " from " << text[i] << ")\n";
		bytes[out_pos] += value;
		if( ++in_pos == 4 ){
			in_pos = 0;
			++out_pos;
		}
	}
//	unsigned long long n_bytes = len/4;
//	if( 4*n_bytes < len ){
//		++n_bytes;
//	}
//	cout << "CompactedText - len: " << len << ", n_bytes: " << n_bytes << "\n";
//	for(unsigned long long i = 0; i < n_bytes; ++i){
//		cout << "CompactedText - byte[" << i << "]: " << (unsigned int)bytes[i] << "\n";
//	}
}

CompactedText::~CompactedText(){
	len = 0;
	if(bytes != NULL){
		delete [] bytes;
		bytes = NULL;
	}
}

char CompactedText::operator[](unsigned long long pos){
//	cout << "CompactedText::operator[" << pos << " / " << len << "]\n";
	if( pos >= len ){
		return 0;
	}
//	unsigned char byte = bytes[pos/4];
//	cout << "CompactedText::operator[] - byte (start): " << (unsigned int)byte << " (byte_pos: " << (pos/4) << ", r_pos: " << (pos&0x3) << ")\n";
//	byte >>= ((pos&0x3)*2);
//	byte &= 0x3;
//	cout << "CompactedText::operator[] - byte (end): " << (unsigned int)byte << "\n";
//	return table[byte];
	return table[ ( bytes[pos>>2] >> ((pos&0x3)<<1) ) & 0x3 ];
}

char CompactedText::at(unsigned long long pos){
	if( pos >= len ){
		return 0;
	}
	return table[ ( bytes[pos>>2] >> ((pos&0x3)<<1) ) & 0x3 ];
}

unsigned long long CompactedText::length(){
	return len;
}

void CompactedText::save(const string &file){
	fstream writer(file, fstream::out | fstream::trunc);
	save(&writer);
	writer.close();
}

void CompactedText::load(const string &file){
	fstream reader(file, fstream::in);
	load(&reader);
	reader.close();
}

void CompactedText::save(fstream *writer){
	if(writer == NULL || !(writer->good()) ){
		cerr << "CompactedText::save - Error in fstream\n";
		return;
	}
//	cout << "CompactedText::save - len: " << len << "\n";
	writer->write((char*)&len, sizeof(long long));
	unsigned long long n_bytes = len/4;
	if( n_bytes*4 < len ){
		++n_bytes;
	}
//	cout << "CompactedText::save - n_bytes: " << n_bytes << "\n";
	writer->write((char*)bytes, n_bytes);
}

void CompactedText::load(fstream *reader){
	if(reader == NULL || !(reader->good()) ){
		cerr << "CompactedText::load - Error in fstream\n";
		return;
	}
	len = 0;
	if( bytes != NULL ){
		delete [] bytes;
	}
	reader->read((char*)&len, sizeof(long long));
//	cout << "CompactedText::load - len: " << len << "\n";
	unsigned long long n_bytes = len/4;
	if( n_bytes*4 < len ){
		++n_bytes;
	}
//	cout << "CompactedText::load - n_bytes: " << n_bytes << "\n";
	bytes = new unsigned char[n_bytes];
	reader->read((char*)bytes, n_bytes);
}










