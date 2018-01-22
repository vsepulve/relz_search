#ifndef _TEXT_FILTER_NO_NL_H
#define _TEXT_FILTER_NO_NL_H

#include "NanoTimer.h"
#include "TextFilter.h"

using namespace std;

class TextFilterNoNL : public TextFilter{

protected: 
	
public: 
	
	TextFilterNoNL();
	
	virtual ~TextFilterNoNL();
	
	// Retorna true para los chars validos para este TextFilterNoNL (No filtrados)
	// Retoran false para todo char invalido (que sera filtrado)
	virtual bool validChar(char c);
	
	// Almacena el alfabeto valido completo en alphabet (si es NULL lo omite)
	// Retorna el total de chars validos (incluso si alphabet == NULL)
	virtual unsigned int getAlphabet(vector<char> *alphabet = NULL);
	
	// Con validChar y el TextFilter::readText basta, no es necesario redefinir en este caso
//	virtual char *readText(const char *in_file, unsigned int &text_length, vector< pair<unsigned int, unsigned int> > *lowcase_runs = NULL);
	
};







#endif //_TEXT_FILTER_NO_NL_H





