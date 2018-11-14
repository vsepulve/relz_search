#include "HashTrie.h"

HashTrie::HashTrie(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

HashTrie::HashTrie(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	build(full_text, len_text, factors_start, arr_y, _karp_rabin, _kr_factors);
}

HashTrie::~HashTrie(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

void HashTrie::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	
	cout << "HashTrie::build - Start (full text of " << len_text << ", " << factors_start.size() << " factors)\n";
	
	root.build(full_text, len_text, factors_start, arr_y, karp_rabin, kr_factors, 0, factors_start.size()-1, 0);
	
	cout << "HashTrie::build - End\n";
}

HashTrieNode::HashTrieNode(){
	len = 0;
	min = 0;
	max = 0;
	hash = 0;
	min_factor_pos = 0;
}

HashTrieNode::~HashTrieNode(){

}

void HashTrieNode::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len){
	if( min == max ){
//		cout << "HashTrieNode::build - Leaf\n";
		return;
	}
	if( full_text == NULL || processed_len == len_text ){
//		cout << "HashTrieNode::build - Null text\n";
		return;
	}
	
//	cout << "HashTrieNode::build - Start (full text of " << len_text << ", range[" << min << ", " << max << "], processed_len: " << processed_len << ")\n";
	
	unsigned int min_pos = min;
	unsigned int min_text_start = factors_start[ arr_y[min_pos] ] + processed_len;
	unsigned int min_text_len = 0;
	
	unsigned int cur_pos = 0;
	unsigned int cur_text_start = 0;
	
	unsigned long long hash;
	
	while( min_pos <= max ){
		
		// Primero busco el cur_pos MAYOR que tenga *primer* char igual al de min_pos
		// Despues reviso el largo mayor
		
//		cout << "HashTrieNode::build - Iniciando busqueda binaria en [" << min_pos << ", " << max << "]\n";
		
		unsigned int l = min_pos;
		unsigned int h = max;
		while(l < h){
			cur_pos = l + ((h-l)>>1);
			cur_text_start = factors_start[ arr_y[cur_pos] ] + processed_len;
			if( full_text[cur_text_start] <= full_text[min_text_start] ){
				l = cur_pos+1;
			}
			else{
				h = cur_pos;
			}
		}
		cur_pos = h;
		cur_text_start = factors_start[ arr_y[cur_pos] ] + processed_len;
		if( (cur_pos > min_pos) && full_text[cur_text_start] != full_text[min_text_start] ){
			--cur_pos;
			cur_text_start = factors_start[ arr_y[cur_pos] ] + processed_len;
		}
//		cout << "HashTrieNode::build - cur_pos: " << cur_pos << "\n";
		
		// Si cur_pos == min_pos, el largo es min_text_len
		// Si no, buscar el largo maximo entre min y cur
		
		if( cur_pos > min_pos ){
//			cout << "HashTrieNode::build - Case 1\n";
			// Buscar el mayor
			// Aqui uso el hecho de que full_text termina en '\0'
			min_text_len = 1;
			while( full_text[min_text_start + min_text_len] == full_text[cur_text_start + min_text_len] ){
				++min_text_len;
			}
		}
		else{
//			cout << "HashTrieNode::build - Case 2\n";
			min_text_len = len_text - min_text_start;
		}
//		cout << "HashTrieNode::build - min_text_len: " << min_text_len << "\n";
//		cout << "HashTrieNode::build - Preparing string s (min_text_start: " << min_text_start << " / " << len_text << ", len: " << min_text_len << ")\n";
		
		// El string es olo para el mensaje de debug
//		string s(full_text + min_text_start, min_text_len);
		char first_char = full_text[min_text_start];
		
		// Notar que aqui necesito la posicion del factor en la coleccion, no en el arreglo Y
		hash = kr_factors->hash(arr_y[min_pos], processed_len, min_text_len);
//		cout << "HashTrieNode::build - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_text_len << ", hash: " << hash << " / " << karp_rabin->hash(s) << " (\"" << s << "\")\n";
//		cout << "HashTrieNode::build - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_text_len << ", hash: " << hash << "\n";
		// Preparar el hijo, ejecutar la llamda sobre esa instancia
		if( min_text_len == 0 ){
//			cout << "HashTrieNode::build - Omiting child of len 0\n";
		}
		else{
			childs[first_char] = std::make_shared<HashTrieNode>();
			childs[first_char]->len = min_text_len;
			childs[first_char]->min = min_pos;
			childs[first_char]->max = cur_pos;
			childs[first_char]->hash = hash;
			childs[first_char]->min_factor_pos = arr_y[min_pos];
			childs[first_char]->build(full_text, len_text, factors_start, arr_y, karp_rabin, kr_factors, min_pos, cur_pos, processed_len + min_text_len);
		}
		
//		cout << "HashTrieNode::build - Preparing min_pos = " << cur_pos+1 << " / " << arr_y.size() << "\n";
		min_pos = cur_pos+1;
		if(min_pos >= arr_y.size()){
			break;
		}
		min_text_start = factors_start[ arr_y[min_pos] ] + processed_len;
		
//		cout << "HashTrieNode::build - min_text_start: " << min_text_start << "\n";
	
	}
	
//	cout << "----- \n";
	
}

void HashTrie::print(){
	root.print(0);
	cout << "Root: " << root.childs.size() << " childs\n";
	cout << "Lengths: ";
	for( auto it_child : root.childs ){
		cout << it_child.second->len << " - ";
	}
	cout << "\n";
}

void HashTrieNode::print(unsigned int level){
	for(unsigned int i = 0; i < level; ++i){
		cout << "- ";
	}
	cout << "hash: " << hash << ", len: " << len << " , range [" << min << ", " << max << "]\n";
	for( auto it : childs ){
		cout << it.first << " ";
		it.second->print(level+1);
	}
}

unsigned int HashTrieNode::totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height){
	unsigned int ret = childs.size();
	if( height > max_height ){
		max_height = height;
	}
	if( len > max_len ){
		max_len = len;
	}
	if( childs.size() > max_childs ){
		max_childs = childs.size();
	}
	for( auto it : childs ){
		ret += it.second->totalChilds(max_len, max_childs, max_height, 1 + height);
	}
	return ret;
}

void HashTrie::printSize(){
	unsigned int max_len = 0;
	unsigned int max_childs = 0;
	unsigned int max_height = 0;
	unsigned int total_childs = root.totalChilds(max_len, max_childs, max_height, 0);
	cout << "HashTrie::printSize - totalChilds: " << total_childs << ", max_len: " << max_len << ", max_childs: " << max_childs << ", max_height: " << max_height << "\n";
}

pair<unsigned int, unsigned int> HashTrie::getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, const string &pattern){
	return root.getRange(kr_pat_vector, pos, 0, karp_rabin, kr_factors, pattern, &hash_nano);
}

pair<unsigned int, unsigned int> HashTrieNode::getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, const string &pattern, unsigned long long *hash_nano){
	
	cout << "HashTrieNode::getRange - Start (prefixes: " << kr_pat_vector.size() << ", pos: " << pos << ", processed: " << processed << ")\n";
	
	if( pos + processed >= kr_pat_vector.size() ){
		cout << "HashTrieNode::getRange - [" << min << ", " << max << "]\n";
		return pair<unsigned int, unsigned int>(min, max);
	}
	
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	unsigned int pat_len = kr_pat_vector.size() - pos - processed;
	char first_char_pat = pattern[pos + processed];
//	cout << "HashTrieNode::getRange - pat_len: " << pat_len << "\n";
	string pat = pattern.substr(pos + processed, pat_len);
	cout << "HashTrieNode::getRange - pat: " << pat << "\n";
	
	auto it_child = childs.find(first_char_pat);
	if(it_child != childs.end()){
		child_len = it_child->second->len;
		if( child_len <= pat_len ){
//			cout << "HashTrieNode::getRange - Case 1, child_len: " << child_len << "\n";
			hash_pat = karp_rabin->subtract_prefix(kr_pat_vector[pos + processed + child_len - 1], kr_pat_vector[pos + processed - 1], child_len);
			string pat_cut = pattern.substr(pos + processed, child_len);
//			cout << "HashTrieNode::getRange - pat_cut: " << pat_cut << " hash_pat: " << hash_pat << " / " << karp_rabin->hash(pat_cut) << " (processed: " << processed << ")\n";
			if( hash_pat == it_child->second->hash ){
//				cout << "HashTrieNode::getRange - Child found -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return it_child->second->getRange(kr_pat_vector, pos, processed + child_len, karp_rabin, kr_factors, pattern, hash_nano);
			}
		}
		else{
			cout << "HashTrieNode::getRange - Case 2, child_len: " << child_len << "\n";
			
//			cout << "Node: \n";
//			it_child->second->print(0);
			// Caso de borde detectado a veces en kr_factors->hashFast(it_child->second->min_factor_pos, processed, pat_len
			
			unsigned long long hash = kr_factors->hashFast(it_child->second->min_factor_pos, processed, pat_len);
			
			
			hash_pat = karp_rabin->subtract_prefix(kr_pat_vector[kr_pat_vector.size() - 1], kr_pat_vector[pos + processed - 1], kr_pat_vector.size() - pos - processed);
//			cout << "HashTrieNode::getRange - hash: " << hash << ", hash_pat: " << hash_pat << "\n";
			if( hash == hash_pat ){
				cout << "HashTrieNode::getRange - Child found -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return pair<unsigned int, unsigned int>(it_child->second->min, it_child->second->max);
			}
			
		}
	}
	
	cout << "HashTrieNode::getRange - Pattern NOT found\n";
	return pair<unsigned int, unsigned int>((unsigned int)(-1), (unsigned int)(-1));
}

void HashTrieNode::save(fstream &writer){
	
	writer.write((char*)&len, sizeof(int));
	writer.write((char*)&min, sizeof(int));
	writer.write((char*)&max, sizeof(int));
	writer.write((char*)&hash, sizeof(int));
	writer.write((char*)&min_factor_pos, sizeof(int));
	
	unsigned int n_childs = childs.size();
	writer.write((char*)&n_childs, sizeof(int));
	for( auto it : childs ){
		char child_first_char = it.first;
		writer.write((char*)&child_first_char, 1);
		it.second->save(writer);
	}
	
}

void HashTrieNode::load(fstream &reader, unsigned int processed, vector<unsigned int> &factors_start, const char *full_text){
	
	reader.read((char*)&len, sizeof(int));
	reader.read((char*)&min, sizeof(int));
	reader.read((char*)&max, sizeof(int));
	reader.read((char*)&hash, sizeof(int));
	reader.read((char*)&min_factor_pos, sizeof(int));
	
	unsigned int n_childs = 0;
	reader.read((char*)&n_childs, sizeof(int));
	for(unsigned int i = 0; i < n_childs; ++i){
		char child_first_char = 0;
		reader.read((char*)&child_first_char, 1);
		childs[child_first_char] = std::make_shared<HashTrieNode>();
		childs[child_first_char]->load(reader, processed, factors_start, full_text);
	}

}

void HashTrie::save(const string &file){
	cout << "HashTrie::save - Start (" << file << ")\n";
	fstream writer(file, fstream::out | fstream::trunc);
	root.save(writer);
	writer.close();
	cout << "HashTrie::save - End\n";
}

void HashTrie::load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, const string &file, vector<unsigned int> &factors_start, const char *full_text){
	cout << "HashTrie::load - Start (" << file << ")\n";
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	fstream reader(file, fstream::in);
	root.load(reader, 0, factors_start, full_text);
	reader.close();
	cout << "HashTrie::load - End\n";
}
	
// VERSIONES REV

HashTrieRev::HashTrieRev(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

HashTrieRev::HashTrieRev(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	build(full_text, len_text, factors_start, arr_x, _karp_rabin, _kr_factors);
}

HashTrieRev::~HashTrieRev(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

void HashTrieRev::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	
	cout << "HashTrieRev::build - Start (full text of " << len_text << ", " << factors_start.size() << " factors)\n";
	
	root.build(full_text, len_text, factors_start, arr_x, karp_rabin, kr_factors, 0, factors_start.size()-1, 0);
	
	cout << "HashTrieRev::build - End\n";
}

HashTrieRevNode::HashTrieRevNode(){
	len = 0;
	min = 0;
	max = 0;
	hash = 0;
	min_factor_pos = 0;
}

HashTrieRevNode::~HashTrieRevNode(){

}

void HashTrieRevNode::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_x, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len){
	if( min == max ){
//		cout << "HashTrieRevNode::build - Leaf\n";
		return;
	}
	if( full_text == NULL || processed_len == len_text ){
//		cout << "HashTrieRevNode::build - Null text\n";
		return;
	}
	
//	cout << "HashTrieRevNode::build - Start (full text of " << len_text << ", range[" << min << ", " << max << "], processed_len: " << processed_len << ")\n";
	
	unsigned int min_pos = min;
	
	unsigned int min_text_len = 0;
	if( arr_x[min_pos] > 0 ){
		min_text_len = factors_start[ arr_x[min_pos] ] - factors_start[ arr_x[min_pos] - 1 ] - processed_len;
	}
	while( min_text_len == 0 && min_pos < max ){
		++min_pos;
		if( arr_x[min_pos] > 0 ){
			min_text_len = factors_start[ arr_x[min_pos] ] - factors_start[ arr_x[min_pos] - 1 ] - processed_len;
		}
	}
	
	unsigned int min_text_start = 0;
	if( arr_x[min_pos] > 0 ){
		min_text_start = factors_start[ arr_x[min_pos] ] - processed_len - 1;
	}
	
	unsigned int cur_pos = 0;
	unsigned int cur_text_start = 0;
	unsigned int cur_text_len = 0;
	unsigned int min_common_text = 0;
	
	unsigned long long hash;
	
	while( min_pos <= max ){
		
		// Primero busco el cur_pos MAYOR que tenga *primer* char igual al de min_pos
		// Despues reviso el largo mayor
		
//		cout << "HashTrieRevNode::build - Iniciando busqueda binaria en [" << min_pos << ", " << max << "]\n";
		
		unsigned int l = min_pos;
		unsigned int h = max;
		while(l < h){
			cur_pos = l + ((h-l)>>1);
//			cur_text_start = factors_start[ arr_y[cur_pos] ] + processed_len;
			cur_text_start = factors_start[ arr_x[cur_pos] ] - processed_len - 1;
			if( arr_x[cur_pos] == 0 ){
				cur_text_len = 0;
			}
			else{
				cur_text_len = factors_start[ arr_x[cur_pos] ] - factors_start[ arr_x[cur_pos] - 1 ] - processed_len;
			}
			if( full_text[cur_text_start] <= full_text[min_text_start] ){
				l = cur_pos+1;
			}
			else{
				h = cur_pos;
			}
		}
		cur_pos = h;
//		cur_text_start = factors_start[ arr_y[cur_pos] ] + processed_len;
		cur_text_start = factors_start[ arr_x[cur_pos] ] - processed_len - 1;
		if( arr_x[cur_pos] == 0 ){
			cur_text_len = 0;
		}
		else{
			cur_text_len = factors_start[ arr_x[cur_pos] ] - factors_start[ arr_x[cur_pos] - 1 ] - processed_len;
		}
		if( (cur_pos > min_pos) && full_text[cur_text_start] != full_text[min_text_start] ){
			--cur_pos;
//			cur_text_start = factors_start[ arr_y[cur_pos] ] + processed_len;
			cur_text_start = factors_start[ arr_x[cur_pos] ] - processed_len - 1;
			if( arr_x[cur_pos] == 0 ){
				cur_text_len = 0;
			}
			else{
				cur_text_len = factors_start[ arr_x[cur_pos] ] - factors_start[ arr_x[cur_pos] - 1 ] - processed_len;
			}
		}
//		cout << "HashTrieRevNode::build - cur_pos: " << cur_pos << "\n";
		
		// Si cur_pos == min_pos, el largo es min_text_len
		// Si no, buscar el largo maximo entre min y cur
		
		if( cur_pos > min_pos ){
//			cout << "HashTrieRevNode::build - Case 1\n";
			// Buscar el mayor
			// Aqui uso el hecho de que full_text termina en '\0'
			min_common_text = 0;
			while( min_common_text < ((min_text_len<cur_text_len)?min_text_len:cur_text_len) 
				&& full_text[min_text_start - min_common_text] == full_text[cur_text_start - min_common_text] 
				){
				++min_common_text;
			}
		}
		else{
//			cout << "HashTrieRevNode::build - Case 2\n";
			min_common_text = min_text_len;
		}
//		cout << "HashTrieRevNode::build - min_common_text: " << min_common_text << "\n";
//		cout << "HashTrieRevNode::build - Preparing string s (min_text_start: " << min_text_start << " / " << len_text << ", len: " << min_common_text << ")\n";
		
		// En el caso de Rev, uso s para el hash directo
		string s;
		for( unsigned int i = 0; i < min_common_text; ++i ){
			s += *(full_text + min_text_start - i);
		}
		hash = karp_rabin->hash(s);
//		cout << "HashTrieRevNode::build - s: " << s << "\n";
		char first_char = s[0];
//		cout << "HashTrieRevNode::build - hash: " << hash << "\n";
		
//		cout << "HashTrieRevNode::build - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_common_text << ", hash: " << hash << "\n";
		// Preparar el hijo, ejecutar la llamda sobre esa instancia
		if( min_common_text == 0 ){
//			cout << "HashTrieRevNode::build - Omiting child of len 0\n";
		}
		else{
			childs[first_char] = std::make_shared<HashTrieRevNode>();
			childs[first_char]->len = min_common_text;
			childs[first_char]->min = min_pos;
			childs[first_char]->max = cur_pos;
			childs[first_char]->text = s;
			childs[first_char]->hash = hash;
			childs[first_char]->min_factor_pos = arr_x[min_pos];
			childs[first_char]->build(full_text, len_text, factors_start, arr_x, karp_rabin, kr_factors, min_pos, cur_pos, processed_len + min_common_text);
		}
		
//		cout << "HashTrieRevNode::build - Preparing min_pos = " << cur_pos+1 << " / " << arr_x.size() << "\n";
		min_pos = cur_pos+1;
		if(min_pos >= arr_x.size()){
			break;
		}
		min_text_start = 0;
		min_text_len = 0;
		if( arr_x[min_pos] > 0 ){
			min_text_start = factors_start[ arr_x[min_pos] ] - processed_len - 1;
			min_text_len = factors_start[ arr_x[min_pos] ] - factors_start[ arr_x[min_pos] - 1 ] - processed_len;
		}
//		cout << "HashTrieRevNode::build - min_text_start: " << min_text_start << " de len " << min_text_len << "\n";
	
	}
	
//	cout << "----- \n";
	
}


void HashTrieRev::print(){
	root.print(0);
	cout << "Root: " << root.childs.size() << " childs\n";
	cout << "Lengths: ";
	for( auto it_child : root.childs ){
		cout << it_child.second->len << " - ";
	}
	cout << "\n";
}

void HashTrieRevNode::print(unsigned int level){
	for(unsigned int i = 0; i < level; ++i){
		cout << "- ";
	}
//	cout << "hash: " << hash << ", len: " << len << ", factor_base " << factor << " , range [" << min << ", " << max << "]\n";
	cout << "hash: " << hash << ", len: " << len << " , range [" << min << ", " << max << "]\n";
	for( auto it : childs ){
		cout << it.first << " ";
		it.second->print(level+1);
	}
}

unsigned int HashTrieRevNode::totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height){
	unsigned int ret = childs.size();
	if( height > max_height ){
		max_height = height;
	}
	if( len > max_len ){
		max_len = len;
	}
	if( childs.size() > max_childs ){
		max_childs = childs.size();
	}
	for( auto it : childs ){
		ret += it.second->totalChilds(max_len, max_childs, max_height, 1 + height);
	}
	return ret;
}

void HashTrieRev::printSize(){
	unsigned int max_len = 0;
	unsigned int max_childs = 0;
	unsigned int max_height = 0;
	unsigned int total_childs = root.totalChilds(max_len, max_childs, max_height, 0);
	cout << "HashTrieRev::printSize - totalChilds: " << total_childs << ", max_len: " << max_len << ", max_childs: " << max_childs << ", max_height: " << max_height << "\n";
}

pair<unsigned int, unsigned int> HashTrieRev::getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, const string &pattern_rev){
	return root.getRange(kr_pat_rev_vector, pos, 0, karp_rabin, kr_factors, pattern_rev);
}

pair<unsigned int, unsigned int> HashTrieRevNode::getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, const string &pattern_rev){
	
	cout << "HashTrieRevNode::getRange - Start (prefixes: " << kr_pat_rev_vector.size() << ", pos: " << pos << ", processed: " << processed << ", pattern_rev: " << pattern_rev << ")\n";
	
	if( processed >= pos ){
		cout << "HashTrieRevNode::getRange - [" << min << ", " << max << "]\n";
		return pair<unsigned int, unsigned int>(min, max);
	}
	
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	unsigned int pat_len = pos - processed;
	char first_char_pat = pattern_rev[kr_pat_rev_vector.size() - 1 - pos + processed];
//	cout << "HashTrieRevNode::getRange - pat = pattern_rev.substr(" << (kr_pat_rev_vector.size() - 1 - pos + processed) << ", " << pat_len << ");\n";
	string pat = pattern_rev.substr(kr_pat_rev_vector.size() - 1 - pos + processed, pat_len);
	cout << "HashTrieRevNode::getRange - pat: " << pat << ", first_char_pat: " << first_char_pat << "\n";
	
	auto it_child = childs.find(first_char_pat);
	if(it_child != childs.end()){
		child_len = it_child->second->len;
		if( child_len <= pat_len ){
//			cout << "HashTrieRevNode::getRange - Case 1, child_len: " << child_len << "\n";
			hash_pat = karp_rabin->hash(pattern_rev.c_str() + pattern_rev.length() - pos + processed, child_len);
			string pat_cut = pattern_rev.substr(kr_pat_rev_vector.size() - 1 - pos + processed, child_len);
//			cout << "HashTrieRevNode::getRange - pat_cut: " << pat_cut << " hash_pat: " << hash_pat << " / " << karp_rabin->hash(pat_cut) << " (processed: " << processed << ")\n";
			if( hash_pat == it_child->second->hash ){
//				cout << "HashTrieRevNode::getRange - Child found -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return it_child->second->getRange(kr_pat_rev_vector, pos, processed + child_len, karp_rabin, kr_factors, pattern_rev);
			}
		}
		else{
			cout << "HashTrieRevNode::getRange - Case 2, child_len: " << child_len << "\n";
			unsigned long long hash = karp_rabin->hash(it_child->second->text.c_str(), pat_len);
			hash_pat = karp_rabin->hash(pattern_rev.c_str() + pattern_rev.length() - pos + processed, pat_len);
//			cout << "HashTrieRevNode::getRange - hash: " << hash << ", hash_pat: " << hash_pat << "\n";
			if( hash == hash_pat ){
				cout << "HashTrieRevNode::getRange - Child found -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return pair<unsigned int, unsigned int>(it_child->second->min, it_child->second->max);
			}
			
		}
	}
	
	cout << "HashTrieRevNode::getRange - Pattern NOT found\n";
	return pair<unsigned int, unsigned int>((unsigned int)(-1), (unsigned int)(-1));
}

void HashTrieRevNode::save(fstream &writer){
	
	writer.write((char*)&len, sizeof(int));
	writer.write((char*)&min, sizeof(int));
	writer.write((char*)&max, sizeof(int));
	writer.write((char*)&hash, sizeof(int));
	writer.write((char*)&min_factor_pos, sizeof(int));
	unsigned int text_len = text.length();
	writer.write((char*)&text_len, sizeof(int));
	writer.write((char*)(text.c_str()), text_len);
	
	unsigned int n_childs = childs.size();
	writer.write((char*)&n_childs, sizeof(int));
	for( auto it : childs ){
		char child_first_char = it.first;
		writer.write((char*)&child_first_char, 1);
		it.second->save(writer);
	}
	
}

void HashTrieRevNode::load(fstream &reader){
	
	reader.read((char*)&len, sizeof(int));
	reader.read((char*)&min, sizeof(int));
	reader.read((char*)&max, sizeof(int));
	reader.read((char*)&hash, sizeof(int));
	reader.read((char*)&min_factor_pos, sizeof(int));
	unsigned int text_len = 0;
	reader.read((char*)&text_len, sizeof(int));
	char buff[text_len + 1];
	reader.read(buff, text_len);
	buff[text_len] = 0;
	text = string(buff);
	
	unsigned int n_childs = 0;
	reader.read((char*)&n_childs, sizeof(int));
	for(unsigned int i = 0; i < n_childs; ++i){
		char child_first_char = 0;
		reader.read((char*)&child_first_char, 1);
		childs[child_first_char] = std::make_shared<HashTrieRevNode>();
		childs[child_first_char]->load(reader);
	}
	
}

void HashTrieRev::save(const string &file){
	cout << "HashTrieRev::save - Start (" << file << ")\n";
	fstream writer(file, fstream::out | fstream::trunc);
	root.save(writer);
	writer.close();
	cout << "HashTrieRev::save - End\n";
}
	
void HashTrieRev::load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, const string &file){
	cout << "HashTrieRev::load - Start (" << file << ")\n";
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	fstream reader(file, fstream::in);
	root.load(reader);
	reader.close();
	cout << "HashTrieRev::load - End\n";
}











