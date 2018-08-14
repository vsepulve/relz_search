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
	min_factor_pos = 0;
}

HashTrieNode::~HashTrieNode(){

}

unsigned int HashTrieNode::getMaxComp(const char *str_1, unsigned int len_1, const char *str_2, unsigned int len_2){
//	cout << "HashTrieNode::getMaxComp - Start\n";
	unsigned int ret = 0;
	for( ; ret < ( (len_1<len_2)?len_1:len_2 ); ++ret ){
		if( str_1[ret] != str_2[ret] ){
			break;
		}
	}
//	cout << "HashTrieNode::getMaxComp - End (" << ret << " / " << ((len_1<len_2)?len_1:len_2) << ")\n";
	return ret;
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
		
		if( max - min_pos < 0 ){
//			cout << "HashTrieNode::build - Busqueda simple\n";
			cur_pos = min_pos;
			cur_text_start = factors_start[ arr_y[cur_pos] ] + processed_len;
			
//			cout << "HashTrieNode::build - cur_pos: " << cur_pos<< " / " << max << ", " << full_text[cur_text_start] << " vs " << full_text[min_text_start] << "\n";
			while( cur_pos <= max && full_text[cur_text_start] == full_text[min_text_start] ){
				++cur_pos;
				cur_text_start = factors_start[ arr_y[cur_pos] ] + processed_len;
//				cout << "HashTrieNode::build - cur_pos: " << cur_pos<< " / " << max << ", " << full_text[cur_text_start] << " vs " << full_text[min_text_start] << "\n";
			}
			if( (cur_pos > min_pos) && full_text[cur_text_start] != full_text[min_text_start] ){
				--cur_pos;
				cur_text_start = factors_start[ arr_y[cur_pos] ] + processed_len;
			}
		}
		else{
//			cout << "HashTrieNode::build - Busqueda binaria\n";
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
		}
//		cout << "HashTrieNode::build - cur_pos: " << cur_pos << "\n";
		
		// Si cur_pos == min_pos, el largo es min_text_len
		// Si no, buscar el largo maximo entre min y cur
		
		if( cur_pos > min_pos ){
			// Buscar el mayor
			// Aqui uso el hecho de que full_text termina en '\0'
			min_text_len = 0;
			while( full_text[min_text_start + min_text_len] == full_text[cur_text_start + min_text_len] ){
				++min_text_len;
			}
		}
		else{
			min_text_len = len_text - min_text_start;
		}
//		cout << "HashTrieNode::build - min_text_len: " << min_text_len << "\n";
//		cout << "HashTrieNode::build - Preparing string s (min_text_start: " << min_text_start << " / " << len_text << ", len: " << min_text_len << ")\n";
		
		// El string es olo para el mensaje de debug
//		string s(full_text + min_text_start, min_text_len);
		
		// Notar que aqui necesito la posicion del factor en la coleccion, no en el arreglo Y
		hash = kr_factors->hash(arr_y[min_pos], processed_len, min_text_len);
//		cout << "HashTrieNode::build - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_text_len << ", hash: " << hash << " / " << karp_rabin->hash(s) << " (\"" << s << "\")\n";
//		cout << "HashTrieNode::build - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_text_len << ", hash: " << hash << "\n";
		// Preparar el hijo, ejecutar la llamda sobre esa instancia
		if( min_text_len == 0 ){
//			cout << "HashTrieNode::build - Omiting child of len 0\n";
		}
		else{
			childs[hash] = std::make_shared<HashTrieNode>();
			childs[hash]->len = min_text_len;
			childs[hash]->min = min_pos;
			childs[hash]->max = cur_pos;
			childs[hash]->min_factor_pos = arr_y[min_pos];
			childs[hash]->build(full_text, len_text, factors_start, arr_y, karp_rabin, kr_factors, min_pos, cur_pos, processed_len + min_text_len);
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
//	for( auto it_child : root.childs ){
//		cout << "child->childs: " << it_child.second->childs.size() << " (";
//		for( auto it : it_child.second->childs ){
//			cout << it.second->len << " - ";
//		}
//		cout << ")\n";
//	}
}

void HashTrieNode::print(unsigned int level){
	for(unsigned int i = 0; i < level; ++i){
		cout << "- ";
	}
//	cout << "hash: " << hash << ", len: " << len << ", factor_base " << factor << " , range [" << min << ", " << max << "]\n";
	cout << "len: " << len << " , range [" << min << ", " << max << "]\n";
	for( auto it : childs ){
		it.second->print(level+1);
	}
}

pair<unsigned int, unsigned int> HashTrie::getRange(const string &pattern){
	return root.getRange(pattern.c_str(), pattern.length(), 0, karp_rabin, kr_factors);
}

pair<unsigned int, unsigned int> HashTrie::getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, const string &pattern){
	return root.getRange(kr_pat_vector, pos, 0, karp_rabin, kr_factors, pattern, &hash_nano);
}

pair<unsigned int, unsigned int> HashTrieNode::getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, const string &pattern, unsigned long long *hash_nano){
	
//	cout << "HashTrieNode::getRange - Start (prefixes: " << kr_pat_vector.size() << ", pos: " << pos << ", processed: " << processed << ")\n";
	
	if( pos + processed >= kr_pat_vector.size() ){
//		cout << "HashTrieNode::getRange - [" << min << ", " << max << "]\n";
		return pair<unsigned int, unsigned int>(min, max);
	}
	
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	unsigned int pat_len = kr_pat_vector.size() - pos - processed;
//	cout << "HashTrieNode::getRange - pat_len: " << pat_len << "\n";
//	string pat = pattern.substr(pos + processed, pat_len);
//	cout << "HashTrieNode::getRange - pat: " << pat << "\n";
	
	// Version childs length
	for( unsigned int i = 0; i < childs_pairs.size(); ++i ){
		child_len = childs_pairs[i].first;
		if( child_len <= pat_len ){
//			cout << "HashTrieNode::getRange - Quick Search for child_len " << child_len << "\n";
//			cout << "HashTrieNode::getRange - Restando hash [" << (pos + processed + child_len - 1) <<"] = " << kr_pat_vector[pos + processed + child_len - 1] << " - hash[" << (pos + processed - 1) << "] = " << kr_pat_vector[pos + processed - 1] << "\n";
			hash_pat = karp_rabin->subtract_prefix(kr_pat_vector[pos + processed + child_len - 1], kr_pat_vector[pos + processed - 1], child_len);
//			string pat_cut = pattern.substr(pos + processed, child_len);
//			cout << "HashTrieNode::getRange - pat_cut: " << pat_cut << " hash_pat: " << hash_pat << " / " << karp_rabin->hash(pat_cut) << " (processed: " << processed << ")\n";
			auto it_child = childs.find(hash_pat);
			if( it_child != childs.end() ){
//				cout << "HashTrieNode::getRange - Child found -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return it_child->second->getRange(kr_pat_vector, pos, processed + child_len, karp_rabin, kr_factors, pattern, hash_nano);
			}
		}
		else{
			vector<shared_ptr<HashTrieNode>> childs_vector = childs_pairs[i].second;
//			cout << "HashTrieNode::getRange - Slow search for child_len " << child_len << " (" << childs_vector.size() << " childs)\n";
			for(unsigned int j = 0; j < childs_vector.size(); ++j){
				// Si este nodo contiene el texto que queda del patron, este es el rango
				
				NanoTimer timer;
//				unsigned long long hash = kr_factors->hash(childs_vector[j]->min_factor_pos, processed, pat_len);
				unsigned long long hash = kr_factors->hashFast(childs_vector[j]->min_factor_pos, processed, pat_len);
				(*hash_nano) += timer.getNanosec();
				
				hash_pat = karp_rabin->subtract_prefix(kr_pat_vector[kr_pat_vector.size() - 1], kr_pat_vector[pos + processed - 1], kr_pat_vector.size() - pos - processed);
//				cout << "HashTrieNode::getRange - hash: " << hash << ", hash_pat: " << hash_pat << "\n";
				if( hash == hash_pat ){
//					cout << "HashTrieNode::getRange - Child found -> [" << childs_vector[j]->min << ", " << childs_vector[j]->max << "]\n";
					return pair<unsigned int, unsigned int>(childs_vector[j]->min, childs_vector[j]->max);
				}
				
			}
		}
	}
	
	
//	cout << "HashTrieNode::getRange - Patron NO encontrado\n";
	return pair<unsigned int, unsigned int>(0, 0);
}

pair<unsigned int, unsigned int> HashTrieNode::getRange(const char *pattern, unsigned int pat_len, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors){
	
	cout << "HashTrieNode::getRange - Start (\"" << pattern << "\", " << pat_len << ", processed: " << processed << ")\n";
	
	if( pat_len == 0 ){
		cout << "HashTrieNode::getRange - [" << min << ", " << max << "]\n";
		return pair<unsigned int, unsigned int>(min, max);
	}
	
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	unsigned long long hash_pat_full = 0;
	bool computed_full = false;
	
	// Version childs length
	for( unsigned int i = 0; i < childs_pairs.size(); ++i ){
		child_len = childs_pairs[i].first;
		if( child_len <= pat_len ){
			cout << "HashTrieNode::getRange - Quick Search for child_len " << child_len << "\n";
			hash_pat = karp_rabin->hash(pattern, child_len);
			auto it_child = childs.find(hash_pat);
			if( it_child != childs.end() ){
				cout << "HashTrieNode::getRange - Child found -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return it_child->second->getRange(pattern + child_len, pat_len - child_len, processed + child_len, karp_rabin, kr_factors);
			}
		}
		else{
			vector<shared_ptr<HashTrieNode>> childs_vector = childs_pairs[i].second;
			cout << "HashTrieNode::getRange - Slow search for child_len " << child_len << " (" << childs_vector.size() << " childs)\n";
			for(unsigned int j = 0; j < childs_vector.size(); ++j){
				// Si este nodo contiene el texto que queda del patron, este es el rango
				unsigned long long hash = kr_factors->hash(childs_vector[j]->min_factor_pos, processed, pat_len);
				if( computed_full ){
					hash_pat = hash_pat_full;
				}
				else{
					hash_pat = hash_pat_full = karp_rabin->hash(pattern, pat_len);
					computed_full = true;
				}
				if( hash == hash_pat ){
					cout << "HashTrieNode::getRange - Child found -> [" << childs_vector[j]->min << ", " << childs_vector[j]->max << "]\n";
					return pair<unsigned int, unsigned int>(childs_vector[j]->min, childs_vector[j]->max);
				}
				
			}
		}
	}
	
	cout << "HashTrieNode::getRange - Patron NO encontrado\n";
	return pair<unsigned int, unsigned int>(0, 0);
}

void HashTrieNode::save(fstream &writer){
	
	writer.write((char*)&len, sizeof(int));
	writer.write((char*)&min, sizeof(int));
	writer.write((char*)&max, sizeof(int));
	writer.write((char*)&min_factor_pos, sizeof(int));
	
	unsigned int n_childs = childs.size();
	writer.write((char*)&n_childs, sizeof(int));
	for( auto it : childs ){
		unsigned int hash = (unsigned int)(it.first);
		writer.write((char*)&hash, sizeof(int));
		it.second->save(writer);
	}
	
}

void HashTrieNode::load(fstream &reader){
	
	reader.read((char*)&len, sizeof(int));
	reader.read((char*)&min, sizeof(int));
	reader.read((char*)&max, sizeof(int));
	reader.read((char*)&min_factor_pos, sizeof(int));
	
	unsigned int n_childs = 0;
	reader.read((char*)&n_childs, sizeof(int));
	for(unsigned int i = 0; i < n_childs; ++i){
		unsigned int hash = 0;
		reader.read((char*)&hash, sizeof(int));
		childs[hash] = std::make_shared<HashTrieNode>();
		childs[hash]->load(reader);
	}
	
}

void HashTrie::save(const string &file){
	cout << "HashTrie::save - Start (" << file << ")\n";
	fstream writer(file, fstream::out | fstream::trunc);
	root.save(writer);
	writer.close();
	cout << "HashTrie::save - End\n";
}

void HashTrie::load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, const string &file){
	cout << "HashTrie::load - Start (" << file << ")\n";
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	fstream reader(file, fstream::in);
	root.load(reader);
	reader.close();
	cout << "HashTrie::load - End\n";
}

void HashTrieNode::prepareChilds(){
//	cout << "HashTrieNode::prepareChilds - Start\n";
//	root.prepareChilds();
//	unordered_map<unsigned int, shared_ptr<HashTrieNode>> childs;
//	vector< pair<unsigned int, vector<shared_ptr<HashTrieNode>>> > childs_pairs;
	childs_pairs.clear();
	map< unsigned int, vector<shared_ptr<HashTrieNode>> > childs_map;
	for( auto it_child : childs ){
		shared_ptr<HashTrieNode> ptr = it_child.second;
		unsigned int len = ptr->len;
		auto it_map = childs_map.find(len);
		if( it_map == childs_map.end() ){
			vector<shared_ptr<HashTrieNode>> childs_vector;
			childs_vector.push_back(ptr);
			childs_map[len] = childs_vector;
		}
		else{
			it_map->second.push_back(ptr);
		}
	}
	for( auto it_map : childs_map ){
		unsigned int len = it_map.first;
		vector<shared_ptr<HashTrieNode>> childs_vector = it_map.second;
		childs_pairs.push_back(
			pair<unsigned int, vector<shared_ptr<HashTrieNode>>>(
				len, childs_vector
			)
		);
	}
//	cout << "HashTrieNode::prepareChilds - Sorting pairs\n";
	sort(childs_pairs.begin(), childs_pairs.end());
	
//	cout << "HashTrieNode::prepareChilds - Resulting pairs:\n";
//	for( auto it_pair : childs_pairs ){
//		cout << "HashTrieNode::prepareChilds - (len: " << it_pair.first << ", n_childs: " << it_pair.second.size() << ")\n";
//	}
	
//	cout << "HashTrieNode::prepareChilds - Continuing with " << childs.size() << " childs\n";
	for( auto it_child : childs ){
		it_child.second->prepareChilds();
	}
	
//	cout << "HashTrieNode::prepareChilds - End\n";
}

void HashTrie::prepareChilds(){
	cout << "HashTrie::prepareChilds - Start\n";
	root.prepareChilds();
	cout << "HashTrie::prepareChilds - End\n";
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
	min_factor_pos = 0;
}

HashTrieRevNode::~HashTrieRevNode(){

}

unsigned int HashTrieRevNode::getMaxComp(const char *str_1, unsigned int len_1, const char *str_2, unsigned int len_2){
//	cout << "HashTrieRevNode::getMaxComp - Start\n";
	unsigned int ret = 0;
	for( ; ret < ( (len_1<len_2)?len_1:len_2 ); ++ret ){
//		cout << "HashTrieRevNode::getMaxComp - " << *(str_1 - ret) << " vs " << *(str_2 - ret) << "\n";
		if( *(str_1 - ret) != *(str_2 - ret) ){
			break;
		}
	}
//	cout << "HashTrieRevNode::getMaxComp - End (" << ret << " / " << ((len_1<len_2)?len_1:len_2) << ")\n";
	return ret;
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
	unsigned int min_text_start = 0;
	unsigned int min_text_len = 0;
	if( arr_x[min_pos] > 0 ){
		min_text_start = factors_start[ arr_x[min_pos] ] - processed_len - 1;
		min_text_len = factors_start[ arr_x[min_pos] ] - factors_start[ arr_x[min_pos] - 1 ] - processed_len;
	}
//	cout << "HashTrieRevNode::build - min_text_len: " << min_text_len << " (f " << arr_x[min_pos] << ")\n";
	
	unsigned long long hash;
	
	while( min_pos <= max ){
		
//		cout << "HashTrieRevNode::build - Preparing cur_pos = " << (min_pos + 1) << "\n";
		unsigned int cur_pos = min_pos + 1;
		unsigned int cur_text_start = 0;
		unsigned int cur_text_len = 0;
		if( arr_x[cur_pos] > 0 ){
			cur_text_start = factors_start[ arr_x[cur_pos] ] - processed_len - 1;
			cur_text_len = factors_start[ arr_x[cur_pos] ] - factors_start[ arr_x[cur_pos] - 1 ] - processed_len;
		}
//		cout << "HashTrieRevNode::build - cur_text_len: " << cur_text_len << " (f " << arr_x[cur_pos] << ")\n";
		
		unsigned int max_comp = 0;
		if( cur_pos < factors_start.size() && cur_pos <= max ){
//			cout << "HashTrieRevNode::build - min_text_start: " << min_text_start << ", cur_text_start: " << cur_text_start << "\n";
			max_comp = getMaxComp(full_text + min_text_start, min_text_len, full_text + cur_text_start, cur_text_len);
		}
		
		while( max_comp > 0 && cur_pos <= max ){
			min_text_len = max_comp;
			++cur_pos;
			if( cur_pos == factors_start.size() ){
				max_comp = 0;
				break;
			}
			cur_text_start = 0;
			cur_text_len = 0;
			if( arr_x[cur_pos] > 0 ){
				cur_text_start = factors_start[ arr_x[cur_pos] ] - processed_len - 1;
				cur_text_len = factors_start[ arr_x[cur_pos] ] - factors_start[ arr_x[cur_pos] - 1 ] - processed_len;
			}
//			cout << "HashTrieRevNode::build - cur_text_len: " << cur_text_len << " (f " << arr_x[cur_pos] << ")\n";
//			cout << "HashTrieRevNode::build - min_text_start: " << min_text_start << ", cur_text_start: " << cur_text_start << "\n";
			max_comp = getMaxComp(full_text + min_text_start, min_text_len, full_text + cur_text_start, cur_text_len);
			
		}
		
//		cout << "HashTrieRevNode::build - Preparing string s (min_text_start: " << min_text_start << " / " << len_text << ", len: " << min_text_len << ")\n";
		
		// En el caso de Rev, uso s para el hash directo
		string s;
		for( unsigned int i = 0; i < min_text_len; ++i ){
			s += *(full_text + min_text_start - i);
		}
		hash = karp_rabin->hash(s);
//		cout << "HashTrieRevNode::build - s: " << s << " (hash: " << hash << ")\n";
		
		// Notar que aqui necesito la posicion del factor en la coleccion, no en el arreglo Y
//		hash = kr_factors->hash(arr_x[min_pos], processed_len, min_text_len);
//		cout << "HashTrieRevNode::build - Adding range [" << min_pos << ", " << cur_pos << "), len " << min_text_len << ", hash: " << hash << " / " << karp_rabin->hash(s) << "\n";
//		cout << "HashTrieRevNode::build - Adding range [" << min_pos << ", " << cur_pos << "), len " << min_text_len << ", hash: " << hash << "\n";
		// Preparar el hijo, ejecutar la llamda sobre esa instancia
		if( min_text_len == 0 ){
//			cout << "HashTrieNode::build - Omiting child of len 0\n";
		}
		else{
			childs[hash] = std::make_shared<HashTrieRevNode>();
			childs[hash]->len = min_text_len;
			childs[hash]->min = min_pos;
			childs[hash]->max = cur_pos-1;
			childs[hash]->min_factor_pos = arr_x[min_pos];
			childs[hash]->text = s;
			childs[hash]->build(full_text, len_text, factors_start, arr_x, karp_rabin, kr_factors, min_pos, cur_pos-1, processed_len + min_text_len);
		}
		
		min_pos = cur_pos;
		if(min_pos >= arr_x.size()){
			break;
		}
		min_text_start = 0;
		min_text_len = 0;
		if( arr_x[min_pos] > 0 ){
			min_text_start = factors_start[ arr_x[min_pos] ] - processed_len - 1;
			min_text_len = factors_start[ arr_x[min_pos] ] - factors_start[ arr_x[min_pos] - 1 ] - processed_len;
		}
//		cout << "HashTrieRevNode::build - min_text_len: " << min_text_len << " (f " << arr_x[min_pos] << ")\n";
	
	}
	
	
}

void HashTrieRev::print(){
	root.print(0);
	cout << "Root: " << root.childs.size() << " childs\n";
	cout << "Lengths: ";
	for( auto it_child : root.childs ){
		cout << it_child.second->len << " - ";
	}
	cout << "\n";
//	for( auto it_child : root.childs ){
//		cout << "child->childs: " << it_child.second->childs.size() << " (";
//		for( auto it : it_child.second->childs ){
//			cout << it.second->len << " - ";
//		}
//		cout << ")\n";
//	}
}

void HashTrieRevNode::print(unsigned int level){
	for(unsigned int i = 0; i < level; ++i){
		cout << "- ";
	}
//	cout << "hash: " << hash << ", len: " << len << ", factor_base " << factor << " , range [" << min << ", " << max << "]\n";
	cout << "len: " << len << " , range [" << min << ", " << max << "]\n";
	for( auto it : childs ){
		it.second->print(level+1);
	}
}

pair<unsigned int, unsigned int> HashTrieRev::getRange(const string &pattern){
	return root.getRange(pattern.c_str(), pattern.length(), 0, karp_rabin, kr_factors);
}

pair<unsigned int, unsigned int> HashTrieRev::getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, const string &pattern_rev){
	return root.getRange(kr_pat_rev_vector, pos, 0, karp_rabin, kr_factors, pattern_rev);
}

pair<unsigned int, unsigned int> HashTrieRevNode::getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, const string &pattern_rev){
	
//	cout << "HashTrieRevNode::getRange - Start (prefixes: " << kr_pat_rev_vector.size() << ", pos: " << pos << ", processed: " << processed << ")\n";
	
//	if( pos + processed >= kr_pat_rev_vector.size() ){
	if( processed >= pos ){
//		cout << "HashTrieRevNode::getRange - [" << min << ", " << max << "]\n";
		return pair<unsigned int, unsigned int>(min, max);
	}
	
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	unsigned long long hash_pat_full = 0;
	bool computed_full = false;
//	unsigned int pat_len = kr_pat_rev_vector.size() - pos - processed;
	unsigned int pat_len = pos - processed;
//	cout << "HashTrieRevNode::getRange - pat_len: " << pat_len << "\n";
//	string pat = pattern_rev.substr(kr_pat_rev_vector.size() - 1 - pos + processed, pat_len);
//	cout << "HashTrieRevNode::getRange - pat: " << pat << "\n";
	
	for( unsigned int i = 0; i < childs_pairs.size(); ++i ){
		child_len = childs_pairs[i].first;
		if( child_len <= pat_len ){
//			cout << "HashTrieRevNode::getRange - Quick Search for child_len " << child_len << "\n";
//			hash_pat = karp_rabin->hash(pattern, child_len);

//			cout << "HashTrieRevNode::getRange - Restando hash [" << (pos) <<"] = " << kr_pat_rev_vector[pos] << " - hash[" << (pos - processed - 1) << "] = " << kr_pat_rev_vector[pos - processed - 1] << "\n";
//			hash_pat = karp_rabin->subtract_prefix(kr_pat_rev_vector[pos], kr_pat_rev_vector[pos - processed - 1], child_len);
			hash_pat = karp_rabin->hash(pattern_rev.c_str() + pattern_rev.length() - pos + processed, child_len);
			
//			string pat_cut = pattern_rev.substr(kr_pat_rev_vector.size() - 1 - pos + processed, child_len);
//			cout << "HashTrieRevNode::getRange - pat_cut: " << pat_cut << " hash_pat: " << hash_pat << " / " << karp_rabin->hash(pat_cut) << " (processed: " << processed << ")\n";
			auto it_child = childs.find(hash_pat);
			if( it_child != childs.end() ){
//				cout << "HashTrieRevNode::getRange - Child found -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return it_child->second->getRange(kr_pat_rev_vector, pos, processed + child_len, karp_rabin, kr_factors, pattern_rev);
			}
		}
		else{
			vector<shared_ptr<HashTrieRevNode>> childs_vector = childs_pairs[i].second;
//			cout << "HashTrieRevNode::getRange - Slow search for child_len " << child_len << " (" << childs_vector.size() << " childs)\n";
			for(unsigned int j = 0; j < childs_vector.size(); ++j){
				// Si este nodo contiene el texto que queda del patron, este es el rango
//				unsigned long long hash = kr_factors->hash(childs_vector[j]->min_factor_pos, processed, pat_len);
				unsigned long long hash = karp_rabin->hash(childs_vector[j]->text.c_str(), pat_len);
//				cout << "HashTrieRevNode::getRange - Node text: " << childs_vector[j]->text << "\n";
				if( computed_full ){
					hash_pat = hash_pat_full;
				}
				else{
//					string pat = pattern_rev.substr(kr_pat_rev_vector.size() - 1 - pos + processed, pat_len);
					hash_pat = hash_pat_full = karp_rabin->hash(pattern_rev.c_str() + pattern_rev.length() - pos + processed, pat_len);
//					hash_pat = hash_pat_full = karp_rabin->subtract_prefix(kr_pat_rev_vector[kr_pat_rev_vector.size() - 1], kr_pat_rev_vector[pos + processed - 1], kr_pat_rev_vector.size() - pos - processed);
					computed_full = true;
				}
//				cout << "HashTrieRevNode::getRange - hash: " << hash << ", hash_pat: " << hash_pat << "\n";
				if( hash == hash_pat ){
//					cout << "HashTrieRevNode::getRange - Child found -> [" << childs_vector[j]->min << ", " << childs_vector[j]->max << "]\n";
					return pair<unsigned int, unsigned int>(childs_vector[j]->min, childs_vector[j]->max);
				}
				
			}
		}
	}
	
//	cout << "HashTrieRevNode::getRange - Patron NO encontrado\n";
	return pair<unsigned int, unsigned int>(0, 0);
}

pair<unsigned int, unsigned int> HashTrieRevNode::getRange(const char *pattern, unsigned int pat_len, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors){
	
	cout << "HashTrieRevNode::getRange - Start (\"" << pattern << "\", " << pat_len << ", processed: " << processed << ")\n";
	
	if( pat_len == 0 ){
		cout << "HashTrieRevNode::getRange - [" << min << ", " << max << "]\n";
		return pair<unsigned int, unsigned int>(min, max);
	}
	
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	unsigned long long hash_pat_full = 0;
	bool computed_full = false;
	
	for( unsigned int i = 0; i < childs_pairs.size(); ++i ){
		child_len = childs_pairs[i].first;
		if( child_len <= pat_len ){
			cout << "HashTrieRevNode::getRange - Quick Search for child_len " << child_len << "\n";
			hash_pat = karp_rabin->hash(pattern, child_len);
			auto it_child = childs.find(hash_pat);
			if( it_child != childs.end() ){
				cout << "HashTrieRevNode::getRange - Child found -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return it_child->second->getRange(pattern + child_len, pat_len - child_len, processed + child_len, karp_rabin, kr_factors);
			}
		}
		else{
			vector<shared_ptr<HashTrieRevNode>> childs_vector = childs_pairs[i].second;
			cout << "HashTrieRevNode::getRange - Slow search for child_len " << child_len << " (" << childs_vector.size() << " childs)\n";
			for(unsigned int j = 0; j < childs_vector.size(); ++j){
				// Si este nodo contiene el texto que queda del patron, este es el rango
				// Notar que en el caso reverso almaceno temporalente el texto 
				// La otra opcion es usar datos del factor (tu, pu, lu) y datos de la referencia
//				unsigned long long hash = kr_factors->hash(childs_vector[j]->min_factor_pos, processed, pat_len);
				unsigned long long hash = karp_rabin->hash(childs_vector[j]->text.c_str(), pat_len);
				cout << "HashTrieRevNode::getRange - Node text: " << childs_vector[j]->text << "\n";
				if( computed_full ){
					hash_pat = hash_pat_full;
				}
				else{
					hash_pat = hash_pat_full = karp_rabin->hash(pattern, pat_len);
					computed_full = true;
				}
				cout << "HashTrieRevNode::getRange - hash: " << hash<< ", hash_pat: " << hash_pat << "\n";
				if( hash == hash_pat ){
					cout << "HashTrieRevNode::getRange - Child found -> [" << childs_vector[j]->min << ", " << childs_vector[j]->max << "]\n";
					return pair<unsigned int, unsigned int>(childs_vector[j]->min, childs_vector[j]->max);
				}
				
			}
		}
	}
	
	cout << "HashTrieRevNode::getRange - Patron NO encontrado\n";
	return pair<unsigned int, unsigned int>(0, 0);
}

void HashTrieRevNode::save(fstream &writer){
	
	writer.write((char*)&len, sizeof(int));
	writer.write((char*)&min, sizeof(int));
	writer.write((char*)&max, sizeof(int));
	writer.write((char*)&min_factor_pos, sizeof(int));
	
	unsigned int n_childs = childs.size();
	writer.write((char*)&n_childs, sizeof(int));
	for( auto it : childs ){
		unsigned int hash = it.first;
		writer.write((char*)&hash, sizeof(int));
		it.second->save(writer);
	}
	
}

void HashTrieRevNode::load(fstream &reader){
	
	reader.read((char*)&len, sizeof(int));
	reader.read((char*)&min, sizeof(int));
	reader.read((char*)&max, sizeof(int));
	reader.read((char*)&min_factor_pos, sizeof(int));
	
	unsigned int n_childs = 0;
	reader.read((char*)&n_childs, sizeof(int));
	for(unsigned int i = 0; i < n_childs; ++i){
		unsigned int hash = 0;
		reader.read((char*)&hash, sizeof(int));
		childs[hash] = std::make_shared<HashTrieRevNode>();
		childs[hash]->load(reader);
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

void HashTrieRevNode::prepareChilds(){
//	cout << "HashTrieRevNode::prepareChilds - Start\n";
//	root.prepareChilds();
//	unordered_map<unsigned int, shared_ptr<HashTrieRevNode>> childs;
//	vector< pair<unsigned int, vector<shared_ptr<HashTrieRevNode>>> > childs_pairs;
	childs_pairs.clear();
	map< unsigned int, vector<shared_ptr<HashTrieRevNode>> > childs_map;
	for( auto it_child : childs ){
		shared_ptr<HashTrieRevNode> ptr = it_child.second;
		unsigned int len = ptr->len;
		auto it_map = childs_map.find(len);
		if( it_map == childs_map.end() ){
			vector<shared_ptr<HashTrieRevNode>> childs_vector;
			childs_vector.push_back(ptr);
			childs_map[len] = childs_vector;
		}
		else{
			it_map->second.push_back(ptr);
		}
	}
	for( auto it_map : childs_map ){
		unsigned int len = it_map.first;
		vector<shared_ptr<HashTrieRevNode>> childs_vector = it_map.second;
		childs_pairs.push_back(
			pair<unsigned int, vector<shared_ptr<HashTrieRevNode>>>(
				len, childs_vector
			)
		);
	}
//	cout << "HashTrieRevNode::prepareChilds - Sorting pairs\n";
	sort(childs_pairs.begin(), childs_pairs.end());
	
//	cout << "HashTrieRevNode::prepareChilds - Resulting pairs:\n";
//	for( auto it_pair : childs_pairs ){
//		cout << "HashTrieRevNode::prepareChilds - (len: " << it_pair.first << ", n_childs: " << it_pair.second.size() << ")\n";
//	}
	
//	cout << "HashTrieRevNode::prepareChilds - Continuing with " << childs.size() << " childs\n";
	for( auto it_child : childs ){
		it_child.second->prepareChilds();
	}
	
//	cout << "HashTrieRevNode::prepareChilds - End\n";
}

void HashTrieRev::prepareChilds(){
	cout << "HashTrieRev::prepareChilds - Start\n";
	root.prepareChilds();
	cout << "HashTrieRev::prepareChilds - End\n";
}











