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

HashTrie::HashTrie(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
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
	unsigned int min_text_len = len_text - min_text_start;
	unsigned long long hash;
	
	set<unsigned int> lengths;
	
	while( min_pos <= max ){
		
//		cout << "HashTrieNode::build - Preparing cur_pos = " << (min_pos + 1) << "\n";
		unsigned int cur_pos = min_pos + 1;
		unsigned int cur_text_start = factors_start[ arr_y[cur_pos] ] + processed_len;
		unsigned int cur_text_len = len_text - cur_text_start;
		
		unsigned int max_comp = 0;
		if( cur_pos < factors_start.size() && cur_pos <= max ){
			max_comp = getMaxComp(full_text + min_text_start, min_text_len, full_text + cur_text_start, cur_text_len);
		}
		
		while( max_comp > 0 && cur_pos <= max ){
			min_text_len = max_comp;
			++cur_pos;
			if( cur_pos == factors_start.size() ){
				max_comp = 0;
				break;
			}
			cur_text_start = factors_start[ arr_y[cur_pos] ] + processed_len;
			cur_text_len = len_text - cur_text_start - processed_len;
			max_comp = getMaxComp(full_text + min_text_start, min_text_len, full_text + cur_text_start, cur_text_len);
			
		}
	
//		cout << "HashTrieNode::build - Preparing string s (min_text_start: " << min_text_start << " / " << len_text << ", len: " << min_text_len << ")\n";
		
		// El string es olo para el mensaje de debug
//		string s(full_text + min_text_start, min_text_len);
//		hash = karp_rabin->hash(full_text, min_text_start, min_text_len, factors_start);
		
		// Notar que aqui necesito la posicion del factor en la coleccion, no en el arreglo Y
		hash = kr_factors->hash(arr_y[min_pos], processed_len, min_text_len);
//		cout << "HashTrieNode::build - Adding range [" << min_pos << ", " << cur_pos << "), len " << min_text_len << ", hash: " << hash << " / " << karp_rabin->hash(s) << "\n";
//		cout << "HashTrieNode::build - Adding range [" << min_pos << ", " << cur_pos << "), len " << min_text_len << ", hash: " << hash << "\n";
		// Preparar el hijo, ejecutar la llamda sobre esa instancia
		if( min_text_len == 0 ){
//			cout << "HashTrieNode::build - Omiting child of len 0\n";
		}
		else{
			childs[hash] = std::make_shared<HashTrieNode>();
			childs[hash]->len = min_text_len;
			childs[hash]->min = min_pos;
			childs[hash]->max = cur_pos-1;
			childs[hash]->min_factor_pos = arr_y[min_pos];
			lengths.insert(min_text_len);
			childs[hash]->build(full_text, len_text, factors_start, arr_y, karp_rabin, kr_factors, min_pos, cur_pos-1, processed_len + min_text_len);
		}
		
//		cout << "HashTrieNode::build - Preparing min_pos = " << cur_pos << " / " << arr_y.size() << "\n";
		min_pos = cur_pos;
		if(min_pos >= arr_y.size()){
			break;
		}
		min_text_start = factors_start[ arr_y[min_pos] ] + processed_len;
		min_text_len = len_text - min_text_start;
//		cout << "HashTrieNode::build - min_text_start: " << min_text_start << ", min_text_len: " << min_text_len << "\n";
	
	}
	
//	cout << "HashTrieNode::build - Adding childs " << lengths.size() << " lentghs\n";
	for( unsigned int len : lengths ){
		childs_lenghts.push_back(len);
	}
	
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
//	return pair<unsigned int, unsigned int>(0, 0);
	return root.getRange(pattern.c_str(), pattern.length(), 0, karp_rabin, kr_factors);
}

pair<unsigned int, unsigned int> HashTrieNode::getRange(const char *pattern, unsigned int pat_len, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors){
	
	cout << "HashTrieNode::getRange - Start (\"" << pattern << "\", " << pat_len << ", processed: " << processed << ")\n";
	
	if( pat_len == 0 ){
		return pair<unsigned int, unsigned int>(min, max);
	}
	
	map<unsigned int, unsigned long long> hash_pat_map;
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	
	// Si todos los hijos del nodo actual son menores que el largo del patron, se puede hacer el descenso rapido
	// Los largos de los hijos ESTAN ordenados de menor a mayor
	cout << "HashTrieNode::getRange - Childs lengths: \n";
	for( unsigned int len : childs_lenghts ){
		cout << len << " - ";
	}
	cout << "\n";
	
	if( childs_lenghts.back() <= pat_len ){
		cout << "HashTrieNode::getRange - Quick Search\n";
		for( unsigned int len : childs_lenghts ){
			cout << "HashTrieNode::getRange - Probando len " << len << "\n";
			hash_pat = karp_rabin->hash(pattern, len);
			auto it_child = childs.find(hash_pat);
			if( it_child != childs.end() ){
				child_len = it_child->second->len;
				cout << "HashTrieNode::getRange - Hijo encontrado de len " << child_len << " -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return it_child->second->getRange(pattern + child_len, pat_len - child_len, processed + child_len, karp_rabin, kr_factors);
			}
		}
		cout << "HashTrieNode::getRange - Patron NO encontrado\n";
		return pair<unsigned int, unsigned int>(0, 0);
	}
	
	for( auto it_child : childs ){
		child_len = it_child.second->len;
		cout << "HashTrieNode::getRange - hijo de len " << child_len << " -> [" << it_child.second->min << ", " << it_child.second->max << "]\n";
		if( child_len < pat_len ){
			cout << "HashTrieNode::getRange - Caso Simple\n";
			auto it_hash = hash_pat_map.find(child_len);
			if( it_hash == hash_pat_map.end() ){
				hash_pat = karp_rabin->hash(pattern, child_len);
				hash_pat_map[child_len] = hash_pat;
			}
			else{
				hash_pat = it_hash->second;
			}
			if( hash_pat == it_child.first ){
				cout << "HashTrieNode::getRange - Hijo encontrado (hash_pat: " << hash_pat << ", hash_child: " << it_child.first << ")\n";
				return it_child.second->getRange(pattern + child_len, pat_len - child_len, processed + child_len, karp_rabin, kr_factors);
			}
		}
		else{
			cout << "HashTrieNode::getRange - Caso complejo\n";
			// Si este nodo contiene el texto que queda del patron, este es el rango
			unsigned long long hash = kr_factors->hash(it_child.second->min_factor_pos, processed, pat_len);
			hash_pat = karp_rabin->hash(pattern, pat_len);
			if( hash == hash_pat ){
				cout << "HashTrieNode::getRange - Respuesta encontrada\n";
				return pair<unsigned int, unsigned int>(it_child.second->min, it_child.second->max);
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
	
	unsigned int n_lens = childs_lenghts.size();
	writer.write((char*)&n_lens, sizeof(int));
	for( unsigned int len : childs_lenghts ){
		writer.write((char*)&len, sizeof(int));
	}
	
	unsigned int n_childs = childs.size();
	writer.write((char*)&n_childs, sizeof(int));
	for( auto it : childs ){
		unsigned long long hash = it.first;
		writer.write((char*)&hash, sizeof(long long));
		it.second->save(writer);
	}
	
}

void HashTrieNode::load(fstream &reader){
	
	reader.read((char*)&len, sizeof(int));
	reader.read((char*)&min, sizeof(int));
	reader.read((char*)&max, sizeof(int));
	reader.read((char*)&min_factor_pos, sizeof(int));
	
	unsigned int n_lens = 0;
	reader.read((char*)&n_lens, sizeof(int));
	for(unsigned int i = 0; i < n_lens; ++i){
		unsigned int len = 0;
		reader.read((char*)&len, sizeof(int));
		childs_lenghts.push_back(len);
	}
	
	unsigned int n_childs = 0;
	reader.read((char*)&n_childs, sizeof(int));
	for(unsigned int i = 0; i < n_childs; ++i){
		unsigned long long hash = 0;
		reader.read((char*)&hash, sizeof(long long));
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

void HashTrie::load(const string &file){
	cout << "HashTrie::load - Start (" << file << ")\n";
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
	
	// Verificacion
	if( childs_pairs.size() != childs_lenghts.size() ){
		cout << "HashTrieNode::prepareChilds - Error in childs_pairs size\n";
	}
	for(unsigned int i = 0; i < childs_pairs.size(); ++i){
		if( childs_pairs[i].first != childs_lenghts[i] ){
			cout << "HashTrieNode::prepareChilds - Error in childs_pairs (pos " << i << ", " << childs_pairs[i].first << " != " << childs_lenghts[i] << ")\n";
		}
	}
	
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

HashTrieRev::HashTrieRev(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
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
	
	set<unsigned int> lengths;
	
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
		
		// El string es olo para el mensaje de debug
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
			lengths.insert(min_text_len);
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
	
//	cout << "HashTrieRevNode::build - Adding childs " << lengths.size() << " lentghs\n";
	for( unsigned int len : lengths ){
		childs_lenghts.push_back(len);
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
//	return pair<unsigned int, unsigned int>(0, 0);
	return root.getRange(pattern.c_str(), pattern.length(), 0, karp_rabin, kr_factors);
}

pair<unsigned int, unsigned int> HashTrieRevNode::getRange(const char *pattern, unsigned int pat_len, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors){
	
	cout << "HashTrieRevNode::getRange - Start (\"" << pattern << "\", " << pat_len << ", processed: " << processed << ")\n";
	
	if( pat_len == 0 ){
		return pair<unsigned int, unsigned int>(min, max);
	}
	
	map<unsigned int, unsigned long long> hash_pat_map;
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	
	// Si todos los hijos del nodo actual son menores que el largo del patron, se puede hacer el descenso rapido
	// Los largos de los hijos ESTAN ordenados de menor a mayor
//	cout << "HashTrieRevNode::getRange - Childs lengths: \n";
//	for( unsigned int len : childs_lenghts ){
//		cout << len << " - ";
//	}
//	cout << "\n";
	
	if( childs_lenghts.back() <= pat_len ){
		cout << "HashTrieRevNode::getRange - Quick Search\n";
		for( unsigned int len : childs_lenghts ){
			cout << "HashTrieRevNode::getRange - Probando len " << len << "\n";
			hash_pat = karp_rabin->hash(pattern, len);
			auto it_child = childs.find(hash_pat);
			if( it_child != childs.end() ){
				child_len = it_child->second->len;
				cout << "HashTrieRevNode::getRange - Hijo encontrado de len " << child_len << " -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return it_child->second->getRange(pattern + child_len, pat_len - child_len, processed + child_len, karp_rabin, kr_factors);
			}
		}
		cout << "HashTrieRevNode::getRange - Patron NO encontrado\n";
		return pair<unsigned int, unsigned int>(0, 0);
	}
	
	for( auto it_child : childs ){
		child_len = it_child.second->len;
		cout << "HashTrieRevNode::getRange - hijo de len " << child_len << " -> [" << it_child.second->min << ", " << it_child.second->max << "]\n";
		if( child_len < pat_len ){
			cout << "HashTrieRevNode::getRange - Caso Simple\n";
			auto it_hash = hash_pat_map.find(child_len);
			if( it_hash == hash_pat_map.end() ){
				hash_pat = karp_rabin->hash(pattern, child_len);
				hash_pat_map[child_len] = hash_pat;
			}
			else{
				hash_pat = it_hash->second;
			}
			if( hash_pat == it_child.first ){
				cout << "HashTrieRevNode::getRange - Hijo encontrado\n";
				return it_child.second->getRange(pattern + child_len, pat_len - child_len, processed + child_len, karp_rabin, kr_factors);
			}
		}
		else{
			cout << "HashTrieRevNode::getRange - Caso complejo\n";
			// Si este nodo contiene el texto que queda del patron, este es el rango
			unsigned long long hash = kr_factors->hash(it_child.second->min_factor_pos, processed, pat_len);
			hash_pat = karp_rabin->hash(pattern, pat_len);
			if( hash == hash_pat ){
				cout << "HashTrieRevNode::getRange - Respuesta encontrada\n";
				return pair<unsigned int, unsigned int>(it_child.second->min, it_child.second->max);
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
	
	unsigned int n_lens = childs_lenghts.size();
	writer.write((char*)&n_lens, sizeof(int));
	for( unsigned int len : childs_lenghts ){
		writer.write((char*)&len, sizeof(int));
	}
	
	unsigned int n_childs = childs.size();
	writer.write((char*)&n_childs, sizeof(int));
	for( auto it : childs ){
		unsigned long long hash = it.first;
		writer.write((char*)&hash, sizeof(long long));
		it.second->save(writer);
	}
	
}

void HashTrieRevNode::load(fstream &reader){
	
	reader.read((char*)&len, sizeof(int));
	reader.read((char*)&min, sizeof(int));
	reader.read((char*)&max, sizeof(int));
	reader.read((char*)&min_factor_pos, sizeof(int));
	
	unsigned int n_lens = 0;
	reader.read((char*)&n_lens, sizeof(int));
	for(unsigned int i = 0; i < n_lens; ++i){
		unsigned int len = 0;
		reader.read((char*)&len, sizeof(int));
		childs_lenghts.push_back(len);
	}
	
	unsigned int n_childs = 0;
	reader.read((char*)&n_childs, sizeof(int));
	for(unsigned int i = 0; i < n_childs; ++i){
		unsigned long long hash = 0;
		reader.read((char*)&hash, sizeof(long long));
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

void HashTrieRev::load(const string &file){
	cout << "HashTrieRev::load - Start (" << file << ")\n";
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
	
	// Verificacion
	if( childs_pairs.size() != childs_lenghts.size() ){
		cout << "HashTrieRevNode::prepareChilds - Error in childs_pairs size\n";
	}
	for(unsigned int i = 0; i < childs_pairs.size(); ++i){
		if( childs_pairs[i].first != childs_lenghts[i] ){
			cout << "HashTrieRevNode::prepareChilds - Error in childs_pairs (pos " << i << ", " << childs_pairs[i].first << " != " << childs_lenghts[i] << ")\n";
		}
	}
	
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











