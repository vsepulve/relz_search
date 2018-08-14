#include "HashTrieChar.h"

HashTrieChar::HashTrieChar(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

HashTrieChar::HashTrieChar(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	build(full_text, len_text, factors_start, arr_y, _karp_rabin, _kr_factors);
}

HashTrieChar::~HashTrieChar(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

void HashTrieChar::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	
	cout << "HashTrieChar::build - Start (full text of " << len_text << ", " << factors_start.size() << " factors)\n";
	
	root.build(full_text, len_text, factors_start, arr_y, karp_rabin, kr_factors, 0, factors_start.size()-1, 0);
	
	cout << "HashTrieChar::build - End\n";
}

HashTrieCharNode::HashTrieCharNode(){
	len = 0;
	min = 0;
	max = 0;
	hash = 0;
	min_factor_pos = 0;
}

HashTrieCharNode::~HashTrieCharNode(){

}

//unsigned int HashTrieCharNode::getMaxComp(const char *str_1, unsigned int len_1, const char *str_2, unsigned int len_2){
////	cout << "HashTrieCharNode::getMaxComp - Start\n";
//	unsigned int ret = 0;
//	for( ; ret < ( (len_1<len_2)?len_1:len_2 ); ++ret ){
//		if( str_1[ret] != str_2[ret] ){
//			break;
//		}
//	}
////	cout << "HashTrieCharNode::getMaxComp - End (" << ret << " / " << ((len_1<len_2)?len_1:len_2) << ")\n";
//	return ret;
//}

void HashTrieCharNode::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len){
	if( min == max ){
//		cout << "HashTrieCharNode::build - Leaf\n";
		return;
	}
	if( full_text == NULL || processed_len == len_text ){
//		cout << "HashTrieCharNode::build - Null text\n";
		return;
	}
	
//	cout << "HashTrieCharNode::build - Start (full text of " << len_text << ", range[" << min << ", " << max << "], processed_len: " << processed_len << ")\n";
	
	unsigned int min_pos = min;
	unsigned int min_text_start = factors_start[ arr_y[min_pos] ] + processed_len;
	unsigned int min_text_len = 0;
	
	unsigned int cur_pos = 0;
	unsigned int cur_text_start = 0;
	
	unsigned long long hash;
	
	while( min_pos <= max ){
		
		// Primero busco el cur_pos MAYOR que tenga *primer* char igual al de min_pos
		// Despues reviso el largo mayor
		
//		cout << "HashTrieCharNode::build - Iniciando busqueda binaria en [" << min_pos << ", " << max << "]\n";
		
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
//		cout << "HashTrieCharNode::build - cur_pos: " << cur_pos << "\n";
		
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
//		cout << "HashTrieCharNode::build - min_text_len: " << min_text_len << "\n";
//		cout << "HashTrieCharNode::build - Preparing string s (min_text_start: " << min_text_start << " / " << len_text << ", len: " << min_text_len << ")\n";
		
		// El string es olo para el mensaje de debug
//		string s(full_text + min_text_start, min_text_len);
		char first_char = full_text[min_text_start];
		
		// Notar que aqui necesito la posicion del factor en la coleccion, no en el arreglo Y
		hash = kr_factors->hash(arr_y[min_pos], processed_len, min_text_len);
//		cout << "HashTrieCharNode::build - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_text_len << ", hash: " << hash << " / " << karp_rabin->hash(s) << " (\"" << s << "\")\n";
//		cout << "HashTrieCharNode::build - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_text_len << ", hash: " << hash << "\n";
		// Preparar el hijo, ejecutar la llamda sobre esa instancia
		if( min_text_len == 0 ){
//			cout << "HashTrieCharNode::build - Omiting child of len 0\n";
		}
		else{
			childs[first_char] = std::make_shared<HashTrieCharNode>();
			childs[first_char]->len = min_text_len;
			childs[first_char]->min = min_pos;
			childs[first_char]->max = cur_pos;
			childs[first_char]->hash = hash;
			childs[first_char]->min_factor_pos = arr_y[min_pos];
			childs[first_char]->build(full_text, len_text, factors_start, arr_y, karp_rabin, kr_factors, min_pos, cur_pos, processed_len + min_text_len);
		}
		
//		cout << "HashTrieCharNode::build - Preparing min_pos = " << cur_pos+1 << " / " << arr_y.size() << "\n";
		min_pos = cur_pos+1;
		if(min_pos >= arr_y.size()){
			break;
		}
		min_text_start = factors_start[ arr_y[min_pos] ] + processed_len;
		
//		cout << "HashTrieCharNode::build - min_text_start: " << min_text_start << "\n";
	
	}
	
//	cout << "----- \n";
	
}

void HashTrieChar::print(){
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

void HashTrieCharNode::print(unsigned int level){
	for(unsigned int i = 0; i < level; ++i){
		cout << "- ";
	}
	cout << "hash: " << hash << ", len: " << len << " , range [" << min << ", " << max << "]\n";
	for( auto it : childs ){
		it.second->print(level+1);
	}
}

pair<unsigned int, unsigned int> HashTrieChar::getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, const string &pattern){
	return root.getRange(kr_pat_vector, pos, 0, karp_rabin, kr_factors, pattern, &hash_nano);
}

pair<unsigned int, unsigned int> HashTrieCharNode::getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, const string &pattern, unsigned long long *hash_nano){
	
//	cout << "HashTrieCharNode::getRange - Start (prefixes: " << kr_pat_vector.size() << ", pos: " << pos << ", processed: " << processed << ")\n";
	
	if( pos + processed >= kr_pat_vector.size() ){
//		cout << "HashTrieCharNode::getRange - [" << min << ", " << max << "]\n";
		return pair<unsigned int, unsigned int>(min, max);
	}
	
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	unsigned int pat_len = kr_pat_vector.size() - pos - processed;
	char first_char_pat = pattern[pos + processed];
//	cout << "HashTrieCharNode::getRange - pat_len: " << pat_len << "\n";
//	string pat = pattern.substr(pos + processed, pat_len);
//	cout << "HashTrieCharNode::getRange - pat: " << pat << "\n";
	
	auto it_child = childs.find(first_char_pat);
	if(it_child != childs.end()){
		child_len = it_child->second->len;
		if( child_len <= pat_len ){
//			cout << "HashTrieCharNode::getRange - Case 1 " << child_len << "\n";
			hash_pat = karp_rabin->subtract_prefix(kr_pat_vector[pos + processed + child_len - 1], kr_pat_vector[pos + processed - 1], child_len);
			string pat_cut = pattern.substr(pos + processed, child_len);
//			cout << "HashTrieCharNode::getRange - pat_cut: " << pat_cut << " hash_pat: " << hash_pat << " / " << karp_rabin->hash(pat_cut) << " (processed: " << processed << ")\n";
			if( hash_pat == it_child->second->hash ){
//				cout << "HashTrieCharNode::getRange - Child found -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return it_child->second->getRange(kr_pat_vector, pos, processed + child_len, karp_rabin, kr_factors, pattern, hash_nano);
			}
		}
		else{
//			cout << "HashTrieCharNode::getRange - Case 2 " << child_len << "\n";
			unsigned long long hash = kr_factors->hashFast(it_child->second->min_factor_pos, processed, pat_len);
			hash_pat = karp_rabin->subtract_prefix(kr_pat_vector[kr_pat_vector.size() - 1], kr_pat_vector[pos + processed - 1], kr_pat_vector.size() - pos - processed);
//			cout << "HashTrieCharNode::getRange - hash: " << hash << ", hash_pat: " << hash_pat << "\n";
			if( hash == hash_pat ){
//				cout << "HashTrieCharNode::getRange - Child found -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return pair<unsigned int, unsigned int>(it_child->second->min, it_child->second->max);
			}
			
		}
	}
	
//	cout << "HashTrieCharNode::getRange - Patron NO encontrado\n";
	return pair<unsigned int, unsigned int>((unsigned int)(-1), (unsigned int)(-1));
}

void HashTrieCharNode::save(fstream &writer){
	
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

void HashTrieCharNode::load(fstream &reader, unsigned int processed, vector<unsigned int> &factors_start, const char *full_text){
	/*
	reader.read((char*)&len, sizeof(int));
	reader.read((char*)&min, sizeof(int));
	reader.read((char*)&max, sizeof(int));
	reader.read((char*)&min_factor_pos, sizeof(int));
	
	unsigned int n_childs = 0;
	reader.read((char*)&n_childs, sizeof(int));
	for(unsigned int i = 0; i < n_childs; ++i){
		unsigned int hash = 0;
		reader.read((char*)&hash, sizeof(int));
		childs[hash] = std::make_shared<HashTrieCharNode>();
		childs[hash]->load(reader);
	}
	*/
	
	// Version antigua
//	cout << "HashTrieCharNode::load - Original Version\n";
	reader.read((char*)&len, sizeof(int));
	reader.read((char*)&min, sizeof(int));
	reader.read((char*)&max, sizeof(int));
	reader.read((char*)&min_factor_pos, sizeof(int));
//	cout << "HashTrieCharNode::load - len: " << len << ", [" << min << ", " << max << "], min_factor_pos: " << min_factor_pos << "\n";
	
	unsigned int n_lens = 0;
	reader.read((char*)&n_lens, sizeof(int));
//	cout << "HashTrieCharNode::load - Omitiendo " << n_lens << " childs_lenghts\n";
	for(unsigned int i = 0; i < n_lens; ++i){
		unsigned int len = 0;
		reader.read((char*)&len, sizeof(int));
//		childs_lenghts.push_back(len);
	}
	
	unsigned int n_childs = 0;
	reader.read((char*)&n_childs, sizeof(int));
//	cout << "HashTrieCharNode::load - Cargando " << n_childs << " childs\n";
	vector<shared_ptr<HashTrieCharNode>> childs_vector;
	for(unsigned int i = 0; i < n_childs; ++i){
		unsigned long long hash = 0;
		reader.read((char*)&hash, sizeof(long long));
//		childs[hash] = std::make_shared<HashTrieCharNode>();
//		childs[hash]->load(reader);
		childs_vector.push_back(std::make_shared<HashTrieCharNode>());
		childs_vector.back()->load(reader, processed + len, factors_start, full_text);
		childs_vector.back()->hash = (unsigned int)hash;
	}
	
	// Lo que necesito aqui es el primer caracter de cada hijo, sacado de su min_factor_pos
	for( auto it_child : childs_vector ){
		// Calcular primer caracter
//		cout << "HashTrieCharNode::load - Preparando first_char desde it_child->min_factor_pos: " << it_child->min_factor_pos << ", pos_text: " << factors_start[ it_child->min_factor_pos ] << " + " << processed << " + " << len << "\n";
		unsigned int text_start = factors_start[ it_child->min_factor_pos ] + processed + len;
		char first_char = full_text[text_start];
//		cout << "HashTrieCharNode::load - first_char: " << first_char << "\n";
		if( childs.find(first_char) != childs.end() ){
			cout << "HashTrieCharNode::load - Warning\n";
		}
		childs[first_char] = it_child;
	}
//	cout << "-----\n";
}

void HashTrieChar::save(const string &file){
	cout << "HashTrieChar::save - Start (" << file << ")\n";
	fstream writer(file, fstream::out | fstream::trunc);
	root.save(writer);
	writer.close();
	cout << "HashTrieChar::save - End\n";
}

void HashTrieChar::load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, const string &file, vector<unsigned int> &factors_start, const char *full_text){
	cout << "HashTrieChar::load - Start (" << file << ")\n";
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	fstream reader(file, fstream::in);
	root.load(reader, 0, factors_start, full_text);
	reader.close();
	cout << "HashTrieChar::load - End\n";
}

/*
void HashTrieCharNode::prepareChilds(){
//	cout << "HashTrieCharNode::prepareChilds - Start\n";
//	root.prepareChilds();
//	unordered_map<unsigned int, shared_ptr<HashTrieCharNode>> childs;
//	vector< pair<unsigned int, vector<shared_ptr<HashTrieCharNode>>> > childs_pairs;
	childs_pairs.clear();
	map< unsigned int, vector<shared_ptr<HashTrieCharNode>> > childs_map;
	for( auto it_child : childs ){
		shared_ptr<HashTrieCharNode> ptr = it_child.second;
		unsigned int len = ptr->len;
		auto it_map = childs_map.find(len);
		if( it_map == childs_map.end() ){
			vector<shared_ptr<HashTrieCharNode>> childs_vector;
			childs_vector.push_back(ptr);
			childs_map[len] = childs_vector;
		}
		else{
			it_map->second.push_back(ptr);
		}
	}
	for( auto it_map : childs_map ){
		unsigned int len = it_map.first;
		vector<shared_ptr<HashTrieCharNode>> childs_vector = it_map.second;
		childs_pairs.push_back(
			pair<unsigned int, vector<shared_ptr<HashTrieCharNode>>>(
				len, childs_vector
			)
		);
	}
//	cout << "HashTrieCharNode::prepareChilds - Sorting pairs\n";
	sort(childs_pairs.begin(), childs_pairs.end());
	
//	cout << "HashTrieCharNode::prepareChilds - Resulting pairs:\n";
//	for( auto it_pair : childs_pairs ){
//		cout << "HashTrieCharNode::prepareChilds - (len: " << it_pair.first << ", n_childs: " << it_pair.second.size() << ")\n";
//	}
	
//	cout << "HashTrieCharNode::prepareChilds - Continuing with " << childs.size() << " childs\n";
	for( auto it_child : childs ){
		it_child.second->prepareChilds();
	}
	
//	cout << "HashTrieCharNode::prepareChilds - End\n";
}
*/

/*
void HashTrieChar::prepareChilds(){
	cout << "HashTrieChar::prepareChilds - Start\n";
	root.prepareChilds();
	cout << "HashTrieChar::prepareChilds - End\n";
}
*/
	
// VERSIONES REV

HashTrieCharRev::HashTrieCharRev(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

HashTrieCharRev::HashTrieCharRev(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	build(full_text, len_text, factors_start, arr_x, _karp_rabin, _kr_factors);
}

HashTrieCharRev::~HashTrieCharRev(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

void HashTrieCharRev::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	
	cout << "HashTrieCharRev::build - Start (full text of " << len_text << ", " << factors_start.size() << " factors)\n";
	
	root.build(full_text, len_text, factors_start, arr_x, karp_rabin, kr_factors, 0, factors_start.size()-1, 0);
	
	cout << "HashTrieCharRev::build - End\n";
}

HashTrieCharRevNode::HashTrieCharRevNode(){
	len = 0;
	min = 0;
	max = 0;
	hash = 0;
	min_factor_pos = 0;
}

HashTrieCharRevNode::~HashTrieCharRevNode(){

}

unsigned int HashTrieCharRevNode::getMaxComp(const char *str_1, unsigned int len_1, const char *str_2, unsigned int len_2){
//	cout << "HashTrieCharRevNode::getMaxComp - Start\n";
	unsigned int ret = 0;
	for( ; ret < ( (len_1<len_2)?len_1:len_2 ); ++ret ){
//		cout << "HashTrieCharRevNode::getMaxComp - " << *(str_1 - ret) << " vs " << *(str_2 - ret) << "\n";
		if( *(str_1 - ret) != *(str_2 - ret) ){
			break;
		}
	}
//	cout << "HashTrieCharRevNode::getMaxComp - End (" << ret << " / " << ((len_1<len_2)?len_1:len_2) << ")\n";
	return ret;
}

void HashTrieCharRevNode::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_x, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len){
	if( min == max ){
//		cout << "HashTrieCharRevNode::build - Leaf\n";
		return;
	}
	if( full_text == NULL || processed_len == len_text ){
//		cout << "HashTrieCharRevNode::build - Null text\n";
		return;
	}
	
//	cout << "HashTrieCharRevNode::build - Start (full text of " << len_text << ", range[" << min << ", " << max << "], processed_len: " << processed_len << ")\n";
	
	unsigned int min_pos = min;
	unsigned int min_text_start = 0;
	unsigned int min_text_len = 0;
	if( arr_x[min_pos] > 0 ){
		min_text_start = factors_start[ arr_x[min_pos] ] - processed_len - 1;
		min_text_len = factors_start[ arr_x[min_pos] ] - factors_start[ arr_x[min_pos] - 1 ] - processed_len;
	}
//	cout << "HashTrieCharRevNode::build - min_text_len: " << min_text_len << " (f " << arr_x[min_pos] << ")\n";
	
	unsigned long long hash;
	
	while( min_pos <= max ){
		
//		cout << "HashTrieCharRevNode::build - Preparing cur_pos = " << (min_pos + 1) << "\n";
		unsigned int cur_pos = min_pos + 1;
		unsigned int cur_text_start = 0;
		unsigned int cur_text_len = 0;
		if( arr_x[cur_pos] > 0 ){
			cur_text_start = factors_start[ arr_x[cur_pos] ] - processed_len - 1;
			cur_text_len = factors_start[ arr_x[cur_pos] ] - factors_start[ arr_x[cur_pos] - 1 ] - processed_len;
		}
//		cout << "HashTrieCharRevNode::build - cur_text_len: " << cur_text_len << " (f " << arr_x[cur_pos] << ")\n";
		
		unsigned int max_comp = 0;
		if( cur_pos < factors_start.size() && cur_pos <= max ){
//			cout << "HashTrieCharRevNode::build - min_text_start: " << min_text_start << ", cur_text_start: " << cur_text_start << "\n";
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
//			cout << "HashTrieCharRevNode::build - cur_text_len: " << cur_text_len << " (f " << arr_x[cur_pos] << ")\n";
//			cout << "HashTrieCharRevNode::build - min_text_start: " << min_text_start << ", cur_text_start: " << cur_text_start << "\n";
			max_comp = getMaxComp(full_text + min_text_start, min_text_len, full_text + cur_text_start, cur_text_len);
			
		}
		
//		cout << "HashTrieCharRevNode::build - Preparing string s (min_text_start: " << min_text_start << " / " << len_text << ", len: " << min_text_len << ")\n";
		
		// En el caso de Rev, uso s para el hash directo
		string s;
		for( unsigned int i = 0; i < min_text_len; ++i ){
			s += *(full_text + min_text_start - i);
		}
		hash = karp_rabin->hash(s);
		char first_char = s[0];
//		cout << "HashTrieCharRevNode::build - s: " << s << " (hash: " << hash << ")\n";
		
		// Notar que aqui necesito la posicion del factor en la coleccion, no en el arreglo Y
//		hash = kr_factors->hash(arr_x[min_pos], processed_len, min_text_len);
//		cout << "HashTrieCharRevNode::build - Adding range [" << min_pos << ", " << cur_pos << "), len " << min_text_len << ", hash: " << hash << " / " << karp_rabin->hash(s) << "\n";
//		cout << "HashTrieCharRevNode::build - Adding range [" << min_pos << ", " << cur_pos << "), len " << min_text_len << ", hash: " << hash << "\n";
		// Preparar el hijo, ejecutar la llamda sobre esa instancia
		if( min_text_len == 0 ){
//			cout << "HashTrieCharNode::build - Omiting child of len 0\n";
		}
		else{
			childs[first_char] = std::make_shared<HashTrieCharRevNode>();
			childs[first_char]->len = min_text_len;
			childs[first_char]->min = min_pos;
			childs[first_char]->max = cur_pos-1;
			childs[first_char]->min_factor_pos = arr_x[min_pos];
			childs[first_char]->text = s;
			childs[first_char]->hash = hash;
			childs[first_char]->build(full_text, len_text, factors_start, arr_x, karp_rabin, kr_factors, min_pos, cur_pos-1, processed_len + min_text_len);
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
//		cout << "HashTrieCharRevNode::build - min_text_len: " << min_text_len << " (f " << arr_x[min_pos] << ")\n";
	
	}
	
	
}

void HashTrieCharRev::print(){
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

void HashTrieCharRevNode::print(unsigned int level){
	for(unsigned int i = 0; i < level; ++i){
		cout << "- ";
	}
//	cout << "hash: " << hash << ", len: " << len << ", factor_base " << factor << " , range [" << min << ", " << max << "]\n";
	cout << "hash: " << hash << ", len: " << len << " , range [" << min << ", " << max << "]\n";
	for( auto it : childs ){
		it.second->print(level+1);
	}
}

pair<unsigned int, unsigned int> HashTrieCharRev::getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, const string &pattern_rev){
	return root.getRange(kr_pat_rev_vector, pos, 0, karp_rabin, kr_factors, pattern_rev);
}

pair<unsigned int, unsigned int> HashTrieCharRevNode::getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, const string &pattern_rev){
	
//	cout << "HashTrieCharRevNode::getRange - Start (prefixes: " << kr_pat_rev_vector.size() << ", pos: " << pos << ", processed: " << processed << ", pattern_rev: " << pattern_rev << ")\n";
	
	if( processed >= pos ){
//		cout << "HashTrieCharRevNode::getRange - [" << min << ", " << max << "]\n";
		return pair<unsigned int, unsigned int>(min, max);
	}
	
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	unsigned int pat_len = pos - processed;
	char first_char_pat = pattern_rev[kr_pat_rev_vector.size() - 1 - pos + processed];
//	cout << "HashTrieCharRevNode::getRange - pat = pattern_rev.substr(" << (kr_pat_rev_vector.size() - 1 - pos + processed) << ", " << pat_len << ");\n";
//	string pat = pattern_rev.substr(kr_pat_rev_vector.size() - 1 - pos + processed, pat_len);
//	cout << "HashTrieCharRevNode::getRange - pat: " << pat << ", first_char_pat: " << first_char_pat << "\n";
	
	auto it_child = childs.find(first_char_pat);
	if(it_child != childs.end()){
		child_len = it_child->second->len;
		if( child_len <= pat_len ){
//			cout << "HashTrieCharRevNode::getRange - Case 1 " << child_len << "\n";
			hash_pat = karp_rabin->hash(pattern_rev.c_str() + pattern_rev.length() - pos + processed, child_len);
			string pat_cut = pattern_rev.substr(kr_pat_rev_vector.size() - 1 - pos + processed, child_len);
//			cout << "HashTrieCharRevNode::getRange - pat_cut: " << pat_cut << " hash_pat: " << hash_pat << " / " << karp_rabin->hash(pat_cut) << " (processed: " << processed << ")\n";
			if( hash_pat == it_child->second->hash ){
//				cout << "HashTrieCharRevNode::getRange - Child found -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return it_child->second->getRange(kr_pat_rev_vector, pos, processed + child_len, karp_rabin, kr_factors, pattern_rev);
			}
		}
		else{
//			cout << "HashTrieCharRevNode::getRange - Case 2 " << child_len << "\n";
			unsigned long long hash = karp_rabin->hash(it_child->second->text.c_str(), pat_len);
			hash_pat = karp_rabin->hash(pattern_rev.c_str() + pattern_rev.length() - pos + processed, pat_len);
//			cout << "HashTrieCharRevNode::getRange - hash: " << hash << ", hash_pat: " << hash_pat << "\n";
			if( hash == hash_pat ){
//				cout << "HashTrieCharRevNode::getRange - Child found -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return pair<unsigned int, unsigned int>(it_child->second->min, it_child->second->max);
			}
			
		}
	}
	
//	cout << "HashTrieCharRevNode::getRange - Patron NO encontrado\n";
	return pair<unsigned int, unsigned int>((unsigned int)(-1), (unsigned int)(-1));
}

void HashTrieCharRevNode::save(fstream &writer){
	
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

void HashTrieCharRevNode::load(fstream &reader){
	
	reader.read((char*)&len, sizeof(int));
	reader.read((char*)&min, sizeof(int));
	reader.read((char*)&max, sizeof(int));
	reader.read((char*)&min_factor_pos, sizeof(int));
	
	unsigned int n_childs = 0;
	reader.read((char*)&n_childs, sizeof(int));
	for(unsigned int i = 0; i < n_childs; ++i){
		unsigned int hash = 0;
		reader.read((char*)&hash, sizeof(int));
		childs[hash] = std::make_shared<HashTrieCharRevNode>();
		childs[hash]->load(reader);
	}
	
}

void HashTrieCharRev::save(const string &file){
	cout << "HashTrieCharRev::save - Start (" << file << ")\n";
	fstream writer(file, fstream::out | fstream::trunc);
	root.save(writer);
	writer.close();
	cout << "HashTrieCharRev::save - End\n";
}
	
void HashTrieCharRev::load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, const string &file){
	cout << "HashTrieCharRev::load - Start (" << file << ")\n";
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	fstream reader(file, fstream::in);
	root.load(reader);
	reader.close();
	cout << "HashTrieCharRev::load - End\n";
}
/*
void HashTrieCharRevNode::prepareChilds(){
//	cout << "HashTrieCharRevNode::prepareChilds - Start\n";
//	root.prepareChilds();
//	unordered_map<unsigned int, shared_ptr<HashTrieCharRevNode>> childs;
//	vector< pair<unsigned int, vector<shared_ptr<HashTrieCharRevNode>>> > childs_pairs;
	childs_pairs.clear();
	map< unsigned int, vector<shared_ptr<HashTrieCharRevNode>> > childs_map;
	for( auto it_child : childs ){
		shared_ptr<HashTrieCharRevNode> ptr = it_child.second;
		unsigned int len = ptr->len;
		auto it_map = childs_map.find(len);
		if( it_map == childs_map.end() ){
			vector<shared_ptr<HashTrieCharRevNode>> childs_vector;
			childs_vector.push_back(ptr);
			childs_map[len] = childs_vector;
		}
		else{
			it_map->second.push_back(ptr);
		}
	}
	for( auto it_map : childs_map ){
		unsigned int len = it_map.first;
		vector<shared_ptr<HashTrieCharRevNode>> childs_vector = it_map.second;
		childs_pairs.push_back(
			pair<unsigned int, vector<shared_ptr<HashTrieCharRevNode>>>(
				len, childs_vector
			)
		);
	}
//	cout << "HashTrieCharRevNode::prepareChilds - Sorting pairs\n";
	sort(childs_pairs.begin(), childs_pairs.end());
	
//	cout << "HashTrieCharRevNode::prepareChilds - Resulting pairs:\n";
//	for( auto it_pair : childs_pairs ){
//		cout << "HashTrieCharRevNode::prepareChilds - (len: " << it_pair.first << ", n_childs: " << it_pair.second.size() << ")\n";
//	}
	
//	cout << "HashTrieCharRevNode::prepareChilds - Continuing with " << childs.size() << " childs\n";
	for( auto it_child : childs ){
		it_child.second->prepareChilds();
	}
	
//	cout << "HashTrieCharRevNode::prepareChilds - End\n";
}
*/

/*
void HashTrieCharRev::prepareChilds(){
	cout << "HashTrieCharRev::prepareChilds - Start\n";
	root.prepareChilds();
	cout << "HashTrieCharRev::prepareChilds - End\n";
}
*/











