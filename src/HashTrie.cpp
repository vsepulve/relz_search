#include "HashTrie.h"

HashTrie::HashTrie(){
	karp_rabin = NULL;
	factor_only = false;
}

HashTrie::HashTrie(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, bool _factor_only, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	factor_only = false;
	build(full_text, len_text, factors_start, arr_y, factor_only, _karp_rabin, _kr_factors);
}

HashTrie::~HashTrie(){
	karp_rabin = NULL;
}

void HashTrie::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, bool _factor_only, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	factor_only = _factor_only;
	
	cout << "HashTrie::build - Start (full text of " << len_text << ", " << factors_start.size() << " factors)\n";
	
	root.build(full_text, len_text, factors_start, arr_y, factor_only, karp_rabin, kr_factors, 0, factors_start.size()-1, 0);
	
	cout << "HashTrie::build - End\n";
}

HashTrieNode::HashTrieNode(){
//	hash = 0;
	len = 0;
	min = 0;
	max = 0;
//	factor = 0;
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

void HashTrieNode::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, bool factor_only, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len){
	if( min == max ){
		cout << "HashTrieNode::build - Leaf\n";
		return;
	}
	if( full_text == NULL || processed_len == len_text ){
		cout << "HashTrieNode::build - Null text\n";
		return;
	}
	
	cout << "HashTrieNode::build - Start (full text of " << len_text << ", range[" << min << ", " << max << "], processed_len: " << processed_len << ")\n";
	
	unsigned int min_pos = min;
	unsigned int min_text_start = factors_start[ arr_y[min_pos] ] + processed_len;
	unsigned int min_text_len = len_text - min_text_start;
	
	set<unsigned int> lengths;
	
	while( min_pos <= max ){
		
//		cout << "Preparing cur_pos = " << (min_pos + 1) << "\n";
		unsigned int cur_pos = min_pos + 1;
		unsigned int cur_text_start = factors_start[ arr_y[cur_pos] ] + processed_len;
		unsigned int cur_text_len = len_text - cur_text_start;
		
		unsigned int max_comp = 0;
		if( cur_pos < factors_start.size() ){
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
//			cout << "1\n";
			max_comp = getMaxComp(full_text + min_text_start, min_text_len, full_text + cur_text_start, cur_text_len);
//			cout << "2\n";
			
		}
	
		cout << "HashTrieNode::build - Preparing string s (min_text_start: " << min_text_start << " / " << len_text << ", len: " << min_text_len << ")\n";
		// El string es olo par ael mensaje de debug
		string s(full_text + min_text_start, min_text_len);
//		string s(full_text + min_text_start, (min_text_len<10)?min_text_len:10);
//		cout << "Preparing Hash\n";
//		unsigned long long hash = karp_rabin->hash(full_text, min_text_start, min_text_len, factors_start);
//		cout << "Preparing Hash for \"" << s << "\" (factor_ini: " << min_pos << " -> " << arr_y[min_pos] << ", offset: " << processed_len << ", length: " << min_text_len << ")\n";
		
		// Notar que aqui necesito la posicion del factor en la coleccion, no en el arreglo Y
//		cout << "3\n";
		unsigned long long hash = kr_factors->hash(arr_y[min_pos], processed_len, min_text_len);
//		cout << "Adding range [" << min_pos << ", " << cur_pos << "), len " << min_text_len << ", \"" << s << "\", hash: " << hash << " / " << karp_rabin->hash(s) << "\n";
//		cout << "Adding range [" << min_pos << ", " << cur_pos << "), len " << min_text_len << ", hash: " << hash << " / " << karp_rabin->hash(s) << "\n";
		cout << "HashTrieNode::build - Adding range [" << min_pos << ", " << cur_pos << "), len " << min_text_len << ", hash: " << hash << "\n";
		// Preparar el hijo, ejecutar la llamda sobre esa instancia
		childs[hash] = std::make_shared<HashTrieNode>();
//		HashTrieNode child;
//		childs[hash]->hash = hash;
		childs[hash]->len = min_text_len;
		childs[hash]->min = min_pos;
		childs[hash]->max = cur_pos-1;
		childs[hash]->min_factor_pos = arr_y[min_pos];
//		childs[hash] = child;
		lengths.insert(min_text_len);
		
		childs[hash]->build(full_text, len_text, factors_start, arr_y, factor_only, karp_rabin, kr_factors, min_pos, cur_pos-1, processed_len + min_text_len);
		
//		cout << "Preparing min_pos = " << cur_pos << "\n";
		min_pos = cur_pos;
		min_text_start = factors_start[ arr_y[min_pos] ] + processed_len;
		min_text_len = len_text - min_text_start;
//		cout << "Ok \n";
	
	}
	
	cout << "HashTrieNode::build - Adding childs " << lengths.size() << " lentghs\n";
	for( unsigned int len : lengths ){
		childs_lenghts.push_back(len);
	}
	
	
}

void HashTrie::print(){
	root.print(0);
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
	
	map<unsigned int, unsigned long long> hash_pat_map;
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	
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
				cout << "HashTrieNode::getRange - Hijo encontrado\n";
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
















