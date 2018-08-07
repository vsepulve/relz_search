#include "HashTrie.h"

HashTrie::HashTrie(){
	karp_rabin = NULL;
	factor_only = false;
}

HashTrie::HashTrie(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, bool _factor_only, KarpRabin *_karp_rabin){
	karp_rabin = NULL;
	factor_only = false;
	build(full_text, len_text, factors_start, arr_y, factor_only, _karp_rabin);
}

HashTrie::~HashTrie(){
	karp_rabin = NULL;
}

void HashTrie::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, bool _factor_only, KarpRabin *_karp_rabin){
	karp_rabin = _karp_rabin;
	factor_only = _factor_only;
	
	cout << "HashTrie::build - Start (full text of " << len_text << ", " << factors_start.size() << " factors)\n";
	
	root.build(full_text, len_text, factors_start, arr_y, factor_only, karp_rabin, 0, factors_start.size()-1, 0);
	
//	unsigned int min_pos = 0;
//	unsigned int min_text_start = factors_start[ arr_y[min_pos] ];
//	unsigned int min_text_len = len_text - min_text_start;
//	
//	
//	while( min_pos < factors_start.size() ){
//	
//		unsigned int cur_pos = min_pos + 1;
//		unsigned int cur_text_start = factors_start[ arr_y[cur_pos] ];
//		unsigned int cur_text_len = len_text - cur_text_start;
//	
//		unsigned int max_comp = getMaxComp(full_text + min_text_start, min_text_len, full_text + cur_text_start, cur_text_len);
//	
//		while( max_comp > 0 ){
//			min_text_len = max_comp;
//			++cur_pos;
//			cur_text_start = factors_start[ arr_y[cur_pos] ];
//			cur_text_len = len_text - cur_text_start;
//			max_comp = getMaxComp(full_text + min_text_start, min_text_len, full_text + cur_text_start, cur_text_len);
//		}
//	
//		string s(full_text + min_text_start, min_text_len);
//		cout << "Adding range [" << min_pos << ", " << cur_pos << "), len " << min_text_len << ", " << s << "\n";
//		
//		
//		
//		min_pos = cur_pos;
//		min_text_start = factors_start[ arr_y[min_pos] ];
//		min_text_len = len_text - min_text_start;
//	
//	}
	
	cout << "HashTrie::build - End\n";
}

HashTrieNode::HashTrieNode(){
	hash = 0;
	len = 0;
	min = 0;
	max = 0;
	factor = 0;
}

HashTrieNode::~HashTrieNode(){

}

unsigned int HashTrieNode::getMaxComp(const char *str_1, unsigned int len_1, const char *str_2, unsigned int len_2){
	cout << "HashTrieNode::getMaxComp - Start\n";
	unsigned int ret = 0;
	for( ; ret < ( (len_1<len_2)?len_1:len_2 ); ++ret ){
		if( str_1[ret] != str_2[ret] ){
			break;
		}
	}
	cout << "HashTrieNode::getMaxComp - End (" << ret << " / " << ((len_1<len_2)?len_1:len_2) << ")\n";
	return ret;
}



void HashTrieNode::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, bool factor_only, KarpRabin *karp_rabin, unsigned int min, unsigned int max, unsigned int processed_len){
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
	
	while( min_pos <= max ){
		
		cout << "Preparing cur_pos = " << (min_pos + 1) << "\n";
		unsigned int cur_pos = min_pos + 1;
		unsigned int cur_text_start = factors_start[ arr_y[cur_pos] ] + processed_len;
		unsigned int cur_text_len = len_text - cur_text_start;
		
		unsigned int max_comp = 0;
		if( cur_pos < factors_start.size() ){
			max_comp = getMaxComp(full_text + min_text_start, min_text_len, full_text + cur_text_start, cur_text_len);
		}
		
		while( max_comp > 0 ){
			min_text_len = max_comp;
			++cur_pos;
			cur_text_start = factors_start[ arr_y[cur_pos] ] + processed_len;
			cur_text_len = len_text - cur_text_start - processed_len;
			max_comp = getMaxComp(full_text + min_text_start, min_text_len, full_text + cur_text_start, cur_text_len);
		}
	
		cout << "Preparing string s (min_text_start: " << min_text_start << " / " << len_text << ", len: " << min_text_len << ")\n";
		// El string es olo par ael mensaje de debug
		string s(full_text + min_text_start, (min_text_len<10)?min_text_len:10);
		cout << "Preparing Hash\n";
		unsigned long long hash = karp_rabin->hash(full_text, min_text_start, min_text_len, factors_start);
		cout << "Adding range [" << min_pos << ", " << cur_pos << "), len " << min_text_len << ", \"" << s << "\", hash: " << hash << "\n";
		// Preparar el hijo, ejecutar la llamda sobre esa instancia
		childs[hash] = std::make_shared<HashTrieNode>();
//		HashTrieNode child;
		childs[hash]->hash = hash;
		childs[hash]->len = min_text_len;
		childs[hash]->min = min_pos;
		childs[hash]->max = cur_pos-1;
		childs[hash]->factor = arr_y[min_pos];
//		childs[hash] = child;
		
		childs[hash]->build(full_text, len_text, factors_start, arr_y, factor_only, karp_rabin, min_pos, cur_pos-1, processed_len + min_text_len);
		
		cout << "Preparing min_pos = " << cur_pos << "\n";
		min_pos = cur_pos;
		min_text_start = factors_start[ arr_y[min_pos] ] + processed_len;
		min_text_len = len_text - min_text_start;
		cout << "Ok \n";
	
	}
	
}

void HashTrie::print(){
	root.print(0);
}

void HashTrieNode::print(unsigned int level){
	for(unsigned int i = 0; i < level; ++i){
		cout << "- ";
	}
	cout << "hash: " << hash << ", len: " << len << ", factor_base " << factor << " , range [" << min << ", " << max << "]\n";
	for( auto it : childs ){
		it.second->print(level+1);
	}
}













