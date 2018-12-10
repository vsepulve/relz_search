#include "HashTriev2.h"

HashTriev2::HashTriev2(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

HashTriev2::HashTriev2(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	arr_y = _arr_y;
	build(full_text, len_text, factors_start, arr_y, _karp_rabin, _kr_factors);
}

HashTriev2::~HashTriev2(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

void HashTriev2::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	arr_y = _arr_y;
	
	cout << "HashTriev2::build - Start (full text of " << len_text << ", " << factors_start.size() << " factors)\n";
	
	root.build(full_text, len_text, factors_start, arr_y, karp_rabin, kr_factors, 0, factors_start.size()-1, 0);
	
	cout << "HashTriev2::build - End\n";
}

HashTriev2Node::HashTriev2Node(){
	len = 0;
	min = 0;
//	max = 0;
	hash = 0;
//	min_factor_pos = 0;
}

HashTriev2Node::~HashTriev2Node(){

}

void HashTriev2Node::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *arr_y, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len){
	if( min == max ){
//		cout << "HashTriev2Node::build - Leaf\n";
		return;
	}
	if( full_text == NULL || processed_len == len_text ){
//		cout << "HashTriev2Node::build - Null text\n";
		return;
	}
	
//	cout << "HashTriev2Node::build - Start (full text of " << len_text << ", range[" << min << ", " << max << "], processed_len: " << processed_len << ")\n";
	
	unsigned int min_pos = min;
	unsigned int min_text_start = factors_start[ (*arr_y)[min_pos] ] + processed_len;
	unsigned int min_text_len = 0;
	
	unsigned int cur_pos = 0;
	unsigned int cur_text_start = 0;
	
	unsigned long long hash;
	
	while( min_pos <= max ){
		
		// Primero busco el cur_pos MAYOR que tenga *primer* char igual al de min_pos
		// Despues reviso el largo mayor
		
//		cout << "HashTriev2Node::build - Iniciando busqueda binaria en [" << min_pos << ", " << max << "]\n";
		
		unsigned int l = min_pos;
		unsigned int h = max;
		while(l < h){
			cur_pos = l + ((h-l)>>1);
			cur_text_start = factors_start[ (*arr_y)[cur_pos] ] + processed_len;
			if( full_text[cur_text_start] <= full_text[min_text_start] ){
				l = cur_pos+1;
			}
			else{
				h = cur_pos;
			}
		}
		cur_pos = h;
		cur_text_start = factors_start[ (*arr_y)[cur_pos] ] + processed_len;
		if( (cur_pos > min_pos) && full_text[cur_text_start] != full_text[min_text_start] ){
			--cur_pos;
			cur_text_start = factors_start[ (*arr_y)[cur_pos] ] + processed_len;
		}
//		cout << "HashTriev2Node::build - cur_pos: " << cur_pos << "\n";
		
		// Si cur_pos == min_pos, el largo es min_text_len
		// Si no, buscar el largo maximo entre min y cur
		
		if( cur_pos > min_pos ){
//			cout << "HashTriev2Node::build - Case 1\n";
			// Buscar el mayor
			// Aqui uso el hecho de que full_text termina en '\0'
			min_text_len = 1;
			while( full_text[min_text_start + min_text_len] == full_text[cur_text_start + min_text_len] ){
				++min_text_len;
			}
		}
		else{
//			cout << "HashTriev2Node::build - Case 2\n";
			min_text_len = len_text - min_text_start;
		}
//		cout << "HashTriev2Node::build - min_text_len: " << min_text_len << "\n";
//		cout << "HashTriev2Node::build - Preparing string s (min_text_start: " << min_text_start << " / " << len_text << ", len: " << min_text_len << ")\n";
		
		// El string es solo para el mensaje de debug
//		string s(full_text + min_text_start, min_text_len);
		char first_char = full_text[min_text_start];
		
		// Notar que aqui necesito la posicion del factor en la coleccion, no en el arreglo Y
		hash = kr_factors->hash( (*arr_y)[min_pos], processed_len, min_text_len);
//		cout << "HashTriev2Node::build - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_text_len << ", hash: " << hash << " / " << karp_rabin->hash(s) << " (\"" << s << "\")\n";
//		cout << "HashTriev2Node::build - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_text_len << ", hash: " << hash << "\n";
		// Preparar el hijo, ejecutar la llamda sobre esa instancia
		if( min_text_len == 0 ){
//			cout << "HashTriev2Node::build - Omiting child of len 0\n";
		}
		else{
			childs_vector.push_back(HashTriev2Node());
			childs_vector.back().first = first_char;
			childs_vector.back().len = min_text_len;
			childs_vector.back().min = min_pos;
			childs_vector.back().hash = hash;
			childs_vector.back().build(full_text, len_text, factors_start, arr_y, karp_rabin, kr_factors, min_pos, cur_pos, processed_len + min_text_len);
		}
		
//		cout << "HashTriev2Node::build - Preparing min_pos = " << cur_pos+1 << " / " << arr_y->size() << "\n";
		min_pos = cur_pos+1;
		if(min_pos >= arr_y->size()){
			break;
		}
		min_text_start = factors_start[ (*arr_y)[min_pos] ] + processed_len;
		
//		cout << "HashTriev2Node::build - min_text_start: " << min_text_start << "\n";
	
	}
	
	std::sort(childs_vector.begin(), childs_vector.end(), 
		[](const HashTriev2Node &p1, const HashTriev2Node &p2) -> bool{return p1.first < p2.first;}
	);
	
//	cout << "----- \n";
	
}

void HashTriev2::print(){
	root.print(0);
	cout << "Root: " << root.childs_vector.size() << " childs\n";
	cout << "Lengths: ";
	for( auto it_child : root.childs_vector ){
		cout << it_child.len << " - ";
	}
	cout << "\n";
}

void HashTriev2Node::print(unsigned int level){
	for(unsigned int i = 0; i < level; ++i){
		cout << "- ";
	}
	cout << "hash: " << hash << ", len: " << len << " , range [" << min << "]\n";
	for( auto it : childs_vector ){
		cout << it.first << " ";
		it.print(level+1);
	}
}

unsigned int HashTriev2Node::totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height){
	unsigned int ret = childs_vector.size();
	if( height > max_height ){
		max_height = height;
	}
	if( len > max_len ){
		max_len = len;
	}
	if( childs_vector.size() > max_childs ){
		max_childs = childs_vector.size();
	}
	for( auto it : childs_vector ){
		ret += it.totalChilds(max_len, max_childs, max_height, 1 + height);
	}
	return ret;
}

void HashTriev2::printSize(){
	unsigned int max_len = 0;
	unsigned int max_childs = 0;
	unsigned int max_height = 0;
	unsigned int total_childs = root.totalChilds(max_len, max_childs, max_height, 0);
	cout << "HashTriev2::printSize - totalChilds: " << total_childs << ", max_len: " << max_len << ", max_childs: " << max_childs << ", max_height: " << max_height << "\n";
}

pair<unsigned int, unsigned int> HashTriev2::getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, const string &pattern){
//	return root.getRange(kr_pat_vector, pos, 0, karp_rabin, kr_factors, arr_y, pattern, &hash_nano);
	return root.getRange(kr_pat_vector, pos, 0, karp_rabin, kr_factors, arr_y, arr_y->size() - 1, pattern, &hash_nano);
}

unsigned int HashTriev2Node::findChild(char c){
	
//	if(childs_vector.size() < 1){
//		return NOT_FOUND;
//	}
//	// Busqueda Binaria
	
	for(unsigned int i = 0; i < childs_vector.size(); ++i){
		if( c == childs_vector[i].first ){
			return i;
		}
	}
	return NOT_FOUND;
}

pair<unsigned int, unsigned int> HashTriev2Node::getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_y, unsigned int cur_max, const string &pattern, unsigned long long *hash_nano){

//pair<unsigned int, unsigned int> HashTriev2Node::getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_y, const string &pattern, unsigned long long *hash_nano){
	
//	cout << "HashTriev2Node::getRange - Start (prefixes: " << kr_pat_vector.size() << ", pos: " << pos << ", processed: " << processed << ")\n";
	
	if( pos + processed >= kr_pat_vector.size() ){
//		cout << "HashTriev2Node::getRange - [" << min << ", " << cur_max << "]\n";
//		return pair<unsigned int, unsigned int>(min, max);
		return pair<unsigned int, unsigned int>(min, cur_max);
	}
	
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	unsigned int pat_len = kr_pat_vector.size() - pos - processed;
	char first_char_pat = pattern[pos + processed];
//	cout << "HashTriev2Node::getRange - pat_len: " << pat_len << "\n";
	string pat = pattern.substr(pos + processed, pat_len);
//	cout << "HashTriev2Node::getRange - pat: " << pat << "\n";
	
	unsigned int pos_child = findChild(first_char_pat);
	if( pos_child != NOT_FOUND ){
		child_len = childs_vector[pos_child].len;
		
		// Ajuste a cur_max
		if( pos_child < childs_vector.size() - 1 ){
			cur_max = childs_vector[pos_child+1].min - 1;
		}
		
		if( child_len <= pat_len ){
//			cout << "HashTriev2Node::getRange - Case 1, child_len: " << child_len << "\n";
			hash_pat = karp_rabin->subtract_prefix(kr_pat_vector[pos + processed + child_len - 1], kr_pat_vector[pos + processed - 1], child_len);
			string pat_cut = pattern.substr(pos + processed, child_len);
//			cout << "HashTriev2Node::getRange - pat_cut: " << pat_cut << " hash_pat: " << hash_pat << " / " << karp_rabin->hash(pat_cut) << " (processed: " << processed << ")\n";
			if( hash_pat == childs_vector[pos_child].hash ){
//				cout << "HashTriev2Node::getRange - Child found -> [" << childs_vector[pos_child].min << ", " << cur_max << "]\n";
//				return childs_vector[pos_child].getRange(kr_pat_vector, pos, processed + child_len, karp_rabin, kr_factors, arr_y, pattern, hash_nano);
				return childs_vector[pos_child].getRange(kr_pat_vector, pos, processed + child_len, karp_rabin, kr_factors, arr_y, cur_max, pattern, hash_nano);
			}
		}
		else{
//			cout << "HashTriev2Node::getRange - Case 2, child_len: " << child_len << "\n";
			
			unsigned int min_factor_pos = (*arr_y)[childs_vector[pos_child].min];
			
//			cout << "Node: \n";
//			childs_vector[pos_child].print(0);
			// Caso de borde detectado a veces en kr_factors->hashFast(childs_vector[pos_child].min_factor_pos, processed, pat_len
			
			unsigned long long hash = kr_factors->hashFast(min_factor_pos, processed, pat_len);
			
			hash_pat = karp_rabin->subtract_prefix(kr_pat_vector[kr_pat_vector.size() - 1], kr_pat_vector[pos + processed - 1], kr_pat_vector.size() - pos - processed);
//			cout << "HashTriev2Node::getRange - hash: " << hash << ", hash_pat: " << hash_pat << "\n";
			if( hash == hash_pat ){
//				cout << "HashTriev2Node::getRange - Child found -> [" << childs_vector[pos_child].min << ", " << cur_max << "]\n";
//				return pair<unsigned int, unsigned int>(childs_vector[pos_child].min, childs_vector[pos_child].max);
				return pair<unsigned int, unsigned int>(childs_vector[pos_child].min, cur_max);
			}
		}
	}
	
//	cout << "HashTriev2Node::getRange - Pattern NOT found\n";
	return pair<unsigned int, unsigned int>((unsigned int)(-1), (unsigned int)(-1));
}

void HashTriev2Node::save(fstream &writer){
	
	writer.write((char*)&len, sizeof(int));
	writer.write((char*)&min, sizeof(int));
	unsigned int max = 0;
	writer.write((char*)&max, sizeof(int));
	writer.write((char*)&hash, sizeof(int));
	unsigned int min_factor_pos = 0;
	writer.write((char*)&min_factor_pos, sizeof(int));
//	writer.write((char*)&first, 1);
	
	unsigned int n_childs = childs_vector.size();
	writer.write((char*)&n_childs, sizeof(int));
	for( auto it : childs_vector ){
		writer.write((char*)&(it.first), 1);
		it.save(writer);
	}
	
}

void HashTriev2Node::load(fstream &reader){
	
	reader.read((char*)&len, sizeof(int));
	reader.read((char*)&min, sizeof(int));
	// Lo dejo temporalmente por legacy
	unsigned int max = 0;
	reader.read((char*)&max, sizeof(int));
	reader.read((char*)&hash, sizeof(int));
	// Lo dejo temporalmente por legacy
	unsigned int min_factor_pos = 0;
	reader.read((char*)&min_factor_pos, sizeof(int));
	
	unsigned int n_childs = 0;
	reader.read((char*)&n_childs, sizeof(int));
	for(unsigned int i = 0; i < n_childs; ++i){
		
		// Lo que sigue debe cargarse solo cuando la clase este terminada
		
		char child_first_char = 0;
		reader.read((char*)&child_first_char, 1);
		
//		childs[child_first_char] = std::make_shared<HashTriev2Node>();
//		childs[child_first_char]->load(reader);
		
		childs_vector.push_back(HashTriev2Node());
		childs_vector.back().load(reader);
		childs_vector.back().first = child_first_char;
		
	}
	
	std::sort(childs_vector.begin(), childs_vector.end(), 
		[](const HashTriev2Node &p1, const HashTriev2Node &p2) -> bool{return p1.first < p2.first;}
	);
	
	// Debug
//	cout << "HashTriev2Node::load - Childs: " << childs_vector.size() << " / ";
//	for(auto it : childs_vector){
//		cout << it.first << " ";
//	}
//	cout << "\n";

}

void HashTriev2::save(const string &file){
	cout << "HashTriev2::save - Start (" << file << ")\n";
	fstream writer(file, fstream::out | fstream::trunc);
	root.save(writer);
	writer.close();
	cout << "HashTriev2::save - End\n";
}

void HashTriev2::load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, int_vector<> *_arr_y, const string &file){
	cout << "HashTriev2::load - Start (" << file << ")\n";
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	arr_y = _arr_y;
	fstream reader(file, fstream::in);
	root.load(reader);
	reader.close();
	cout << "HashTriev2::load - End\n";
}
	
// VERSIONES REV

HashTriev2Rev::HashTriev2Rev(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

HashTriev2Rev::HashTriev2Rev(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	arr_x = _arr_x;
	build(full_text, len_text, factors_start, arr_x, _karp_rabin, _kr_factors);
}

HashTriev2Rev::~HashTriev2Rev(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

void HashTriev2Rev::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	arr_x = _arr_x;
	
	cout << "HashTriev2Rev::build - Start (full text of " << len_text << ", " << factors_start.size() << " factors)\n";
	
	root.build(full_text, len_text, factors_start, arr_x, karp_rabin, kr_factors, 0, factors_start.size()-1, 0);
	
	cout << "HashTriev2Rev::build - End\n";
}

HashTriev2RevNode::HashTriev2RevNode(){
	len = 0;
	min = 0;
//	max = 0;
	hash = 0;
//	min_factor_pos = 0;
}

HashTriev2RevNode::~HashTriev2RevNode(){

}

void HashTriev2RevNode::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *arr_x, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len){
	if( min == max ){
//		cout << "HashTriev2RevNode::build - Leaf\n";
		return;
	}
	if( full_text == NULL || processed_len == len_text ){
//		cout << "HashTriev2RevNode::build - Null text\n";
		return;
	}
	
//	cout << "HashTriev2RevNode::build - Start (full text of " << len_text << ", range[" << min << ", " << max << "], processed_len: " << processed_len << ")\n";
	
	unsigned int min_pos = min;
	
	unsigned int min_text_len = 0;
	if( (*arr_x)[min_pos] > 0 ){
		min_text_len = factors_start[ (*arr_x)[min_pos] ] - factors_start[ (*arr_x)[min_pos] - 1 ] - processed_len;
	}
	while( min_text_len == 0 && min_pos < max ){
		++min_pos;
		if( (*arr_x)[min_pos] > 0 ){
			min_text_len = factors_start[ (*arr_x)[min_pos] ] - factors_start[ (*arr_x)[min_pos] - 1 ] - processed_len;
		}
	}
	
	unsigned int min_text_start = 0;
	if( (*arr_x)[min_pos] > 0 ){
		min_text_start = factors_start[ (*arr_x)[min_pos] ] - processed_len - 1;
	}
	
	unsigned int cur_pos = 0;
	unsigned int cur_text_start = 0;
	unsigned int cur_text_len = 0;
	unsigned int min_common_text = 0;
	
	unsigned long long hash;
	
	while( min_pos <= max ){
		
		// Primero busco el cur_pos MAYOR que tenga *primer* char igual al de min_pos
		// Despues reviso el largo mayor
		
//		cout << "HashTriev2RevNode::build - Iniciando busqueda binaria en [" << min_pos << ", " << max << "]\n";
		
		unsigned int l = min_pos;
		unsigned int h = max;
		while(l < h){
			cur_pos = l + ((h-l)>>1);
			cur_text_start = factors_start[ (*arr_x)[cur_pos] ] - processed_len - 1;
			if( (*arr_x)[cur_pos] == 0 ){
				cur_text_len = 0;
			}
			else{
				cur_text_len = factors_start[ (*arr_x)[cur_pos] ] - factors_start[ (*arr_x)[cur_pos] - 1 ] - processed_len;
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
		cur_text_start = factors_start[ (*arr_x)[cur_pos] ] - processed_len - 1;
		if( (*arr_x)[cur_pos] == 0 ){
			cur_text_len = 0;
		}
		else{
			cur_text_len = factors_start[ (*arr_x)[cur_pos] ] - factors_start[ (*arr_x)[cur_pos] - 1 ] - processed_len;
		}
		if( (cur_pos > min_pos) && full_text[cur_text_start] != full_text[min_text_start] ){
			--cur_pos;
			cur_text_start = factors_start[ (*arr_x)[cur_pos] ] - processed_len - 1;
			if( (*arr_x)[cur_pos] == 0 ){
				cur_text_len = 0;
			}
			else{
				cur_text_len = factors_start[ (*arr_x)[cur_pos] ] - factors_start[ (*arr_x)[cur_pos] - 1 ] - processed_len;
			}
		}
//		cout << "HashTriev2RevNode::build - cur_pos: " << cur_pos << "\n";
		
		// Si cur_pos == min_pos, el largo es min_text_len
		// Si no, buscar el largo maximo entre min y cur
		
		if( cur_pos > min_pos ){
//			cout << "HashTriev2RevNode::build - Case 1\n";
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
//			cout << "HashTriev2RevNode::build - Case 2\n";
			min_common_text = min_text_len;
		}
//		cout << "HashTriev2RevNode::build - min_common_text: " << min_common_text << "\n";
//		cout << "HashTriev2RevNode::build - Preparing string s (min_text_start: " << min_text_start << " / " << len_text << ", len: " << min_common_text << ")\n";
		
		// En el caso de Rev, uso s para el hash directo
		string s;
		for( unsigned int i = 0; i < min_common_text; ++i ){
			s += *(full_text + min_text_start - i);
		}
		hash = karp_rabin->hash(s);
//		cout << "HashTriev2RevNode::build - s: " << s << "\n";
		char first_char = s[0];
//		cout << "HashTriev2RevNode::build - hash: " << hash << "\n";
		
//		cout << "HashTriev2RevNode::build - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_common_text << ", hash: " << hash << "\n";
		// Preparar el hijo, ejecutar la llamda sobre esa instancia
		if( min_common_text == 0 ){
//			cout << "HashTriev2RevNode::build - Omiting child of len 0\n";
		}
		else{
//			childs[first_char] = std::make_shared<HashTriev2RevNode>();
//			childs[first_char]->len = min_common_text;
//			childs[first_char]->min = min_pos;
////			childs[first_char]->max = cur_pos;
//			childs[first_char]->hash = hash;
////			childs[first_char]->min_factor_pos = (*arr_x)[min_pos];
//			childs[first_char]->build(full_text, len_text, factors_start, arr_x, karp_rabin, kr_factors, min_pos, cur_pos, processed_len + min_common_text);
			
			childs_vector.push_back(HashTriev2RevNode());
			childs_vector.back().first = first_char;
			childs_vector.back().len = min_common_text;
			childs_vector.back().min = min_pos;
			childs_vector.back().hash = hash;
			childs_vector.back().build(full_text, len_text, factors_start, arr_x, karp_rabin, kr_factors, min_pos, cur_pos, processed_len + min_common_text);
			
		}
		
//		cout << "HashTriev2RevNode::build - Preparing min_pos = " << cur_pos+1 << " / " << arr_x->size() << "\n";
		min_pos = cur_pos+1;
		if(min_pos >= arr_x->size()){
			break;
		}
		min_text_start = 0;
		min_text_len = 0;
		if( (*arr_x)[min_pos] > 0 ){
			min_text_start = factors_start[ (*arr_x)[min_pos] ] - processed_len - 1;
			min_text_len = factors_start[ (*arr_x)[min_pos] ] - factors_start[ (*arr_x)[min_pos] - 1 ] - processed_len;
		}
//		cout << "HashTriev2RevNode::build - min_text_start: " << min_text_start << " de len " << min_text_len << "\n";
	
	}
	
	std::sort(childs_vector.begin(), childs_vector.end(), 
		[](const HashTriev2RevNode &p1, const HashTriev2RevNode &p2) -> bool{return p1.first < p2.first;}
	);
	
//	cout << "----- \n";
	
}


void HashTriev2Rev::print(){
	root.print(0);
	cout << "Root: " << root.childs_vector.size() << " childs\n";
	cout << "Lengths: ";
	for( auto it_child : root.childs_vector ){
		cout << it_child.len << " - ";
	}
	cout << "\n";
}

void HashTriev2RevNode::print(unsigned int level){
	for(unsigned int i = 0; i < level; ++i){
		cout << "- ";
	}
//	cout << "hash: " << hash << ", len: " << len << ", factor_base " << factor << " , range [" << min << ", " << max << "]\n";
	cout << "hash: " << hash << ", len: " << len << " , range [" << min << "]\n";
	for( auto it : childs_vector ){
		cout << it.first << " ";
		it.print(level+1);
	}
}

unsigned int HashTriev2RevNode::totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height){
	unsigned int ret = childs_vector.size();
	if( height > max_height ){
		max_height = height;
	}
	if( len > max_len ){
		max_len = len;
	}
	if( childs_vector.size() > max_childs ){
		max_childs = childs_vector.size();
	}
	for( auto it : childs_vector ){
		ret += it.totalChilds(max_len, max_childs, max_height, 1 + height);
	}
	return ret;
}

void HashTriev2Rev::printSize(){
	unsigned int max_len = 0;
	unsigned int max_childs = 0;
	unsigned int max_height = 0;
	unsigned int total_childs = root.totalChilds(max_len, max_childs, max_height, 0);
	cout << "HashTriev2Rev::printSize - totalChilds: " << total_childs << ", max_len: " << max_len << ", max_childs: " << max_childs << ", max_height: " << max_height << "\n";
}

pair<unsigned int, unsigned int> HashTriev2Rev::getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, const string &pattern_rev){
	return root.getRange(kr_pat_rev_vector, pos, 0, karp_rabin, kr_factors, arr_x, arr_x->size() - 1, pattern_rev);
}

unsigned int HashTriev2RevNode::findChild(char c){
	
//	if(childs_vector.size() < 1){
//		return NOT_FOUND;
//	}
//	// Busqueda Binaria
	
//	cout << "HashTriev2RevNode::findChild - Char " << c << ", " << childs_vector.size() << " childs ( ";
//	for(unsigned int i = 0; i < childs_vector.size(); ++i){
//		cout << childs_vector[i].first << " ";
//	}
//	cout << ")\n";
	
	for(unsigned int i = 0; i < childs_vector.size(); ++i){
		if( c == childs_vector[i].first ){
			return i;
		}
	}
	return NOT_FOUND;
}

pair<unsigned int, unsigned int> HashTriev2RevNode::getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_x, unsigned int cur_max, const string &pattern_rev){
	
//	cout << "HashTriev2RevNode::getRange - Start (prefixes: " << kr_pat_rev_vector.size() << ", pos: " << pos << ", processed: " << processed << ", pattern_rev: " << pattern_rev << ")\n";
	
	if( processed >= pos ){
//		cout << "HashTriev2RevNode::getRange - [" << min << ", " << cur_max << "]\n";
		return pair<unsigned int, unsigned int>(min, cur_max);
	}
	
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	unsigned int pat_len = pos - processed;
	char first_char_pat = pattern_rev[kr_pat_rev_vector.size() - 1 - pos + processed];
//	cout << "HashTriev2RevNode::getRange - pat = pattern_rev.substr(" << (kr_pat_rev_vector.size() - 1 - pos + processed) << ", " << pat_len << ");\n";
	string pat = pattern_rev.substr(kr_pat_rev_vector.size() - 1 - pos + processed, pat_len);
//	cout << "HashTriev2RevNode::getRange - pat: " << pat << ", first_char_pat: " << first_char_pat << "\n";
	
	
	
	unsigned int pos_child = findChild(first_char_pat);
	if( pos_child != NOT_FOUND ){
		child_len = childs_vector[pos_child].len;
		
		// Ajuste a cur_max
//		cout << "HashTriev2RevNode::getRange - Adjusting cur_max (pos_child: " << pos_child << " / " << childs_vector.size() << ")\n";
		if( pos_child < childs_vector.size() - 1 ){
			cur_max = childs_vector[pos_child+1].min - 1;
//			cout << "HashTriev2RevNode::getRange - cur_max: " << cur_max << ")\n";
		}
		
		if( child_len <= pat_len ){
//			cout << "HashTriev2RevNode::getRange - Case 1, child_len: " << child_len << "\n";
			hash_pat = karp_rabin->hash(pattern_rev.c_str() + pattern_rev.length() - pos + processed, child_len);
			string pat_cut = pattern_rev.substr(kr_pat_rev_vector.size() - 1 - pos + processed, child_len);
//			cout << "HashTriev2RevNode::getRange - pat_cut: " << pat_cut << " hash_pat: " << hash_pat << " / " << karp_rabin->hash(pat_cut) << " (processed: " << processed << ")\n";
			if( hash_pat == childs_vector[pos_child].hash ){
//				cout << "HashTriev2RevNode::getRange - Child found -> [" << childs_vector[pos_child].min << ", " << cur_max << "]\n";
				return childs_vector[pos_child].getRange(kr_pat_rev_vector, pos, processed + child_len, karp_rabin, kr_factors, arr_x, cur_max, pattern_rev);
			}
		}
		else{
//			cout << "HashTriev2RevNode::getRange - Case 2, child_len: " << child_len << "\n";
			
			string test_text = "";
			unsigned int min_factor_pos = (*arr_x)[childs_vector[pos_child].min];
			
			if( min_factor_pos > 0 ){
				KarpRabinFactorsSuffixesv2 *ptr = static_cast<KarpRabinFactorsSuffixesv2*>(kr_factors);
				
				unsigned int cur_pi = (*(ptr->pi_inv))[min_factor_pos-1];
				unsigned int tu = ptr->select1_s->operator()(cur_pi + 1) - cur_pi;
				unsigned int pu = ptr->select1_b->operator()(min_factor_pos-1 + 1);
				unsigned int lu = ptr->select1_b->operator()(min_factor_pos-1 + 2) - pu;
				
				if(processed < lu){
					unsigned int len = lu - processed;
					if( childs_vector[pos_child].len < len ){
						len = childs_vector[pos_child].len;
					}
					if( pat_len < len ){
						len = pat_len;
					}
//					cout << "HashTriev2RevNode::getRange - Adding " << len << " chars\n";
					for(unsigned int i = 0; i < len; ++i){
						test_text += *(ptr->ref_text + tu + lu - processed - i - 1);
					}
				}
			}
			
			unsigned long long hash = karp_rabin->hash(test_text.c_str(), pat_len);
			hash_pat = karp_rabin->hash(pattern_rev.c_str() + pattern_rev.length() - pos + processed, pat_len);
//			cout << "HashTriev2RevNode::getRange - hash: " << hash << ", hash_pat: " << hash_pat << "\n";
			if( hash == hash_pat ){
//				cout << "HashTriev2RevNode::getRange - Child found -> [" << childs_vector[pos_child].min << ", " << cur_max << "]\n";
				return pair<unsigned int, unsigned int>(childs_vector[pos_child].min, cur_max);
			}
			
		}
	}
	
	
	
	/*
	auto it_child = childs.find(first_char_pat);
	if(it_child != childs.end()){
		child_len = it_child->second->len;
		if( child_len <= pat_len ){
//			cout << "HashTriev2RevNode::getRange - Case 1, child_len: " << child_len << "\n";
			hash_pat = karp_rabin->hash(pattern_rev.c_str() + pattern_rev.length() - pos + processed, child_len);
			string pat_cut = pattern_rev.substr(kr_pat_rev_vector.size() - 1 - pos + processed, child_len);
//			cout << "HashTriev2RevNode::getRange - pat_cut: " << pat_cut << " hash_pat: " << hash_pat << " / " << karp_rabin->hash(pat_cut) << " (processed: " << processed << ")\n";
			if( hash_pat == it_child->second->hash ){
//				cout << "HashTriev2RevNode::getRange - Child found -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return it_child->second->getRange(kr_pat_rev_vector, pos, processed + child_len, karp_rabin, kr_factors, arr_x, pattern_rev);
			}
		}
		else{
//			cout << "HashTriev2RevNode::getRange - Case 2, child_len: " << child_len << "\n";
			
			
			string test_text = "";
			unsigned int min_factor_pos = (*arr_x)[it_child->second->min];
			
			if( min_factor_pos > 0 ){
				KarpRabinFactorsSuffixesv2 *ptr = static_cast<KarpRabinFactorsSuffixesv2*>(kr_factors);
				
				unsigned int cur_pi = (*(ptr->pi_inv))[min_factor_pos-1];
				unsigned int tu = ptr->select1_s->operator()(cur_pi + 1) - cur_pi;
				unsigned int pu = ptr->select1_b->operator()(min_factor_pos-1 + 1);
				unsigned int lu = ptr->select1_b->operator()(min_factor_pos-1 + 2) - pu;
				
				if(processed < lu){
					unsigned int len = lu - processed;
					if( it_child->second->len < len ){
						len = it_child->second->len;
					}
					if( pat_len < len ){
						len = pat_len;
					}
//					cout << "HashTriev2RevNode::getRange - Adding " << len << " chars\n";
					for(unsigned int i = 0; i < len; ++i){
						test_text += *(ptr->ref_text + tu + lu - processed - i - 1);
					}
				}
			}
			
//			cout << "HashTriev2RevNode::getRange - Factor (rev): " << min_factor_pos << "-1, length " << test_text.length() << " / " << it_child->second->len << "\n";
//			cout << "HashTriev2RevNode::getRange - text: " << test_text << " / " << it_child->second->text << "\n";
//			if( test_text.length() != it_child->second->text.length() ){
//				cout << "HashTriev2RevNode::getRange - Error\n";
//				exit(0);
//			}
//			for(unsigned int i = 0; i < test_text.length(); ++i){
//				if( test_text[i] != it_child->second->text[i] ){
//					cout << "HashTriev2RevNode::getRange - Error\n";
//					exit(0);
//				}
//			}
//			cout << "\n";
			
//			unsigned long long hash = karp_rabin->hash(it_child->second->text.c_str(), pat_len);
			unsigned long long hash = karp_rabin->hash(test_text.c_str(), pat_len);
			hash_pat = karp_rabin->hash(pattern_rev.c_str() + pattern_rev.length() - pos + processed, pat_len);
//			cout << "HashTriev2RevNode::getRange - hash: " << hash << ", hash_pat: " << hash_pat << "\n";
			if( hash == hash_pat ){
//				cout << "HashTriev2RevNode::getRange - Child found -> [" << it_child->second->min << ", " << it_child->second->max << "]\n";
				return pair<unsigned int, unsigned int>(it_child->second->min, it_child->second->max);
			}
			
		}
	}
	*/
	
	
	
	
	
//	cout << "HashTriev2RevNode::getRange - Pattern NOT found\n";
	return pair<unsigned int, unsigned int>((unsigned int)(-1), (unsigned int)(-1));
}

void HashTriev2RevNode::save(fstream &writer){
	
	writer.write((char*)&len, sizeof(int));
	writer.write((char*)&min, sizeof(int));
	unsigned int max = 0;
	writer.write((char*)&max, sizeof(int));
	writer.write((char*)&hash, sizeof(int));
	unsigned int min_factor_pos = 0;
	writer.write((char*)&min_factor_pos, sizeof(int));
	
	unsigned int n_childs = childs_vector.size();
	writer.write((char*)&n_childs, sizeof(int));
	for( auto it : childs_vector ){
		char child_first_char = it.first;
		writer.write((char*)&child_first_char, 1);
		it.save(writer);
	}
	
}

void HashTriev2RevNode::load(fstream &reader){
	
	reader.read((char*)&len, sizeof(int));
	reader.read((char*)&min, sizeof(int));
	unsigned int max = 0;
	reader.read((char*)&max, sizeof(int));
	reader.read((char*)&hash, sizeof(int));
	// Lo dejo temporalmente por legacy
	unsigned int min_factor_pos = 0;
	reader.read((char*)&min_factor_pos, sizeof(int));
	
	unsigned int n_childs = 0;
	reader.read((char*)&n_childs, sizeof(int));
	for(unsigned int i = 0; i < n_childs; ++i){
		char child_first_char = 0;
		reader.read((char*)&child_first_char, 1);
		
//		childs[child_first_char] = std::make_shared<HashTriev2RevNode>();
//		childs[child_first_char]->load(reader);

		childs_vector.push_back(HashTriev2RevNode());
		childs_vector.back().load(reader);
		childs_vector.back().first = child_first_char;
	}
	
	std::sort(childs_vector.begin(), childs_vector.end(), 
		[](const HashTriev2RevNode &p1, const HashTriev2RevNode &p2) -> bool{return p1.first < p2.first;}
	);
	
}

void HashTriev2Rev::save(const string &file){
	cout << "HashTriev2Rev::save - Start (" << file << ")\n";
	fstream writer(file, fstream::out | fstream::trunc);
	root.save(writer);
	writer.close();
	cout << "HashTriev2Rev::save - End\n";
}
	
void HashTriev2Rev::load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, int_vector<> *_arr_x, const string &file){
	cout << "HashTriev2Rev::load - Start (" << file << ")\n";
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	arr_x = _arr_x;
	fstream reader(file, fstream::in);
	root.load(reader);
	reader.close();
	cout << "HashTriev2Rev::load - End\n";
}











