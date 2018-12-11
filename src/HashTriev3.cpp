#include "HashTriev3.h"

HashTriev3::HashTriev3(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

HashTriev3::HashTriev3(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	arr_y = _arr_y;
	build(full_text, len_text, factors_start, arr_y, _karp_rabin, _kr_factors);
}

HashTriev3::~HashTriev3(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

void HashTriev3::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	arr_y = _arr_y;
	
	cout << "HashTriev3::build - Start (full text of " << len_text << ", " << factors_start.size() << " factors)\n";
	
	root.build(full_text, len_text, factors_start, arr_y, karp_rabin, kr_factors, 0, factors_start.size()-1, 0);
	
	// build en un nodo root local
	// Preparacion de arreglos para los datos agrupados
	// Iterar recursivamente por los nodos pasandole los arreglos
	// Desechar el root local
	
//	int_vector<> positions_childs;
//	int_vector<> n_childs;
//	int_vector<> len_childs;
//	int_vector<> min_childs;
//	int_vector<> hash_childs;
//	vector<char> first_childs

	unsigned int max_len = 0;
	unsigned int max_childs = 0;
	unsigned int max_height = 0;
	// including the root
	unsigned int n_nodes = 1 + root.totalChilds(max_len, max_childs, max_height, 0);
	
	cout << "HashTriev3::build - Total Nodes: " << n_nodes << " (max_len: " << max_len << ", max_childs: " << max_childs << ", max_height: " << max_height << ")\n";
	
	positions_childs.resize(n_nodes);
	n_childs.resize(n_nodes);
	len_childs.resize(n_nodes);
	min_childs.resize(n_nodes);
	hash_childs.resize(n_nodes);
	first_childs.resize(n_nodes);
	
	unsigned int next_pos = 0;
	
	// El llamador guarda la posicion de inicio de la raiz
	// Es decir, cada nodo guarda los datos de sus hijos
	positions_childs[next_pos] = 1;
	cout << "HashTriev3::build - positions_childs[" << next_pos << "]: 1\n";
	n_childs[next_pos] = root.childs_vector.size();
	len_childs[next_pos] = root.len;
	min_childs[next_pos] = root.min;
	hash_childs[next_pos] = root.hash;
	first_childs[next_pos] = root.first;
	++next_pos;
	
	root.compactData(next_pos, positions_childs, n_childs, len_childs, min_childs, hash_childs, first_childs);
	
	sdsl::util::bit_compress(positions_childs);
	sdsl::util::bit_compress(n_childs);
	sdsl::util::bit_compress(len_childs);
	sdsl::util::bit_compress(min_childs);
	sdsl::util::bit_compress(hash_childs);
//	sdsl::util::bit_compress(first_childs);
	
	// Debug
	for(unsigned int i = 0; i < n_nodes; ++i){
//		cout << "HashTriev3::build - positions_childs[" << i << "]: " << positions_childs[i] << " / n_childs[" << i << "]: " << n_childs[i] << "\n";
		cout << "HashTriev3::build - node[" << i << "]: (" << positions_childs[i] << ", " << n_childs[i] << ", " << len_childs[i] << ", " << min_childs[i] << ", " << hash_childs[i] << ", " << first_childs[i] << ")\n";
	}
	
	cout << "HashTriev3::build - End\n";
}

HashTriev3Node::HashTriev3Node(){
	len = 0;
	min = 0;
//	max = 0;
	hash = 0;
//	min_factor_pos = 0;
}

HashTriev3Node::~HashTriev3Node(){

}

void HashTriev3Node::compactData(unsigned int &next_pos, int_vector<> &positions_childs, int_vector<> &n_childs, int_vector<> &len_childs, int_vector<> &min_childs, int_vector<> &hash_childs, vector<char> &first_childs){
	
	// Guardo la posicion para los hijos de este nodo
	// El valor de next_pos cambiara con cada hijo
	unsigned int cur_pos = next_pos;
	next_pos += childs_vector.size();
	
//	cout << "HashTriev3Node::compactData - childs_vector.size(): " << childs_vector.size() << "\n";
	
	for(unsigned int i = 0; i < childs_vector.size(); ++i){
		
		positions_childs[cur_pos] = next_pos;
//		cout << "HashTriev3Node::compactData - positions_childs[" << cur_pos << "]: " << next_pos << "\n";
		n_childs[cur_pos] = childs_vector[i].childs_vector.size();
		len_childs[cur_pos] = childs_vector[i].len;
		min_childs[cur_pos] = childs_vector[i].min;
		hash_childs[cur_pos] = childs_vector[i].hash;
		first_childs[cur_pos] = childs_vector[i].first;
		
		childs_vector[i].compactData(next_pos, positions_childs, n_childs, len_childs, min_childs, hash_childs, first_childs);
		
		++cur_pos;
	}
	
}

void HashTriev3Node::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *arr_y, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len){
	if( min == max ){
//		cout << "HashTriev3Node::build - Leaf\n";
		return;
	}
	if( full_text == NULL || processed_len == len_text ){
//		cout << "HashTriev3Node::build - Null text\n";
		return;
	}
	
//	cout << "HashTriev3Node::build - Start (full text of " << len_text << ", range[" << min << ", " << max << "], processed_len: " << processed_len << ")\n";
	
	unsigned int min_pos = min;
	unsigned int min_text_start = factors_start[ (*arr_y)[min_pos] ] + processed_len;
	unsigned int min_text_len = 0;
	
	unsigned int cur_pos = 0;
	unsigned int cur_text_start = 0;
	
	unsigned long long hash;
	
	while( min_pos <= max ){
		
		// Primero busco el cur_pos MAYOR que tenga *primer* char igual al de min_pos
		// Despues reviso el largo mayor
		
//		cout << "HashTriev3Node::build - Iniciando busqueda binaria en [" << min_pos << ", " << max << "]\n";
		
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
//		cout << "HashTriev3Node::build - cur_pos: " << cur_pos << "\n";
		
		// Si cur_pos == min_pos, el largo es min_text_len
		// Si no, buscar el largo maximo entre min y cur
		
		if( cur_pos > min_pos ){
//			cout << "HashTriev3Node::build - Case 1\n";
			// Buscar el mayor
			// Aqui uso el hecho de que full_text termina en '\0'
			min_text_len = 1;
			while( full_text[min_text_start + min_text_len] == full_text[cur_text_start + min_text_len] ){
				++min_text_len;
			}
		}
		else{
//			cout << "HashTriev3Node::build - Case 2\n";
			min_text_len = len_text - min_text_start;
		}
//		cout << "HashTriev3Node::build - min_text_len: " << min_text_len << "\n";
//		cout << "HashTriev3Node::build - Preparing string s (min_text_start: " << min_text_start << " / " << len_text << ", len: " << min_text_len << ")\n";
		
		// El string es solo para el mensaje de debug
//		string s(full_text + min_text_start, min_text_len);
		char first_char = full_text[min_text_start];
		
		// Notar que aqui necesito la posicion del factor en la coleccion, no en el arreglo Y
		hash = kr_factors->hash( (*arr_y)[min_pos], processed_len, min_text_len);
//		cout << "HashTriev3Node::build - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_text_len << ", hash: " << hash << " / " << karp_rabin->hash(s) << " (\"" << s << "\")\n";
//		cout << "HashTriev3Node::build - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_text_len << ", hash: " << hash << "\n";
		// Preparar el hijo, ejecutar la llamda sobre esa instancia
		if( min_text_len == 0 ){
//			cout << "HashTriev3Node::build - Omiting child of len 0\n";
		}
		else{
			childs_vector.push_back(HashTriev3Node());
			childs_vector.back().first = first_char;
			childs_vector.back().len = min_text_len;
			childs_vector.back().min = min_pos;
			childs_vector.back().hash = hash;
			childs_vector.back().build(full_text, len_text, factors_start, arr_y, karp_rabin, kr_factors, min_pos, cur_pos, processed_len + min_text_len);
		}
		
//		cout << "HashTriev3Node::build - Preparing min_pos = " << cur_pos+1 << " / " << arr_y->size() << "\n";
		min_pos = cur_pos+1;
		if(min_pos >= arr_y->size()){
			break;
		}
		min_text_start = factors_start[ (*arr_y)[min_pos] ] + processed_len;
		
//		cout << "HashTriev3Node::build - min_text_start: " << min_text_start << "\n";
	
	}
	
	std::sort(childs_vector.begin(), childs_vector.end(), 
		[](const HashTriev3Node &p1, const HashTriev3Node &p2) -> bool{return p1.first < p2.first;}
	);
	
//	cout << "----- \n";
	
}

void HashTriev3::print(){
	root.print(0);
	cout << "Root: " << root.childs_vector.size() << " childs\n";
	cout << "Lengths: ";
	for( auto it_child : root.childs_vector ){
		cout << it_child.len << " - ";
	}
	cout << "\n";
}

void HashTriev3Node::print(unsigned int level){
	for(unsigned int i = 0; i < level; ++i){
		cout << "- ";
	}
	cout << "hash: " << hash << ", len: " << len << " , range [" << min << "]\n";
	for( auto it : childs_vector ){
		cout << it.first << " ";
		it.print(level+1);
	}
}

unsigned int HashTriev3Node::totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height){
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

void HashTriev3::printSize(){
	unsigned int max_len = 0;
	unsigned int max_childs = 0;
	unsigned int max_height = 0;
	unsigned int total_childs = root.totalChilds(max_len, max_childs, max_height, 0);
	cout << "HashTriev3::printSize - totalChilds: " << total_childs << ", max_len: " << max_len << ", max_childs: " << max_childs << ", max_height: " << max_height << "\n";
}

pair<unsigned int, unsigned int> HashTriev3::getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, const string &pattern){
//	return root.getRange(kr_pat_vector, pos, 0, karp_rabin, kr_factors, arr_y, pattern, &hash_nano);
	return root.getRange(kr_pat_vector, pos, 0, karp_rabin, kr_factors, arr_y, arr_y->size() - 1, pattern, &hash_nano);
}

unsigned int HashTriev3Node::findChild(char c){
	
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

pair<unsigned int, unsigned int> HashTriev3Node::getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_y, unsigned int cur_max, const string &pattern, unsigned long long *hash_nano){

//pair<unsigned int, unsigned int> HashTriev3Node::getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_y, const string &pattern, unsigned long long *hash_nano){
	
//	cout << "HashTriev3Node::getRange - Start (prefixes: " << kr_pat_vector.size() << ", pos: " << pos << ", processed: " << processed << ")\n";
	
	if( pos + processed >= kr_pat_vector.size() ){
//		cout << "HashTriev3Node::getRange - [" << min << ", " << cur_max << "]\n";
//		return pair<unsigned int, unsigned int>(min, max);
		return pair<unsigned int, unsigned int>(min, cur_max);
	}
	
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	unsigned int pat_len = kr_pat_vector.size() - pos - processed;
	char first_char_pat = pattern[pos + processed];
//	cout << "HashTriev3Node::getRange - pat_len: " << pat_len << "\n";
	string pat = pattern.substr(pos + processed, pat_len);
//	cout << "HashTriev3Node::getRange - pat: " << pat << "\n";
	
	unsigned int pos_child = findChild(first_char_pat);
	if( pos_child != NOT_FOUND ){
		child_len = childs_vector[pos_child].len;
		
		// Ajuste a cur_max
		if( pos_child < childs_vector.size() - 1 ){
			cur_max = childs_vector[pos_child+1].min - 1;
		}
		
		if( child_len <= pat_len ){
//			cout << "HashTriev3Node::getRange - Case 1, child_len: " << child_len << "\n";
			hash_pat = karp_rabin->subtract_prefix(kr_pat_vector[pos + processed + child_len - 1], kr_pat_vector[pos + processed - 1], child_len);
			string pat_cut = pattern.substr(pos + processed, child_len);
//			cout << "HashTriev3Node::getRange - pat_cut: " << pat_cut << " hash_pat: " << hash_pat << " / " << karp_rabin->hash(pat_cut) << " (processed: " << processed << ")\n";
			if( hash_pat == childs_vector[pos_child].hash ){
//				cout << "HashTriev3Node::getRange - Child found -> [" << childs_vector[pos_child].min << ", " << cur_max << "]\n";
//				return childs_vector[pos_child].getRange(kr_pat_vector, pos, processed + child_len, karp_rabin, kr_factors, arr_y, pattern, hash_nano);
				return childs_vector[pos_child].getRange(kr_pat_vector, pos, processed + child_len, karp_rabin, kr_factors, arr_y, cur_max, pattern, hash_nano);
			}
		}
		else{
//			cout << "HashTriev3Node::getRange - Case 2, child_len: " << child_len << "\n";
			
			unsigned int min_factor_pos = (*arr_y)[childs_vector[pos_child].min];
			
//			cout << "Node: \n";
//			childs_vector[pos_child].print(0);
			// Caso de borde detectado a veces en kr_factors->hashFast(childs_vector[pos_child].min_factor_pos, processed, pat_len
			
			unsigned long long hash = kr_factors->hashFast(min_factor_pos, processed, pat_len);
			
			hash_pat = karp_rabin->subtract_prefix(kr_pat_vector[kr_pat_vector.size() - 1], kr_pat_vector[pos + processed - 1], kr_pat_vector.size() - pos - processed);
//			cout << "HashTriev3Node::getRange - hash: " << hash << ", hash_pat: " << hash_pat << "\n";
			if( hash == hash_pat ){
//				cout << "HashTriev3Node::getRange - Child found -> [" << childs_vector[pos_child].min << ", " << cur_max << "]\n";
//				return pair<unsigned int, unsigned int>(childs_vector[pos_child].min, childs_vector[pos_child].max);
				return pair<unsigned int, unsigned int>(childs_vector[pos_child].min, cur_max);
			}
		}
	}
	
//	cout << "HashTriev3Node::getRange - Pattern NOT found\n";
	return pair<unsigned int, unsigned int>((unsigned int)(-1), (unsigned int)(-1));
}

void HashTriev3Node::save(fstream &writer){
	
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

void HashTriev3Node::load(fstream &reader){
	
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
		
//		childs[child_first_char] = std::make_shared<HashTriev3Node>();
//		childs[child_first_char]->load(reader);
		
		childs_vector.push_back(HashTriev3Node());
		childs_vector.back().load(reader);
		childs_vector.back().first = child_first_char;
		
	}
	
	std::sort(childs_vector.begin(), childs_vector.end(), 
		[](const HashTriev3Node &p1, const HashTriev3Node &p2) -> bool{return p1.first < p2.first;}
	);
	
	// Debug
//	cout << "HashTriev3Node::load - Childs: " << childs_vector.size() << " / ";
//	for(auto it : childs_vector){
//		cout << it.first << " ";
//	}
//	cout << "\n";

}

void HashTriev3::save(const string &file){
	cout << "HashTriev3::save - Start (" << file << ")\n";
	fstream writer(file, fstream::out | fstream::trunc);
	root.save(writer);
	writer.close();
	cout << "HashTriev3::save - End\n";
}

void HashTriev3::load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, int_vector<> *_arr_y, const string &file){
	cout << "HashTriev3::load - Start (" << file << ")\n";
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	arr_y = _arr_y;
	fstream reader(file, fstream::in);
	root.load(reader);
	reader.close();
	cout << "HashTriev3::load - End\n";
}
	
// VERSIONES REV

HashTriev3Rev::HashTriev3Rev(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

HashTriev3Rev::HashTriev3Rev(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	arr_x = _arr_x;
	build(full_text, len_text, factors_start, arr_x, _karp_rabin, _kr_factors);
}

HashTriev3Rev::~HashTriev3Rev(){
	karp_rabin = NULL;
	kr_factors = NULL;
}

void HashTriev3Rev::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors){
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	arr_x = _arr_x;
	
	cout << "HashTriev3Rev::build - Start (full text of " << len_text << ", " << factors_start.size() << " factors)\n";
	
	root.build(full_text, len_text, factors_start, arr_x, karp_rabin, kr_factors, 0, factors_start.size()-1, 0);
	
	cout << "HashTriev3Rev::build - End\n";
}

HashTriev3RevNode::HashTriev3RevNode(){
	len = 0;
	min = 0;
//	max = 0;
	hash = 0;
//	min_factor_pos = 0;
}

HashTriev3RevNode::~HashTriev3RevNode(){

}

void HashTriev3RevNode::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *arr_x, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len){
	if( min == max ){
//		cout << "HashTriev3RevNode::build - Leaf\n";
		return;
	}
	if( full_text == NULL || processed_len == len_text ){
//		cout << "HashTriev3RevNode::build - Null text\n";
		return;
	}
	
//	cout << "HashTriev3RevNode::build - Start (full text of " << len_text << ", range[" << min << ", " << max << "], processed_len: " << processed_len << ")\n";
	
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
		
//		cout << "HashTriev3RevNode::build - Iniciando busqueda binaria en [" << min_pos << ", " << max << "]\n";
		
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
//		cout << "HashTriev3RevNode::build - cur_pos: " << cur_pos << "\n";
		
		// Si cur_pos == min_pos, el largo es min_text_len
		// Si no, buscar el largo maximo entre min y cur
		
		if( cur_pos > min_pos ){
//			cout << "HashTriev3RevNode::build - Case 1\n";
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
//			cout << "HashTriev3RevNode::build - Case 2\n";
			min_common_text = min_text_len;
		}
//		cout << "HashTriev3RevNode::build - min_common_text: " << min_common_text << "\n";
//		cout << "HashTriev3RevNode::build - Preparing string s (min_text_start: " << min_text_start << " / " << len_text << ", len: " << min_common_text << ")\n";
		
		// En el caso de Rev, uso s para el hash directo
		string s;
		for( unsigned int i = 0; i < min_common_text; ++i ){
			s += *(full_text + min_text_start - i);
		}
		hash = karp_rabin->hash(s);
//		cout << "HashTriev3RevNode::build - s: " << s << "\n";
		char first_char = s[0];
//		cout << "HashTriev3RevNode::build - hash: " << hash << "\n";
		
//		cout << "HashTriev3RevNode::build - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_common_text << ", hash: " << hash << "\n";
		// Preparar el hijo, ejecutar la llamda sobre esa instancia
		if( min_common_text == 0 ){
//			cout << "HashTriev3RevNode::build - Omiting child of len 0\n";
		}
		else{
			childs_vector.push_back(HashTriev3RevNode());
			childs_vector.back().first = first_char;
			childs_vector.back().len = min_common_text;
			childs_vector.back().min = min_pos;
			childs_vector.back().hash = hash;
			childs_vector.back().build(full_text, len_text, factors_start, arr_x, karp_rabin, kr_factors, min_pos, cur_pos, processed_len + min_common_text);
		}
		
//		cout << "HashTriev3RevNode::build - Preparing min_pos = " << cur_pos+1 << " / " << arr_x->size() << "\n";
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
//		cout << "HashTriev3RevNode::build - min_text_start: " << min_text_start << " de len " << min_text_len << "\n";
	
	}
	
	std::sort(childs_vector.begin(), childs_vector.end(), 
		[](const HashTriev3RevNode &p1, const HashTriev3RevNode &p2) -> bool{return p1.first < p2.first;}
	);
	
//	cout << "----- \n";
	
}


void HashTriev3Rev::print(){
	root.print(0);
	cout << "Root: " << root.childs_vector.size() << " childs\n";
	cout << "Lengths: ";
	for( auto it_child : root.childs_vector ){
		cout << it_child.len << " - ";
	}
	cout << "\n";
}

void HashTriev3RevNode::print(unsigned int level){
	for(unsigned int i = 0; i < level; ++i){
		cout << "- ";
	}
//	cout << "hash: " << hash << ", len: " << len << ", factor_base " << factor << " , range [" << min << "]\n";
	cout << "hash: " << hash << ", len: " << len << " , range [" << min << "]\n";
	for( auto it : childs_vector ){
		cout << it.first << " ";
		it.print(level+1);
	}
}

unsigned int HashTriev3RevNode::totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height){
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

void HashTriev3Rev::printSize(){
	unsigned int max_len = 0;
	unsigned int max_childs = 0;
	unsigned int max_height = 0;
	unsigned int total_childs = root.totalChilds(max_len, max_childs, max_height, 0);
	cout << "HashTriev3Rev::printSize - totalChilds: " << total_childs << ", max_len: " << max_len << ", max_childs: " << max_childs << ", max_height: " << max_height << "\n";
}

pair<unsigned int, unsigned int> HashTriev3Rev::getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, const string &pattern_rev){
	return root.getRange(kr_pat_rev_vector, pos, 0, karp_rabin, kr_factors, arr_x, arr_x->size() - 1, pattern_rev);
}

unsigned int HashTriev3RevNode::findChild(char c){
	
//	if(childs_vector.size() < 1){
//		return NOT_FOUND;
//	}
//	// Busqueda Binaria
	
//	cout << "HashTriev3RevNode::findChild - Char " << c << ", " << childs_vector.size() << " childs ( ";
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

pair<unsigned int, unsigned int> HashTriev3RevNode::getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_x, unsigned int cur_max, const string &pattern_rev){
	
//	cout << "HashTriev3RevNode::getRange - Start (prefixes: " << kr_pat_rev_vector.size() << ", pos: " << pos << ", processed: " << processed << ", pattern_rev: " << pattern_rev << ")\n";
	
	if( processed >= pos ){
//		cout << "HashTriev3RevNode::getRange - [" << min << ", " << cur_max << "]\n";
		return pair<unsigned int, unsigned int>(min, cur_max);
	}
	
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	unsigned int pat_len = pos - processed;
	char first_char_pat = pattern_rev[kr_pat_rev_vector.size() - 1 - pos + processed];
//	cout << "HashTriev3RevNode::getRange - pat = pattern_rev.substr(" << (kr_pat_rev_vector.size() - 1 - pos + processed) << ", " << pat_len << ");\n";
	string pat = pattern_rev.substr(kr_pat_rev_vector.size() - 1 - pos + processed, pat_len);
//	cout << "HashTriev3RevNode::getRange - pat: " << pat << ", first_char_pat: " << first_char_pat << "\n";
	
	unsigned int pos_child = findChild(first_char_pat);
	if( pos_child != NOT_FOUND ){
		child_len = childs_vector[pos_child].len;
		
		// Ajuste a cur_max
//		cout << "HashTriev3RevNode::getRange - Adjusting cur_max (pos_child: " << pos_child << " / " << childs_vector.size() << ")\n";
		if( pos_child < childs_vector.size() - 1 ){
			cur_max = childs_vector[pos_child+1].min - 1;
//			cout << "HashTriev3RevNode::getRange - cur_max: " << cur_max << ")\n";
		}
		
		if( child_len <= pat_len ){
//			cout << "HashTriev3RevNode::getRange - Case 1, child_len: " << child_len << "\n";
			hash_pat = karp_rabin->hash(pattern_rev.c_str() + pattern_rev.length() - pos + processed, child_len);
			string pat_cut = pattern_rev.substr(kr_pat_rev_vector.size() - 1 - pos + processed, child_len);
//			cout << "HashTriev3RevNode::getRange - pat_cut: " << pat_cut << " hash_pat: " << hash_pat << " / " << karp_rabin->hash(pat_cut) << " (processed: " << processed << ")\n";
			if( hash_pat == childs_vector[pos_child].hash ){
//				cout << "HashTriev3RevNode::getRange - Child found -> [" << childs_vector[pos_child].min << ", " << cur_max << "]\n";
				return childs_vector[pos_child].getRange(kr_pat_rev_vector, pos, processed + child_len, karp_rabin, kr_factors, arr_x, cur_max, pattern_rev);
			}
		}
		else{
//			cout << "HashTriev3RevNode::getRange - Case 2, child_len: " << child_len << "\n";
			
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
//					cout << "HashTriev3RevNode::getRange - Adding " << len << " chars\n";
					for(unsigned int i = 0; i < len; ++i){
						test_text += *(ptr->ref_text + tu + lu - processed - i - 1);
					}
				}
			}
			
			unsigned long long hash = karp_rabin->hash(test_text.c_str(), pat_len);
			hash_pat = karp_rabin->hash(pattern_rev.c_str() + pattern_rev.length() - pos + processed, pat_len);
//			cout << "HashTriev3RevNode::getRange - hash: " << hash << ", hash_pat: " << hash_pat << "\n";
			if( hash == hash_pat ){
//				cout << "HashTriev3RevNode::getRange - Child found -> [" << childs_vector[pos_child].min << ", " << cur_max << "]\n";
				return pair<unsigned int, unsigned int>(childs_vector[pos_child].min, cur_max);
			}
			
		}
	}
	
//	cout << "HashTriev3RevNode::getRange - Pattern NOT found\n";
	return pair<unsigned int, unsigned int>((unsigned int)(-1), (unsigned int)(-1));
}

void HashTriev3RevNode::save(fstream &writer){
	
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

void HashTriev3RevNode::load(fstream &reader){
	
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
		
//		childs[child_first_char] = std::make_shared<HashTriev3RevNode>();
//		childs[child_first_char]->load(reader);

		childs_vector.push_back(HashTriev3RevNode());
		childs_vector.back().load(reader);
		childs_vector.back().first = child_first_char;
	}
	
	std::sort(childs_vector.begin(), childs_vector.end(), 
		[](const HashTriev3RevNode &p1, const HashTriev3RevNode &p2) -> bool{return p1.first < p2.first;}
	);
	
}

void HashTriev3Rev::save(const string &file){
	cout << "HashTriev3Rev::save - Start (" << file << ")\n";
	fstream writer(file, fstream::out | fstream::trunc);
	root.save(writer);
	writer.close();
	cout << "HashTriev3Rev::save - End\n";
}
	
void HashTriev3Rev::load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, int_vector<> *_arr_x, const string &file){
	cout << "HashTriev3Rev::load - Start (" << file << ")\n";
	karp_rabin = _karp_rabin;
	kr_factors = _kr_factors;
	arr_x = _arr_x;
	fstream reader(file, fstream::in);
	root.load(reader);
	reader.close();
	cout << "HashTriev3Rev::load - End\n";
}











