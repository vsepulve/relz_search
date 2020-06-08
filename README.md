# README #

Repository for the extension of compressed search based on RLZ.

We use some structures from the SDSL (https://github.com/simongog/sdsl-lite) and from Google Sparsehash (https://github.com/sparsehash/sparsehash).

install dependecies
-------------------
SDSL
```
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite/
sudo ./install.sh /usr/local/
cd ..
```
Sparsehash
```
git clone https://github.com/sparsehash/sparsehash
cd sparsehash/
./configure
make
sudo make install
```

Compile
-------------------
```
mkdir build
cd build
cmake ..
make
```

Building (indexing) the Reference
-------------------
The first step to use the system is indexing of the reference text to generate the binary references used for the rest of the processes. The program bin/build\_reference performs the indexing, and its arguments are the reference text file, the output to write the binary and indexed reference, and the number of threads used for construction.
An example could be:
```
./bin/build_reference reference_text reference_output n_threads
```

Data Preparation
-------------------
```
cd ../data
tar -zxvf yeast_ref.tgz
tar -zxvf yeast_queries.tgz
cd ../build
./bin/build_reference ../data/yeast_ref.txt ../data/yeast_ref.bin 4
./bin/decompress ../data/yeast_ref.bin ../data/yeast_data.relz yeast_data.txt 1000000
```

Indexes Construction (each process can take some time)
-------------------
```
mkdir index_yeast
./bin/build_relz_basic ../data/yeast_ref.bin yeast_data.txt index_yeast/
./bin/build_relz_reference ../data/yeast_ref.bin yeast_data.txt index_yeast/
./bin/build_relz_hash ../data/yeast_ref.bin yeast_data.txt index_yeast/
```


Testing Basic Version
-------------------
```
./bin/relz_search_basic index_yeast/ ../data/yeast_queries_10k_10.txt
./bin/relz_search_basic index_yeast/ ../data/yeast_queries_10k_20.txt
./bin/relz_search_basic index_yeast/ ../data/yeast_queries_10k_40.txt
./bin/relz_search_basic index_yeast/ ../data/yeast_queries_10k_80.txt
```

Testing Reference Version
-------------------
```
./bin/relz_search_reference index_yeast/ ../data/yeast_queries_10k_10.txt
./bin/relz_search_reference index_yeast/ ../data/yeast_queries_10k_20.txt
./bin/relz_search_reference index_yeast/ ../data/yeast_queries_10k_40.txt
./bin/relz_search_reference index_yeast/ ../data/yeast_queries_10k_80.txt
```

Testing Hash Version
-------------------
```
./bin/relz_search_hash index_yeast/ ../data/yeast_queries_10k_10.txt
./bin/relz_search_hash index_yeast/ ../data/yeast_queries_10k_20.txt
./bin/relz_search_hash index_yeast/ ../data/yeast_queries_10k_40.txt
./bin/relz_search_hash index_yeast/ ../data/yeast_queries_10k_80.txt
```

Testing Hash with Global Table Version
-------------------
```
./bin/relz_search_hash_table index_yeast/ ../data/yeast_queries_10k_10.txt
./bin/relz_search_hash_table index_yeast/ ../data/yeast_queries_10k_20.txt
./bin/relz_search_hash_table index_yeast/ ../data/yeast_queries_10k_40.txt
./bin/relz_search_hash_table index_yeast/ ../data/yeast_queries_10k_80.txt
```

Testing No Hash Version
-------------------
This version uses the same index of Hash, but omits the usage of the rolling hash. Instead, it receives approximate (greater) ranges from the trees and filter the results with the correctness test before reporting them.
```
./bin/relz_search_hash_simple index_yeast/ ../data/yeast_queries_10k_10.txt
./bin/relz_search_hash_simple index_yeast/ ../data/yeast_queries_10k_20.txt
./bin/relz_search_hash_simple index_yeast/ ../data/yeast_queries_10k_40.txt
./bin/relz_search_hash_simple index_yeast/ ../data/yeast_queries_10k_80.txt
```




