# README #

Repository for the extension of compressed search based on relz_fs.
This version uses the code of https://github.com/hferrada/rmq [3], [4].
 

The main programs of this version of the code are compress, decompress.

Test
-------------------
```
mkdir build
cd build
cmake ..
make
./bin/compress ../data/ref.txt ../data/seq.txt ../data/seq.relz 100000 4 ../data/ref.bin 1 1
./bin/decompress ../data/ref.bin ../data/seq.relz ../data/seq_decomp.txt 100000

```

References:

[1]. J. Fischer and V. Heun. Space-efficient preprocessing schemes for range minimum 
queries on static arrays. SIAM '11.

[2]. K. Sadakane and G. Navarro. Fully-functional succinct trees. In Proceedings of 
the Twenty-First Annual ACM-SIAM Symposium on Discrete Algorithms, SODA '10.

[3]. H. Ferrada and G. Navarro. Improved Range Minimum Queries. In 
Journal of Discrete Algorithms, JDA'17.

[4]. H. Ferrada and G. Navarro. Improved Range Minimum Queries. In Proceedings of 
Data Compression Conference, DCC'16.
