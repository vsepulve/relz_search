# README #

Repository for the extension of compressed search based on relz_fs.

We use some structures from the SDSL (https://github.com/simongog/sdsl-lite).

install dependecies
-------------------
```
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
sudo ./install.sh /usr/local/
```

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


