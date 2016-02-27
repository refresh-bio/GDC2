# Genome Differential Compressor 2
## GDC — What is it?

Genome Differential Compressor is a utility designed for compression of genome collections from the same species. The amount of such collections can be huge, e.g., a few (or tens) of gigabytes, so a need for a robust data compression tool is clear. Universal compression programs like gzip or bzip2 might be used for this purpose, but it is obvious that a specialized tool can work much better, since a universal compressor does not use the properties of such data sets, e.g., long approximate repetitions at long distances.

## GDC 2

### The architecture of GDC 2

GDC 2 is designed as a C++ application. The key features of the software are:
+ compression of collections of genomes in FASTA format,
+ decompression of the whole collection,
+ decompression of only a single genome without decompressing the complete collection,

### How good is GDC 2?

* Compression factor

In terms of compression factor (the ability to reduce the file size), GDC is usually much better than universal compressors and other specialized genome compressors. Its compression factor for some test data sets are:

∼9500 — on a collection of 1092 diploid genomes of H. sapiens from 1000 Genome Project,
∼590 — on a collection of 785 genomes of A. thaliana from 1001 Genome Project.

* Compression and decompression speed

The compression speed of GDC varies depending on data, but for the mentioned data sets is from 95 to 200MB/s at a computer equipped with 6-core Intel i7 4930K 3.4GHz processor. The decompression speed is dependent on the disk speed and is up to 1000 MB/s on the mentioned system.

### Publications

+ Deorowicz, S., Danek, A., Niemiec, M., GDC 2: Compression of large collections of genomes, Scientific Reports, 2015; 5(11565):1–12

### Developers

The GDC 2.x algorithm was invented by Sebastian Deorowicz and Agnieszka Danek. The implementation is by Sebastian Deorowicz, [Agnieszka Danek](https://github.com/agnieszkadanek), and Marcin Niemiec.
