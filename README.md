# Genome Differential Compressor
[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/gdc2/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/GDC2/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/gdc.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/gdc)

Genome Differential Compressor is a utility designed for compression of genome collections from the same species. The amount of such collections can be huge, e.g., a few (or tens) of gigabytes, so a need for a robust data compression tool is clear. Universal compression programs like gzip or bzip2 might be used for this purpose, but it is obvious that a specialized tool can work much better, since a universal compressor does not use the properties of such data sets, e.g., long approximate repetitions at long distances.

Thre two major releases of GDC: 2.1 and 0.3 that have different goals.
Version 0.3 was designed was compression of small collections of genomes, while GDC 2 was designed for larger collections.
Therefore, if you plan to compress a few genomes, GDC 0.x will likely produce smaller archive.
If you plan 100+ genomes, GDC 2 should be used.
GDC assumes that the genomes are given as sets of chromosomes.
In case you have them in sets of contigs, you should use <a href="https://github.com/refresh-bio/agc">AGC</a>.

## GDC 2
GDC 2 is designed as a C++ application. The key features of the software are:
+ compression of collections of genomes in FASTA format,
+ decompression of the whole collection,
+ decompression of only a single genome without decompressing the complete collection,

In terms of compression factor (the ability to reduce the file size), GDC is usually much better than universal compressors and other specialized genome compressors. Its compression factor for some test data sets are:
* ∼9500 &mdash; on a collection of 1092 diploid genomes of H. sapiens from 1000 Genome Project,
* ∼590 &mdash; on a collection of 785 genomes of A. thaliana from 1001 Genome Project.

The compression speed of GDC varies depending on data, but for the mentioned data sets is from 95 to 200MB/s at a computer equipped with 6-core Intel i7 4930K 3.4GHz processor. The decompression speed is dependent on the disk speed and is up to 1000 MB/s on the mentioned system.

## GDC 0.x
GDC 0.x is designed as a C++ library that can be used by various applications. The key features of the library are:
* compression of collections of genomes in FASTA format,
* decompression of the whole collection,
* decompression of only a single genome without decompressing the complete collection,
* decompression of any part of any chromosome from any genome of the collection without decompressing the complete collection and the complete genome (i.e., random * access functionality).

The package contains also two sample applications.

### GDC compressor
GDC is able to:
* compress collections of genomes in FASTA format,
* decompress whole collection,
* decompress only a single genome without decompressing the complete collection.

In terms of compression factor (the ability to reduce the file size), GDC is usually much better than universal compressors. Its compression factor for some test data sets are:
* 180&ndash;240 &mdash; on a collection of 70 genomes of H.sapiens from Complete Genomics Inc. (218,961.98 MB),
* 70&ndash;100 &mdash; on a collection of 39 genomes of S.cerevisiae (493.98 MB),
* 40&ndash;90 &mdash; on a collection of 36 genomes of S.paradoxus (436.43 MB),
* 12&ndash;17 &mdash; on a collaction of 4 genomes of H.sapiens (12,253.14 MB).

Since the idea of the GDC is to compress the genomes relatively to some reference genome, the more important compression factors could be defined as the ability of compact representation of the next genomes (the ones that follow the reference sequence). In such a case the compression factors for H. sapiens achieve from 420 to 1000. This means that a whole human genome can be stored in less than 3.12 MB.

The compression speed of GDC varies depending on data, but for the mentioned data sets is from 17 to 40MB/s at a computer equipped with AMD Opteron 2.4GHz processor. The decompression speed is about 150MB/s which is at least on a par with the I/O system speed.

### TEST_RA program
TEST_RA is an application that performs tests of the random access queries to the compressed archive.

The time necessary to decompress a snippet from the archive depends on the compression mode and the snippet size. However, for snippets of length 100 symbols the time is about 10&mu;s.


## Developers

GDC 2.x algorithm was invented by Sebastian Deorowicz and Agnieszka Danek. The implementation is by Sebastian Deorowicz, Agnieszka Danek, and Marcin Niemiec.

GDC 0.x algorithm was invented by Sebastian Deorowicz and Szymon Grabowski. The implementation is by Sebastian Deorowicz.

## Citing

[Deorowicz, S., Danek, A., Niemiec, M. (2015) GDC 2: Compression of large collections of genomes, Scientific Reports, 5:11565.](https://www.nature.com/articles/srep11565)

[Deorowicz, S., Grabowski, Sz. (2011) Robust relative compression of genomes with random access, Bioinformatics 27(21):2979&ndash;2986.](https://doi.org/10.1093/bioinformatics/btr505)

