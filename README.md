# align
A simple global alignment program for pairs of sequences

## Description
This repository contains a C++ program that performs sequence alignment using the Smith-Waterman algorithm. The program uses multi-threading with POSIX threads (pthreads) to perform calculations in parallel, significantly speeding up the alignment process.

It can handle both standard and gzipped FASTA files as input, thanks to the zlib library.

## Features
- Affine gap penalties with unique parameters for opening, extending, and ending gaps
- Multi-threading for initializing matrices and calculating alignment scores
- Extensible Fasta and Align classes to easily adapt or extend the program

## Dependencies
- C++ Standard Library
- POSIX Threads Library (pthreads)
- zlib

## How to Compile
To compile the program, you'll need to link against the pthreads and zlib libraries. 

Here is a sample command:
```
g++ -o align -lpthread -lz src/align.cc
```

## Usage
```
align test1.fa.gz test2.fa.gz id_table.txt
```

## Limitations
This program does not include any way to handle very large sequences that might not fit into memory. 
Nor does it include any kind of user interface or extensive error handling.
