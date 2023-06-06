# CSE305-Project-Parallel-sequence-alignment
The goal of the project is to implement two recent parallel algorithms (Needleman-Wunsch & Smith-Waterman) for sequence alignment and test their performance on real-life data.

To compile on mac (example):

```
g++ nw.cpp -o nw -std=c++11 -lpthread
./nw tests/big/midacin/H2QTE8.fasta tests/big/midacin/Q9NU22.fasta
```

AND 

