# CSE305-Project-Parallel-sequence-alignment
The goal of the project is to implement two recent parallel algorithms (Needleman-Wunsch & Smith-Waterman) for sequence alignment and test their performance on real-life data.

### To Run the Needleman-Wunsch Algorithm:

```
g++ nw.cpp -o nw -std=c++11 -lpthread
./nw ../tests/big/midacin/H2QTE8.fasta ../tests/big/midacin/Q9NU22.fasta
```

where file1 and file2 are .fasta files that contain a sequence of amino acids and int is the number of threads you want to run the algorithm on.

AND 

### To Run the Smith-Waterman Algorithm:

To run the file, you must work only on the final_sw.cpp file (the again_sw.cpp file was used to test various functions and structure for the algorithm):

```
g++ Smith-Waterman/final_sw.cpp -o final_sw
./final_sw file1 file2 int
```

For example:

```
g++ Smith-Waterman/final_sw.cpp -o final_sw
./final_sw "tests/big/fibronectin/P02751.fasta" "tests/big/fibronectin/P07589.fasta" 64
``` 
### Other Files

Graphs.ipynb: to plot graphs to compare the different speeds for the algorithms.

sw.csv and nw.csv: Are the csv files that contain the results of running each of the commands from the command text files. The values in the csv files are: n (size of sequence 1), m (size of sequence 2), chunksize, runtime, threads
