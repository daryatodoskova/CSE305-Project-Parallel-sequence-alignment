# CSE305-Project-Parallel-sequence-alignment
The goal of the project will be to implement one of the recent parallel algorithms for sequence alignment and test its performance on real-life data.


Mid Project Report:

We started by analyzing the different algorithms we could implement and chose the Needleman-Wunsch algorithm. We might implement the Gotoh algorithm or Smith-Waterman too if we have time. We then thoroughly studied its main data structures (for ex the scoring matrix) and steps to compute the alignment score and traceback path.

We then focused on which parallelization strategy we were going to use. We decided to divide the problem into smaller subproblems ,by dividing the list of indices (representing diagonal elements of the matrix) into blocks. We then chose how we wanted to distribute the work across multiple threads (we assigned each block to a separate thread) and how to handle dependencies and communication between different parts of the algorithm (choosing which parts of the algorithm could be executed concurrently).

Elsa and I are currently looking for real-life datasets from bioinformatics databases (.fasta format) that we could use for testing and validation. We are looking for some of different lengths (big vs small) and different characteristics. Once we’re done, we’ll compare the results and execution times with a sequential implementation to check how much faster it got after setting up the parallelization.

In the meantime, Jean-Sébastien is comparing the Gotoh and Smith-Waterman (and potentially another one that we might prefer) to start implementing it as we did with NW.