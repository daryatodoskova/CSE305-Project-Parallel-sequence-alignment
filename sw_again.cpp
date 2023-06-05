#include <iostream>
#include <vector>
#include <omp.h>
#include <fstream>
#include <string>

// Function to initialize the matrix
std::vector<std::vector<int>> initializeMatrix(int rows, int cols) {
    std::vector<std::vector<int>> matrix(rows, std::vector<int>(cols, 0));
    return matrix;
}

// Function to calculate the maximum of three values
int maximum(int a, int b, int c) {
    return std::max(std::max(a, b), c);
}

// Function to perform the Smith-Waterman algorithm
std::vector<std::vector<int>> smithWaterman(const std::string& sequence1, const std::string& sequence2,
                                            int matchScore, int mismatchScore, int gapPenalty) {
    // Get the lengths of the sequences
    int m = sequence1.length();
    int n = sequence2.length();

    // Initialize the matrix
    std::vector<std::vector<int>> score(m + 1, std::vector<int>(n + 1, 0));

    // Set the number of threads for parallelization
    int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);

    // Iterate over each diagonal in the matrix
    for (int d = 1 - m; d < n; ++d) {
        #pragma omp parallel for
        for (int i = std::max(1, d + 1); i <= std::min(m, d + n); ++i) {
            int j = i - d;

            // Calculate the score for the current cell
            int scoreDiagonal = score[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? matchScore : mismatchScore);
            int scoreUp = score[i - 1][j] - gapPenalty;
            int scoreLeft = score[i][j - 1] - gapPenalty;

            // Take the maximum score
            int maxScore = maximum(0, scoreDiagonal, scoreUp, scoreLeft);
            score[i][j] = maxScore;
        }
    }

    return score;
}

// Function to find the cell with the maximum score
std::pair<int, int> findMaxScoreCell(const std::vector<std::vector<int>>& score) {
    int maxScore = 0;
    int maxScoreRow = 0;
    int maxScoreCol = 0;

    for (int i = 0; i < score.size(); ++i) {
        for (int j = 0; j < score[i].size(); ++j) {
            if (score[i][j] > maxScore) {
                maxScore = score[i][j];
                maxScoreRow = i;
                maxScoreCol = j;
            }
        }
    }

    return std::make_pair(maxScoreRow, maxScoreCol);
}

// Function to perform traceback and find the optimal local alignment
std::pair<std::string, std::string> traceback(const std::string& sequence1, const std::string& sequence2,
                                              const std::vector<std::vector<int>>& score,
                                              int matchScore, int mismatchScore, int gapPenalty) {
    // Find the cell with the maximum score
    std::pair<int, int> maxScoreCell = findMaxScoreCell(score);
    int i = maxScoreCell.first;
    int j = maxScoreCell.second;

    // Traceback until reaching a cell with a score of 0
    std::string alignedSequence1, alignedSequence2;
    while (score[i][j] != 0) {
        if (score[i][j] == score[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? matchScore : mismatchScore)) {
            alignedSequence1 = sequence1[i - 1] + alignedSequence1;
            alignedSequence2 = sequence2[j - 1] + alignedSequence2;
            --i;
            --j;
        } else if (score[i][j] == score[i - 1][j] - gapPenalty) {
            alignedSequence1 = sequence1[i - 1] + alignedSequence1;
            alignedSequence2 = '-' + alignedSequence2;
            --i;
        } else {
            alignedSequence1 = '-' + alignedSequence1;
            alignedSequence2 = sequence2[j - 1] + alignedSequence2;
            --j;
        }
    }

    return std::make_pair(alignedSequence1, alignedSequence2);
}


// Function to read sequences from a .fasta file
std::vector<std::string> readSequencesFromFile(const std::string& filename) {
    std::vector<std::string> sequences;
    std::ifstream file(filename);

    if (file.is_open()) {
        std::string line;
        std::string sequence;
        while (std::getline(file, line)) {
            if (!line.empty() && line[0] != '>') {
                sequence += line;
            }
        }

        sequences.push_back(sequence);
        file.close();
    }

    return sequences;
}

int main() {
    std::string filename = "sequences.fasta";
    std::vector<std::string> sequences = readSequencesFromFile(filename);

    if (sequences.size() < 2) {
        std::cout << "Error: Insufficient number of sequences in the file.\n";
        return 1;
    }

    std::string sequence1 = sequences[0];
    std::string sequence2 = sequences[1];

    int matchScore = 2;
    int mismatchScore = -1;
    int gapPenalty = -2;

    // Perform Smith-Waterman algorithm
    std::vector<std::vector<int>> scoreMatrix = smithWaterman(sequence1, sequence2, matchScore, mismatchScore, gapPenalty);

    // Print the score matrix
    for (int i = 0; i < scoreMatrix.size(); ++i) {
        for (int j = 0; j < scoreMatrix[i].size(); ++j) {
            std::cout << scoreMatrix[i][j] << '\t';
        }
        std::cout << '\n';
    }

    // Perform traceback and print the optimal alignment
    std::pair<std::string, std::string> alignment = traceback(sequence1, sequence2, scoreMatrix, matchScore, mismatchScore, gapPenalty);
    std::cout << "Optimal alignment:\n";
    std::cout << alignment.first << '\n';
    std::cout << alignment.second << '\n';

    return 0;
}

