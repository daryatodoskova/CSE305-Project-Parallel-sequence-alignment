#include <iostream>
#include <fstream>
#include <vector>
#include <thread>

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
void smithWaterman(const std::string& sequence1, const std::string& sequence2,
                   int matchScore, int mismatchScore, int gapPenalty,
                   std::vector<std::vector<int>>& score, int startRow, int endRow) {
    // Get the lengths of the sequences
    int m = sequence1.length();
    int n = sequence2.length();

    // Iterate over each cell in the specified range of rows
    for (int i = startRow; i <= endRow; ++i) {
        for (int j = 1; j <= n; ++j) {
            // Calculate the score for the current cell
            int scoreDiagonal = score[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? matchScore : mismatchScore);
            int scoreUp = score[i - 1][j] - gapPenalty;
            int scoreLeft = score[i][j - 1] - gapPenalty;

            // Take the maximum score
            int maxScore = maximum(0, scoreDiagonal, scoreUp, scoreLeft);
            score[i][j] = maxScore;
        }
    }
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

int main() {
    // Read the sequences from a file
    std::ifstream inputFile("sequences.txt");
    std::string sequence1, sequence2;
    if (inputFile.is_open()) {
        std::getline(inputFile, sequence1);
        std::getline(inputFile, sequence2);
        inputFile.close();
    } else {
        std::cout << "Failed to open the input file." << std::endl;
        return 1;
    }

    // Specify the scoring scheme and other parameters
    int matchScore = 2;
    int mismatchScore = -1;
    int gapPenalty = 1;

    // Get the lengths of the sequences
    int m = sequence1.length();
    int n = sequence2.length();

    // Initialize the score matrix
    std::vector<std::vector<int>> score = initializeMatrix(m + 1, n + 1);

    // Manually specify the number of threads
    unsigned int numThreads = 4;

    // Divide the rows evenly among the threads
    unsigned int rowsPerThread = m / numThreads;

    // Create a vector to store the threads
    std::vector<std::thread> threads;

    // Perform the Smith-Waterman algorithm in parallel
    for (unsigned int t = 0; t < numThreads; ++t) {
        int startRow = t * rowsPerThread + 1;
        int endRow = (t == numThreads - 1) ? m : (t + 1) * rowsPerThread;
        threads.emplace_back(smithWaterman, std::ref(sequence1), std::ref(sequence2),
                             matchScore, mismatchScore, gapPenalty,
                             std::ref(score), startRow, endRow);
    }

    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
    }

    // Find the optimal local alignment
    std::pair<std::string, std::string> alignment = traceback(sequence1, sequence2, score, matchScore, mismatchScore, gapPenalty);

    // Print the alignment
    std::cout << "Sequence 1: " << alignment.first << std::endl;
    std::cout << "Sequence 2: " << alignment.second << std::endl;

    return 0;
}
