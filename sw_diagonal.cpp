#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>

const int MATCH_SCORE = 2;
const int MISMATCH_SCORE = -1;
const int GAP_SCORE = -2;

// Define the sequence types
using Sequence = std::vector<char>;
using ScoreMatrix = std::vector<std::vector<int>>;

// Function to initialize the score matrix with zeros
void initializeMatrix(ScoreMatrix& matrix, int rows, int cols) {
    matrix.resize(rows, std::vector<int>(cols, 0));
}


// Function to calculate the maximum score and its position
void findMaxScore(const ScoreMatrix& matrix, int& maxScore, int& maxI, int& maxJ) {
    maxScore = 0;
    maxI = 0;
    maxJ = 0;
  
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix[i].size(); ++j) {
            if (matrix[i][j] > maxScore) {
                maxScore = matrix[i][j];
                maxI = i;
                maxJ = j;
            }
        }
    }
}

// Function to perform sequence alignment using the Smith-Waterman algorithm
void smithWaterman(const Sequence& sequence1, const Sequence& sequence2, ScoreMatrix& matrix) {
    const int rows = sequence1.size() + 1;
    const int cols = sequence2.size() + 1;

    initializeMatrix(matrix, rows, cols);

    for (int i = 1; i < rows; ++i) {
        for (int j = 1; j < cols; ++j) {
            int match = (sequence1[i - 1] == sequence2[j - 1]) ? MATCH_SCORE : MISMATCH_SCORE;

            int diagonalScore = matrix[i - 1][j - 1] + match;
            int topScore = matrix[i - 1][j] + GAP_SCORE;
            int leftScore = matrix[i][j - 1] + GAP_SCORE;

            matrix[i][j] = std::max({0, diagonalScore, topScore, leftScore});
        }
    }
}


// Function to perform sequence alignment using the Smith-Waterman algorithm with parallelization
void parallelSmithWaterman(const Sequence& sequence1, const Sequence& sequence2, ScoreMatrix& matrix) {
    const int rows = sequence1.size() + 1;
    const int cols = sequence2.size() + 1;

    initializeMatrix(matrix, rows, cols);

    #pragma omp parallel for
    for (int d = 1; d < rows + cols - 1; ++d) {
        int start = std::max(1, d - cols + 1);
        int end = std::min(d, rows - 1);

        for (int i = start; i <= end; ++i) {
            int j = d - i;

            int match = (sequence1[i - 1] == sequence2[j - 1]) ? MATCH_SCORE : MISMATCH_SCORE;

            int diagonalScore = matrix[i - 1][j - 1] + match;
            int topScore = matrix[i - 1][j] + GAP_SCORE;
            int leftScore = matrix[i][j - 1] + GAP_SCORE;

            matrix[i][j] = std::max({0, diagonalScore, topScore, leftScore});
        }
    }
}


int main() {
    Sequence sequence1 = {'A', 'G', 'T', 'A', 'C', 'G', 'A', 'C', 'T'};
    Sequence sequence2 = {'T', 'A', 'C', 'G', 'A', 'G', 'T'};

    ScoreMatrix matrix;
    parallelSmithWaterman(sequence1, sequence2, matrix);

    // Print the score matrix
    for (const auto& row : matrix) {
        for (int score : row) {
            std::cout << score << '\t';
        }
        std::cout << '\n';
    }

    // Find the maximum score and its position
    int maxScore, maxI, maxJ;
    findMaxScore(matrix, maxScore, maxI, maxJ);

    std::cout << "Maximum score: " << maxScore << '\n';
    std::cout << "Position: (" << maxI << ", " << maxJ << ")\n";

    return 0;
}
