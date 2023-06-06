#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <chrono>
#include <string>
#include <cmath>

using namespace std;

// Function to initialize the matrix
std::vector<std::vector<int>> initializeMatrix(int rows, int cols) {
    std::vector<std::vector<int>> matrix(rows, std::vector<int>(cols, 0));
    return matrix;
}

int max(int a, int b, int c) {
    return std::max(a, std::max(b, c));
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

int calculateScore(char a, char b) {
    if (a == b) {
        return 2;
    } else {
        return 0;
    }
}

void AuxScoreMatrix(const std::string seq1, const std::string seq2, int k, int startChunk, int endChunk, std::vector<std::vector<int> >& scoreMatrix){

    int m = seq2.length();

    for (int i = startChunk; i < endChunk; ++i) {

        for (int j = k - i + 1; j >= 1 && j <= m ; ++j) {
            int matchScore = scoreMatrix[i - 1][j - 1] + calculateScore(seq1[i - 1], seq2[j - 1]);
            int deleteScore = scoreMatrix[i - 1][j] - 2;
            int insertScore = scoreMatrix[i][j - 1] -2 ;

            scoreMatrix[i][j] = max(0, max(matchScore, deleteScore, insertScore));
        }
    }
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

void smithWaterman(const std::string& seq1, const std::string& seq2,
                   int matchScore, int mismatchScore, int gapPenalty,
                   std::vector<std::vector<int>>& score, int numThreads) {
    // Get the lengths of the sequences
    int n = seq1.length();
    int m = seq2.length();

    // DEFINING THE NUMBER OF THREADS
    int chunkSize = ceil(max(n,m) / numThreads);
    std::cout << "chunkSize " <<  chunkSize << " - numThreads " << numThreads << std::endl;

    score.resize(n + 1, std::vector<int>(m + 1));

    int init_row = n + 1;
    int init_col = m + 1;

    // Create and initialize the score matrix to 0
    for (int i = 0; i <= n; ++i)
        score[i][0] = 0;
    for (int j = 0; j <= m; ++j)
        score[0][j] = 0;

    // PARALLELIZED
    for (int k = 1; k <= n; ++k) {
        int startRow = max(1, k - chunkSize + 1);
        int endRow = min(n, k);
        int startChunk = startRow;
        int endChunk = min(startChunk + chunkSize, endRow + 1);

        std::vector<std::thread> workers(numThreads - 1);

        for (size_t i = 0; i < numThreads - 1; ++i) {
            workers[i] = std::thread(&AuxScoreMatrix, seq1, seq2, k, startChunk, endChunk, std::ref(score));
            startChunk = endChunk;
            endChunk = min(startChunk + chunkSize, endRow + 1);
        }

        AuxScoreMatrix(seq1, seq2, k, startChunk, endChunk, score);

        for (size_t i = 0; i < numThreads - 1; ++i) {
            workers[i].join();
        }
    }

    // Find the cell with the maximum score
    std::pair<int, int> maxScoreCell = findMaxScoreCell(score);
    int i = maxScoreCell.first;
    int j = maxScoreCell.second;

    // Perform traceback to determine the aligned sequences
    string alignedSeq1 = "";
    string alignedSeq2 = "";

    // Find the optimal local alignment
    std::pair<std::string, std::string> alignment = traceback(seq1, seq2, score, matchScore, mismatchScore, gapPenalty);

        // Print the alignment
    std::cout << "Sequence 1: " << alignment.first << std::endl;
    std::cout << "Sequence 2: " << alignment.second << std::endl;
    }


std::string stream_to_seq(std::ifstream& f)
{
    std::string seq;
    char line[10000];
    while(f.good())
    {
        f.getline(line,10000);
        if( line[0] == 0 || line[0]=='#' ) 
        {
            continue;
        }
        for(int i = 1; line[i] != 0; i++)
        {
            int ch = toupper(line[i]);
            if(ch != 'A' && ch != 'G' && ch != 'C' && ch != 'T')
            {
                continue;
            }
            seq.push_back(char(ch));
        }
    }
    return seq;
}

void writeAlignmentToFile(const std::string& alignedSeq1, const std::string& alignedSeq2, int comparisonNumber, const std::string& filename) {
    std::ofstream outputFile(filename, std::ios::app);  // Open the file in append mode
    if (!outputFile) {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    outputFile << "Comparison " << comparisonNumber << ":" << std::endl;
    outputFile << "Aligned Sequence 1: " << alignedSeq1 << std::endl;
    outputFile << "Aligned Sequence 2: " << alignedSeq2 << std::endl;
    outputFile << std::endl;  // Add a newline between alignments

    outputFile.close();
    std::cout << "Alignment result for comparison " << comparisonNumber << " appended to " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file1> <input_file2> <num_threads>" << std::endl;
        return 1;
    }

    std::string inputFilePath1 = argv[1];
    std::string inputFilePath2 = argv[2];
    int numThreads = std::stoi(argv[3]);

    std::ifstream inputFile1(inputFilePath1);
    if (!inputFile1) {
        std::cerr << "Failed to open input file: " << inputFilePath1 << std::endl;
        return 1;
    }

    std::ifstream inputFile2(inputFilePath2);
    if (!inputFile2) {
        std::cerr << "Failed to open input file: " << inputFilePath2 << std::endl;
        return 1;
    }

    std::string seq1, seq2;
    if (!(inputFile1 >> seq1) || !(inputFile2 >> seq2)) {
        std::cerr << "Failed to read sequences from the input files." << std::endl;
        return 1;
    }

    inputFile1.close();
    inputFile2.close();

    std::vector<std::vector<int>> scoreMatrix;

    int n = seq1.length();
    int m = seq2.length();

    // Perform Smith-Waterman algorithm
    smithWaterman(seq1, seq2, scoreMatrix, numThreads);

    std::string alignedSeq1, alignedSeq2;
    tracebackSmithWaterman(seq1, seq2, scoreMatrix, alignedSeq1, alignedSeq2);

    std::string outputFilename = "alignment_results.txt";
    writeAlignmentToFile(alignedSeq1, alignedSeq2, outputFilename);

    std::cout << "Alignment result appended to " << outputFilename << std::endl;

    return 0;
}
