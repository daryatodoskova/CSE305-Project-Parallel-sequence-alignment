#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <chrono>
#include <string>
#include <cmath>

// Function to initialize the matrix
std::vector<std::vector<int>> initializeMatrix(int rows, int cols) {
    std::vector<std::vector<int>> matrix(rows, std::vector<int>(cols, 0));
    return matrix;
}

// Function to calculate the maximum of three values
int maximum(int a, int b, int c) {
    return std::max(a, std::max(b, c));
}

// Function to perform the Smith-Waterman algorithm
void smithWaterman(const std::string& sequence1, const std::string& sequence2,
                   int matchScore, int mismatchScore, int gapPenalty,
                   std::vector<std::vector<int>>& score,
                   int startDiagonal, int endDiagonal, std::mutex& mtx) {
    // Get the lengths of the sequences
    int m = sequence1.length();
    int n = sequence2.length();

    // Iterate over each diagonal in the specified range
    for (int diagonal = startDiagonal; diagonal <= endDiagonal; ++diagonal) {
        for (int i = 1; i <= m; ++i) {
            int j = diagonal - i;

            if (j < 1 || j > n) {
                continue;
            }

            // Calculate the score for the current cell
            int scoreDiagonal = score[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? matchScore : mismatchScore);
            int scoreUp = score[i - 1][j] - gapPenalty;
            int scoreLeft = score[i][j - 1] - gapPenalty;

            // Take the maximum score
            int maxScore = maximum(scoreDiagonal, scoreUp, scoreLeft);

            // Update the score matrix with thread synchronization
            std::lock_guard<std::mutex> lock(mtx);
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

// Function to perform traceback and find the positions of the optimal local alignment
// std::pair<std::pair<int, int>, std::pair<int, int>> traceback(const std::string& sequence1, const std::string& sequence2,
//                                                               const std::vector<std::vector<int>>& score,
//                                                               int matchScore, int mismatchScore, int gapPenalty) {
//     // Find the cell with the maximum score
//     std::pair<int, int> maxScoreCell = findMaxScoreCell(score);
//     int i = maxScoreCell.first;
//     int j = maxScoreCell.second;

//     // Traceback until reaching a cell with a score of 0
//     while (score[i][j] != 0) {
//         if (score[i][j] == score[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? matchScore : mismatchScore)) {
//             --i;
//             --j;
//         } else if (score[i][j] == score[i - 1][j] - gapPenalty) {
//             --i;
//         } else {
//             --j;
//         }
//     }

//     // Return the positions of the optimal local alignment
//     return std::make_pair(std::make_pair(i, j), std::make_pair(maxScoreCell.first - 1, maxScoreCell.second - 1));
// }

// Function to read sequences from a .txt file
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

std::string extractDNASequence(const std::string& filename) {
    std::ifstream file(filename);
    std::string sequence;

    if (file.is_open()) {
        std::string line;
        bool isSequence = false;

        while (std::getline(file, line)) {
            if (line.empty())
                continue;

            if (line[0] == '>') {
                // Skip header line
                continue;
            }

            if (!isSequence) {
                // Start of sequence
                isSequence = true;
                sequence = line;
            } else {
                // Append subsequent lines to the sequence
                sequence += line;
            }
        }

        // Remove non-ACGT characters
        sequence.erase(std::remove_if(sequence.begin(), sequence.end(), [](char c) {
            return (c != 'A' && c != 'C' && c != 'G' && c != 'T');
        }), sequence.end());

        file.close();
    } else {
        std::cerr << "Failed to open the file: " << filename << std::endl;
    }

    return sequence;
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

int main(int argc, char* argv[]) {
    // Check if two file arguments are provided
    if (argc != 4) {
        std::cout << "Two file arguments are required." << std::endl;
        return 1;
    }
    std::string sequence1, sequence2;

    // load seqs from file

    char *nameSeq1 = argv[1];
    char *nameSeq2 = argv[2];
    int num = std::stoi(argv[3]);
    
    //std::string seq1,seq2;
    std::ifstream StrSeq1;
    StrSeq1.open(nameSeq1);
    sequence1 = stream_to_seq(StrSeq1);

    std::ifstream StrSeq2;
    StrSeq2.open(nameSeq2);
    sequence2 = stream_to_seq(StrSeq2);

    // Retrieve the file names from the command-line arguments
    // const char* file1 = argv[1];
    // const char* file2 = argv[2];
    // Start the timer
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    // Read the sequences from a file
    // std::string filename = "sequences.fasta";
    // std::vector<std::string> sequences = readSequencesFromFile(filename);

    // if (sequences.size() < 2) {
    //     std::cout << "Error: Insufficient number of sequences in the file.\n";
    //     return 1;
    // }

    // std::string sequence1 = extractDNASequence(file1);
    // std::string sequence2 = extractDNASequence(file2);

    // std::cout << "Sequence 1: " << sequence1 << std::endl;

    // std::string sequence1 = "AGTACGCA";
    // std::string sequence2 = "TATGC";
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
    unsigned int numThreads = num;

    // Divide the diagonals evenly among the threads
    unsigned int diagonalsPerThread = (m + n - 1) / numThreads;

    // Create a vector to store the threads
    std::vector<std::thread> threads;

    // Create a mutex for thread synchronization
    std::mutex mtx;

    // Perform the Smith-Waterman algorithm in parallel
    for (unsigned int t = 0; t < numThreads; ++t) {
        int startDiagonal = t * diagonalsPerThread + 1;
        int endDiagonal = (t == numThreads - 1) ? (m + n - 1) : (t + 1) * diagonalsPerThread;
        threads.emplace_back(smithWaterman, std::ref(sequence1), std::ref(sequence2),
                             matchScore, mismatchScore, gapPenalty,
                             std::ref(score), startDiagonal, endDiagonal, std::ref(mtx));
    }

    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
    }

    // End the timer
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    // Calculate the elapsed time
    std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    double runtime = duration.count();

    // Print the runtime
    std::cout << "Runtime: " << runtime << " seconds" << std::endl;

    // std::pair<std::pair<int, int>, std::pair<int, int>> alignmentPositions = traceback(sequence1, sequence2, score, matchScore, mismatchScore, gapPenalty);

    // // Print the positions of the optimal local alignment
    // std::cout << "Start position of Sequence 1: " << alignmentPositions.first.first << std::endl;
    // std::cout << "End position of Sequence 1: " << alignmentPositions.second.first << std::endl;
    // std::cout << "Start position of Sequence 2: " << alignmentPositions.first.second << std::endl;
    // std::cout << "End position of Sequence 2: " << alignmentPositions.second.second << std::endl;

    // Find the optimal local alignment
    std::pair<std::string, std::string> alignment = traceback(sequence1, sequence2, score, matchScore, mismatchScore, gapPenalty);

    // Print the alignment
    std::cout << "Sequence 1: " << alignment.first << std::endl;
    std::cout << "Sequence 2: " << alignment.second << std::endl;

    return 0;
}

int max(int a, int b, int c) {
    return std::max(a, std::max(b, c));
}


#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <chrono>
#include <string>
#include <cmath>

using namespace std;


int calculateScore(char a, char b) {
    if (a == b) {
        return 2;
    } else {
        return -1;
    }
}

void AuxScoreMatrix(const std::string& seq1, const std::string& seq2, int k, int startChunk, int endChunk, std::vector<std::vector<int>>& scoreMatrix) {
    int m = seq2.length();

    for (int i = startChunk; i < endChunk; ++i) {
        for (int j = 1; j <= m; ++j) {
            int matchScore = scoreMatrix[i - 1][j - 1] + calculateScore(seq1[i - 1], seq2[j - 1]);
            int deleteScore = scoreMatrix[i - 1][j] + 1;
            int insertScore = scoreMatrix[i][j - 1] + 1;
            int zeroScore = 0;

            scoreMatrix[i][j] = max(matchScore, deleteScore, insertScore, zeroScore);
        }
    }
}

void smithWaterman(const std::string& seq1, const std::string& seq2,
                   int matchScore, int mismatchScore, int gapPenalty,
                   std::vector<std::vector<int>>& score,
                   int startDiagonal, int endDiagonal, std::mutex& mtx, int numThreads) {
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
            workers[i] = std::thread(&AuxScoreMatrix, seq1, seq2, k, startChunk, endChunk, std::ref(scoreMatrix));
            startChunk = endChunk;
            endChunk = min(startChunk + chunkSize, endRow + 1);
        }

        AuxScoreMatrix(seq1, seq2, k, startChunk, endChunk, scoreMatrix);

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