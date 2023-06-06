#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <chrono>
#include <string>
#include <cmath>
#include <sstream>

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
                   std::vector<std::vector<int>>& score, int numThreads) 
    {
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

void writeOutput(const std::string& seq1, const std::string& seq2, int numThreads, double runtime, const std::string& alignment1, const std::string& alignment2) {
    std::ofstream outputFile("sw_comp.txt", std::ios::app);
    if (!outputFile) {
        std::cout << "Failed to open the output file." << std::endl;
        return;
    }

    // outputFile << "Sequences:\n";
    // outputFile << "Sequence 1: " << seq1 << "\n";
    // outputFile << "Sequence 2: " << seq2 << "\n";

    outputFile << "Number of Threads: " << numThreads << "\n";

    outputFile << "Alignment:\n";
    outputFile << "Aligned Sequence 1: " << alignment1 << "\n";
    outputFile << "Aligned Sequence 2: " << alignment2 << "\n";

    outputFile << "Runtime: " << runtime << " seconds\n";

    outputFile << "-------------------------------------------------------------------\n";

    outputFile.close();
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cout << "Two file arguments and one integer are required." << std::endl;
        return 1;
    }

    std::string seq1, seq2;

    // load seqs from file
    char *nameSeq1 = argv[1];
    char *nameSeq2 = argv[2];
    int num = std::stoi(argv[3]);


    std::ifstream StrSeq1;
    StrSeq1.open(nameSeq1);
    seq1 = stream_to_seq(StrSeq1);

    std::ifstream StrSeq2;
    StrSeq2.open(nameSeq2);
    seq2 = stream_to_seq(StrSeq2);

    // Start the timer
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    // Specify the scoring scheme and other parameters
    int matchScore = 2;
    int mismatchScore = -1;
    int gapPenalty = -2;

    // Get the lengths of the sequences
    int n = seq1.length();
    int m = seq2.length();

    int chunkSize = ceil(max(n,m) / num);

    // Initialize the score matrix
    std::vector<std::vector<int>> score = initializeMatrix(m + 1, n + 1);

    // DEFINING THE NUMBER OF THREADS
    int numThreads = num;

    smithWaterman(seq1, seq2, matchScore, mismatchScore, gapPenalty, score, num);

    // End the timer
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    // Calculate the elapsed time
    std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    double runtime = duration.count();

    // Print the runtime
    std::cout << "Runtime: " << runtime << " seconds" << std::endl;

    // Perform traceback to determine the aligned sequences
    std::pair<std::string, std::string> alignment = traceback(seq1, seq2, score, matchScore, mismatchScore, gapPenalty);
    std::string alignment1 = alignment.first;
    std::string alignment2 = alignment.second;

    // Write the output to a file
    writeOutput(seq1, seq2, numThreads, runtime, alignment1, alignment2);

    // const string& filename = "output.txt";

    // ofstream outputFile(filename);

    // // Print alignment scoreMatrix
    // for (const auto& row : score) {
    //     for (int scoreM : row) {
    //         //std::cout << score << "\t";
    //         outputFile << scoreM << "\t";
    //     }
    //     outputFile << endl;
    //     //std::cout << std::endl;
    // }

    // outputFile.close();

    std::ofstream csvFile("sw.csv", std::ios::app);
    csvFile << n << "," << m << "," << chunkSize << "," << runtime <<"," << num << "\n";
    csvFile.close();

    return 0;
}
