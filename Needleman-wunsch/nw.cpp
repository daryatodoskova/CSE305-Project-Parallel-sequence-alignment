#include <iostream>
#include <fstream>
#include <string>
#include <vector>
//#include <omp.h>
// #include "/usr/local/opt/libomp/include/omp.h"
#include <thread>
#include <mutex>
#include <cmath>

using namespace std;

//SCORING SCHEME
#define MATCH_SCORE 2       //for aligning identical characters
#define MISMATCH_SCORE -1   //for aligning diff character
#define GAP_SCORE -2        //for introducing gaps in the alignment

int max(int a, int b, int c) {
    return std::max(a, std::max(b, c));
}

int calculateScore(char a, char b) {
    if (a == b) {
        return MATCH_SCORE;
    } else {
        return MISMATCH_SCORE;
    }
}

void AuxScoreMatrix(const std::string seq1, const std::string seq2, int k, int startChunk, int endChunk, std::vector<std::vector<int> >& scoreMatrix){

    int m = seq2.length();

    for (int i = startChunk; i < endChunk; ++i) {

        for (int j = k - i + 1; j >= 1 && j <= m ; ++j) {
            int matchScore = scoreMatrix[i - 1][j - 1] + calculateScore(seq1[i - 1], seq2[j - 1]);
            int deleteScore = scoreMatrix[i - 1][j] + GAP_SCORE;
            int insertScore = scoreMatrix[i][j - 1] + GAP_SCORE;

            scoreMatrix[i][j] = max(matchScore, deleteScore, insertScore);
        }
    }
}

// Function to perform parallel sequence alignment using NW
void needlemanWunsch(const std::string& seq1, const std::string& seq2, std::vector<std::vector<int> >& scoreMatrix, int numThreads) {
    // DEFINED IN MAIN()
    int n = seq1.length();
    int m = seq2.length();

    // DEFINING THE NUMBER OF THREADS
    int chunkSize = ceil(max(n,m) / numThreads);
    std::cout << "chunkSize " <<  chunkSize << " - numThreads " << numThreads << std::endl;

    scoreMatrix.resize(n + 1, std::vector<int>(m + 1));

    int init_row = n + 1;
    int init_col = m + 1;

    // Create and initialize the score matrix based on the gap penalty
    for (int i = 0; i <= n; ++i)
        scoreMatrix[i][0] = i * GAP_SCORE;
    for (int j = 0; j <= m; ++j)
        scoreMatrix[0][j] = j * GAP_SCORE;

    // PARALLELIZED
    for (int k = 1; k <= n + m - 1; ++k) {
    
        int startRow = max(1, k - m + 1);
        int endRow = min(n, k);
        int startChunk = startRow;
        int endChunk = min(startChunk + chunkSize, endRow + 1);

        std::vector<std::thread> workers(numThreads - 1);

        for (size_t i = 0 ; i < numThreads - 1; ++i){
            workers[i] = std::thread(&AuxScoreMatrix, seq1, seq2, k, startChunk, endChunk, std::ref(scoreMatrix));
            startChunk = endChunk;
            endChunk = min(startChunk + chunkSize, endRow + 1);
        }

        AuxScoreMatrix(seq1, seq2, k, startChunk, endChunk, scoreMatrix);
        
        for (size_t i = 0; i < numThreads - 1; ++i) {
            workers[i].join();
        }   

    }


    // Perform traceback to determine the aligned sequences
    string alignedseq1 = "";
    string alignedseq2 = "";
    int i = n;
    int j = m;

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + ((seq1[i - 1] == seq2[j - 1]) ? MATCH_SCORE : MISMATCH_SCORE)) {
            alignedseq1 = seq1[i - 1] + alignedseq1;
            alignedseq2 = seq2[j - 1] + alignedseq2;
            --i;
            --j;
        } else if (i > 0 && scoreMatrix[i][j] == scoreMatrix[i - 1][j] + GAP_SCORE) {
            alignedseq1 = seq1[i - 1] + alignedseq1;
            alignedseq2 = '-' + alignedseq2;
            --i;
        } else {
            alignedseq1 = '-' + alignedseq1;
            alignedseq2 = seq2[j - 1] + alignedseq2;
            --j;
        }
    }

    cout << "Aligned Sequence 1: " << alignedseq1 << endl;
    cout << "Aligned Sequence 2: " << alignedseq2 << endl;
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

// void saveScoreMatrix(const vector<vector<int>>& score, const string& filename) {
//     ofstream outputFile(filename);
//     if (!outputFile) {
//         cout << "Error creating output file!" << endl;
//         return;
//     }

//     int m = score.size();
//     int n = score[0].size();

//     for (int i = 0; i < m; ++i) {
//         for (int j = 0; j < n; ++j) {
//             outputFile << score[i][j] << "\t";
//         }
//         outputFile << endl;
//     }

//     outputFile.close();
//     cout << "Alignment score matrix saved to " << filename << endl;
// }

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    std::string inputFilePath = argv[1];
    std::ifstream inputFile(inputFilePath);
    if (!inputFile) {
        std::cerr << "Failed to open input file: " << inputFilePath << std::endl;
        return 1;
    }

    std::string seq1, seq2;
    if (!(inputFile >> seq1 >> seq2)) {
        std::cerr << "Failed to read sequences from the input file." << std::endl;
        return 1;
    }

    // load seqs from file

    char *nameSeq1 = argv[1];
    char *nameSeq2 = argv[2];
    
    //std::string seq1,seq2;
    std::ifstream StrSeq1;
    StrSeq1.open(nameSeq1);
    seq1 = stream_to_seq(StrSeq1);

    std::ifstream StrSeq2;
    StrSeq2.open(nameSeq2);
    seq2 = stream_to_seq(StrSeq2);

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<int> > scoreMatrix;

    int n = seq1.length();
    int m = seq2.length();

    // DEFINING THE NUMBER OF THREADS
    int numThreads = 8;

    needlemanWunsch(seq1, seq2, scoreMatrix, numThreads);

    auto end = std::chrono::high_resolution_clock::now();

    const string& filename = "output.txt";

    ofstream outputFile(filename);

    // Print alignment scoreMatrix
    for (const auto& row : scoreMatrix) {
        for (int score : row) {
            // std::cout << score << "\t";
            outputFile << score << "\t";
        }
        outputFile << endl;
        // std::cout << std::endl;
    }

    outputFile.close();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  
    std::cout << "Function duration: " << duration.count() << " ms" << std::endl;

    return 0;
}

