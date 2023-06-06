#include <iostream>
#include <vector>
#include <thread>
#include <mutex>

// Function to calculate the maximum of three values
int maximum(int a, int b, int c) {
    return std::max(std::max(a, b), c);
}

// Function to perform the Smith-Waterman algorithm (parallel version)
void smithWatermanParallel(const std::string& sequence1, const std::string& sequence2) {
    const int m = sequence1.length();
    const int n = sequence2.length();

    std::vector<std::vector<int>> score(m + 1, std::vector<int>(n + 1, 0));
    std::mutex mtx; // Mutex for shared matrix

    // Define the number of threads
    const int numThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads(numThreads);

    // Compute the score matrix (parallel computation)
    for (int k = 0; k < numThreads; ++k) {
        threads[k] = std::thread([&, k]() {
            for (int i = 1 + k; i <= m; i += numThreads) {
                for (int j = 1; j <= n; ++j) {
                    int match = score[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? 1 : -1);
                    int deleteGap = score[i - 1][j] - 1;
                    int insertGap = score[i][j - 1] - 1;
                    score[i][j] = maximum(match, deleteGap, insertGap);
                }
            }
        });
    }

    // Join the threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Find the maximum score
    int maxScore = 0;
    int maxI = 0;
    int maxJ = 0;
    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (score[i][j] > maxScore) {
                maxScore = score[i][j];
                maxI = i;
                maxJ = j;
            }
        }
    }

    // Backtrace to find the aligned sequences
    std::string alignedSeq1, alignedSeq2;
    while (maxI > 0 && maxJ > 0 && score[maxI][maxJ] != 0) {
        if (score[maxI][maxJ] == score[maxI - 1][maxJ - 1] + (sequence1[maxI - 1] == sequence2[maxJ - 1] ? 1 : -1)) {
            alignedSeq1 = sequence1[maxI - 1] + alignedSeq1;
            alignedSeq2 = sequence2[maxJ - 1] + alignedSeq2;
            maxI--;
            maxJ--;
        } else if (score[maxI][maxJ] == score[maxI - 1][maxJ] - 1) {
            alignedSeq1 = sequence1[maxI - 1] + alignedSeq1;
            alignedSeq2 = '-' + alignedSeq2;
            maxI--;
        } else {
            alignedSeq1 = '-' + alignedSeq1;
            alignedSeq2 = sequence2[maxJ - 1] + alignedSeq2;
            maxJ--;
        }
    }

    // Print the aligned sequences and the maximum score
    std::cout << "Aligned Sequence 1: " << alignedSeq1 << std::endl;
    std::cout << "Aligned Sequence 2: " << alignedSeq2 << std::endl;
    std::cout << "Maximum Score: " << maxScore << std::endl;
}

int main() {
    std::string sequence1 = "AGTACGCA";
    std::string sequence2 = "TATGC";

    smithWatermanParallel(sequence1, sequence2);

    return 0;
}
