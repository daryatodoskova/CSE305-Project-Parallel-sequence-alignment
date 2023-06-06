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

int main() {
    std::string filename = "Sequences/seq.txt";
    std::string sequence = extractDNASequence(filename);

    std::cout << "Sequence: " << sequence << std::endl;

    return 0;
}
