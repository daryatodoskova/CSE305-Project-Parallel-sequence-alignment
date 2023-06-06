#include <iostream>
#include <fstream>
#include <cstdlib>

void runCommandsInTerminal(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Failed to open commands file: " << filename << std::endl;
        return;
    }

    std::string command;
    while (std::getline(file, command)) {
        int exitCode = system(command.c_str());
        if (exitCode != 0) {
            std::cerr << "Failed to execute command: " << command << std::endl;
            break;
        }
    }

    file.close();
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Please provide a filename as an argument." << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    runCommandsInTerminal(filename);

    return 0;
}






