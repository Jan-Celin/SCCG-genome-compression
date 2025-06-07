#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <set>
#include <functional>
#include <climits>
#include <filesystem>
#include <unordered_map>

using namespace std;

void read_genomes_from_files(const string& reference_file, const string& target_file, string& reference_genome, string& target_genome) {
    // Open the reference file.
    ifstream ref_stream(reference_file);
    if (!ref_stream.is_open()) {
        cerr << "Error opening reference file: " << reference_file << "\n";
        exit(1);
    }

    stringstream ref_buffer;
    ref_buffer << ref_stream.rdbuf();
    reference_genome = ref_buffer.str();
    ref_stream.close();
    
    reference_genome.erase(remove_if(reference_genome.begin(), reference_genome.end(), ::isspace), reference_genome.end());

    // Open the target file.
    ifstream target_stream(target_file);
    if (!target_stream.is_open()) {
        cerr << "Error opening target file: " << target_file << "\n";
        exit(1);
    }

    // Read the target genome.
    stringstream target_buffer;
    target_buffer << target_stream.rdbuf();
    target_genome = target_buffer.str();
    target_stream.close();

    target_genome.erase(remove_if(target_genome.begin(), target_genome.end(), ::isspace), target_genome.end());
}

void compress_genome(const string& reference_file, const string& target_file, const string& output_folder) {
    string reference_genome;
    string target_genome;

    // Read genomes from files.
    read_genomes_from_files(reference_file, target_file, reference_genome, target_genome);

    cout << reference_genome << "\n";
    cout << target_genome << "\n";
    
    string output_file = output_folder + "/compressed_genome.txt";
}

// Main function.
int main(int argc, char* argv[]) {
    cout << "Received " << argc << " arguments.\n";
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <reference_file> <target_file> <output_folder>\n";
        return 1;
    }
    try {
        // Ensure output folder exists.
        if (!filesystem::exists(argv[3])) {
            filesystem::create_directory(argv[3]);
        }
        cout << "Successfully created output folder: " << argv[3] << "\n";

        compress_genome(argv[1], argv[2], argv[3]);

    } catch (const std::exception& ex) {
        cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
