#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>

using namespace std;
namespace fs = std::filesystem;

void decompress_genome(string compressed_file_path, string reference_genome_path, string output_folder) {
    // Ensure the output folder exists
    if (!fs::exists(output_folder)) {
        fs::create_directory(output_folder);
    }

    // Construct the 7z extraction command and run it
    string command = "7z e \"" + compressed_file_path + "\" -o\"" + output_folder + "\" -y";
    int result = system(command.c_str());
    if (result != 0) {
        cerr << "Error decompressing file: " << compressed_file_path << "\n";
        exit(1);
    } else {
        cout << "File successfully decompressed to: " << output_folder << "\n";
    }

    // Determine the path of the decompressed file (using the stem of compressed_file_path)
    fs::path compPath(compressed_file_path);
    string decompressed_file_path = output_folder + "/" + compPath.stem().string();
    cout << "Found decompressed file: " << decompressed_file_path << "\n";

    // Read the reference genome from the reference file
    ifstream reference_file(reference_genome_path);
    if (!reference_file.is_open()) {
        cerr << "Error opening reference file: " << reference_genome_path << "\n";
        exit(1);
    }
    stringstream refBuffer;
    refBuffer << reference_file.rdbuf();
    string reference_genome = refBuffer.str();
    reference_file.close();

    // Read the target genome from the decompressed file
    ifstream decompressed_file(decompressed_file_path);
    if (!decompressed_file.is_open()) {
        cerr << "Error opening decompressed file: " << decompressed_file_path << "\n";
        exit(1);
    }
    stringstream targetBuffer;
    targetBuffer << decompressed_file.rdbuf();
    string target_genome = targetBuffer.str();
    decompressed_file.close();

    // (Optional) Further processing on the reference and target genomes can be done here.
    cout << "Reference genome loaded (" << reference_genome.size() << " characters)." << "\n";
    cout << "Target genome loaded (" << target_genome.size() << " characters)." << "\n";
}

// Reconstruct the target genome using the reference genome and the compressed file.
string reconstruct_genome(const string& reference_genome, const string& compressed_file_path) {
    ifstream comp_file(compressed_file_path);
    if (!comp_file.is_open()) {
        throw runtime_error("Error opening compressed file: " + compressed_file_path);
    }

    string reconstructed_genome;
    string token;
    while (comp_file >> token) {
        // Check if token is a match token (starts with '(')
        if (!token.empty() && token.front() == '(') {
            // Remove surrounding parentheses. Assumes token ends with ')'
            if (token.back() != ')') {
                throw runtime_error("Invalid token format (missing ')'): " + token);
            }
            string content = token.substr(1, token.size() - 2);
            size_t comma_pos = content.find(',');
            if (comma_pos == string::npos) {
                throw runtime_error("Invalid token format (missing comma): " + token);
            }
            string start_ref_str = content.substr(0, comma_pos);
            string length_str = content.substr(comma_pos + 1);
            int start_ref = stoi(start_ref_str);
            int seg_length = stoi(length_str);
            // Extract matching segment from the reference genome.
            if(start_ref < 0 || start_ref + seg_length > (int)reference_genome.size()){
                throw runtime_error("Token range is out of bounds: " + token);
            }
            string substring = reference_genome.substr(start_ref, seg_length);
            reconstructed_genome.append(substring);
        } else {
            // Token is a literal mismatch; append as-is.
            reconstructed_genome.append(token);
        }
    }
    comp_file.close();
    return reconstructed_genome;
}

int main(int argc, char* argv[]) {
    cout << "Received " << argc - 1 << " arguments.\n";
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <compressed_file> <reference_file> <output_folder>\n";
        return 1;
    }
    try {
        decompress_genome(argv[1], argv[2], argv[3]);
    } catch (const std::exception& ex) {
        cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    // Read entire reference genome from reference file.
    ifstream ref_file(argv[2]);
    if (!ref_file.is_open()) {
        cerr << "Error opening reference file: " << argv[2] << "\n";
        return 1;
    }
    stringstream refBuffer;
    refBuffer << ref_file.rdbuf();
    string reference_genome = refBuffer.str();
    ref_file.close();

    // Reconstruct the target genome using the compressed file.
    string reconstructed_genome;
    try {
        reconstructed_genome = reconstruct_genome(reference_genome, argv[1]);
    } catch (const exception &ex) {
        cerr << "Error during reconstruction: " << ex.what() << "\n";
        return 1;
    }
    
    // Write the reconstructed genome to the output file.
    ofstream out_file(argv[3]);
    if (!out_file.is_open()) {
        cerr << "Error opening output file: " << argv[3] << "\n";
        return 1;
    }
    out_file << reconstructed_genome;
    out_file.close();
    cout << "Reconstructed genome written (" << reconstructed_genome.size() << " characters) to: " << argv[3] << "\n";
    
    return 0;
}