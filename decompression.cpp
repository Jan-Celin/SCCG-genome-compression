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
#include <string_view>

using namespace std;
namespace fs = std::filesystem;

void decompress_genome(string compressed_file_path, 
                       string reference_genome_path, 
                       string output_folder, 
                       string& reference_genome, 
                       string& encoded_genome,
                       string& n_indices,
                       string& lowercase_indices) {
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
    
    string line;
    while (getline(reference_file, line)) {
        // Skip header lines starting with '>'
        if (line.empty() || line[0] == '>') {
            continue;
        }
        reference_genome += line;
    }
    reference_file.close();
    reference_genome.erase(remove_if(reference_genome.begin(), reference_genome.end(), ::isspace), reference_genome.end());

    // Read the target genome from the decompressed file
    ifstream decompressed_file(decompressed_file_path);
    if (!decompressed_file.is_open()) {
        cerr << "Error opening decompressed file: " << decompressed_file_path << "\n";
        exit(1);
    }
    
    line = "";
    while (getline(decompressed_file, line)) {
        // Skip header lines starting with '>'
        if (line.empty() || line[0] == '>') {
            continue;
        } else if (line.rfind("N: ", 0) == 0) {
            n_indices = line.substr(3);
        } else if (line.rfind("lower: ", 0) == 0) {
            lowercase_indices = line.substr(7);
        } else {
            encoded_genome += line;
        }
    }
    decompressed_file.close();

    // (Optional) Further processing on the reference and target genomes can be done here.
    cout << "Reference genome loaded (" << reference_genome.size() << " characters)." << "\n";
    cout << "Target genome loaded (" << encoded_genome.size() << " characters)." << "\n";
}

// Reconstruct the target genome using the reference genome and the compressed file.
string reconstruct_genome(const string& reference_genome, 
                          const string& encoded_genome,
                          const string& n_indices_str,
                          const string& lowercase_indices_str) {

    // Process the n_indices as before.
    vector<int> n_indices;
    stringstream n_stream(n_indices_str);
    string token;
    while (getline(n_stream, token, ',')) {
        if (!token.empty()) {
            n_indices.push_back(stoi(token));
        }
    }
    sort(n_indices.begin(), n_indices.end());

    // Updated processing for lowercase_indices:
    vector<int> lowercase_positions;
    int prev_lower = 0;
    size_t pos = 0;
    while (pos < lowercase_indices_str.size()) {
        if (lowercase_indices_str[pos] == '(') {
            // Token is in the form "(delta,length)"
            size_t close = lowercase_indices_str.find(')', pos);
            string tuple_token = lowercase_indices_str.substr(pos + 1, close - pos - 1);
            size_t comma = tuple_token.find(',');
            int delta = stoi(tuple_token.substr(0, comma));
            int len = stoi(tuple_token.substr(comma + 1));
            int lower_start = prev_lower + delta;
            // Add run of lowercase positions.
            for (int j = 0; j < len; j++) {
                lowercase_positions.push_back(lower_start + j);
            }
            prev_lower = lower_start;
            pos = close + 1;
            // Skip a following comma if present.
            if (pos < lowercase_indices_str.size() && lowercase_indices_str[pos] == ',')
                pos++;
        } else {
            // Token is a single delta value.
            size_t comma = lowercase_indices_str.find(',', pos);
            string num_token;
            if (comma == string::npos) {
                num_token = lowercase_indices_str.substr(pos);
                pos = lowercase_indices_str.size();
            } else {
                num_token = lowercase_indices_str.substr(pos, comma - pos);
                pos = comma + 1;
            }
            if (!num_token.empty()) {
                int delta = stoi(num_token);
                int lower_start = prev_lower + delta;
                lowercase_positions.push_back(lower_start);
                prev_lower = lower_start;
            }
        }
    }
    sort(lowercase_positions.begin(), lowercase_positions.end());

    // Updated processing for N indices
    vector<int> n_positions;
    int prev_n = 0;
    size_t pos_n = 0;
    while (pos_n < n_indices_str.size()) {
        if (n_indices_str[pos_n] == '(') {
            // Token is in the form "(delta,length)"
            size_t close = n_indices_str.find(')', pos_n);
            string tuple_token = n_indices_str.substr(pos_n + 1, close - pos_n - 1);
            size_t comma = tuple_token.find(',');
            int delta = stoi(tuple_token.substr(0, comma));
            int len = stoi(tuple_token.substr(comma + 1));
            int n_start = prev_n + delta;
            // Add run of N positions.
            for (int j = 0; j < len; j++) {
                n_positions.push_back(n_start + j);
            }
            prev_n = n_start;
            pos_n = close + 1;
            // Skip a following comma, if present.
            if (pos_n < n_indices_str.size() && n_indices_str[pos_n] == ',')
                pos_n++;
        } else {
            // Token is a single delta value.
            size_t comma = n_indices_str.find(',', pos_n);
            string num_token;
            if (comma == string::npos) {
                num_token = n_indices_str.substr(pos_n);
                pos_n = n_indices_str.size();
            } else {
                num_token = n_indices_str.substr(pos_n, comma - pos_n);
                pos_n = comma + 1;
            }
            if (!num_token.empty()) {
                int delta = stoi(num_token);
                int n_start = prev_n + delta;
                n_positions.push_back(n_start);
                prev_n = n_start;
            }
        }
    }
    sort(n_positions.begin(), n_positions.end());

    string reconstructed_genome;
    size_t i = 0;
    int prev_abs_start = 0; // Holds the previous absolute start position.
    while (i < encoded_genome.size()) {
        if (encoded_genome[i] == '(') {
            size_t end_pos = encoded_genome.find(')', i);
            size_t comma_pos = encoded_genome.find(',', i);
            // Read the delta value.
            int delta = stoi(encoded_genome.substr(i + 1, comma_pos - i - 1));
            // Read the length.
            int length = stoi(encoded_genome.substr(comma_pos + 1, end_pos - comma_pos - 1));
            // Compute the absolute start index.
            int absolute_start = prev_abs_start + delta;
            prev_abs_start = absolute_start;
            // Append the substring from the reference genome.
            reconstructed_genome += reference_genome.substr(absolute_start, length);
            i = end_pos + 1;
        } else {
            reconstructed_genome.push_back(encoded_genome[i]);
            ++i;
        }
    }

    string result;
    size_t temp_index = 0;
    size_t n_pos = 0;

    for (int i = 0; i < reconstructed_genome.size() + n_indices.size(); ++i) {
        if (n_pos < n_indices.size() && n_indices[n_pos] == i) {
            result += 'N';
            ++n_pos;
        } else {
            result += reconstructed_genome[temp_index++];
        }
    }
    for (int l_index : lowercase_positions) {
        result[l_index] = tolower(result[l_index]);
    }
    
    return result;
}

int main(int argc, char* argv[]) {
    cout << "Received " << argc - 1 << " arguments.\n";
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <compressed_file> <reference_file> <output_folder>\n";
        return 1;
    }

    string reference_genome;
    string encoded_genome;

    string n_indices;
    string lowercase_indices;

    try {
        decompress_genome(argv[1], argv[2], argv[3], reference_genome, encoded_genome, n_indices, lowercase_indices);
    } catch (const std::exception& ex) {
        cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    // Reconstruct the target genome using the compressed file.
    string reconstructed_genome;
    try {
        reconstructed_genome = reconstruct_genome(reference_genome, encoded_genome, n_indices, lowercase_indices);
    } catch (const exception &ex) {
        cerr << "Error during reconstruction: " << ex.what() << "\n";
        return 1;
    }
    
    // Write the reconstructed genome to the output file.
    string output_file_path = string(argv[3]) + "/reconstructed_genome.txt";
    ofstream out_file(output_file_path);
    if (!out_file.is_open()) {
        cerr << "Error opening output file: " << argv[3] << "\n";
        return 1;
    }
    out_file << reconstructed_genome;
    out_file.close();
    cout << "Reconstructed genome written (" << reconstructed_genome.size() << " characters) to: " << argv[3] << "\n";
    
    return 0;
}