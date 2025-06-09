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
#include <chrono>

using namespace std;
namespace fs = std::filesystem;

void decompress_genome(string compressed_file_path, 
                       string reference_genome_path, 
                       string output_folder, 
                       string& reference_genome, 
                       string& encoded_genome,
                       string& n_indices,
                       string& lowercase_indices,
                       string& target_header) {
    if (!fs::exists(output_folder)) {
        fs::create_directory(output_folder);
    }

    // Dekomprimiranje zip arhive.
    string command = "7z e \"" + compressed_file_path + "\" -o\"" + output_folder + "\" -y";
    int result = system(command.c_str());
    if (result != 0) {
        cerr << "Greska pri dekompresiji: " << compressed_file_path << "\n";
        exit(1);
    } else {
        cout << "Datoteka uspjesno dekomprimirana: " << output_folder << "\n";
    }

    fs::path compPath(compressed_file_path);
    string decompressed_file_path = output_folder + "/" + compPath.stem().string();
    cout << "Dekomprimirana datoteka: " << decompressed_file_path << "\n";

    ifstream reference_file(reference_genome_path);
    if (!reference_file.is_open()) {
        cerr << "Greska pri otvaranju reference: " << reference_genome_path << "\n";
        exit(1);
    }
    string line;
    while (getline(reference_file, line)) {
        if (line.empty() || line[0] == '>') continue; // preskoci headere
        reference_genome += line;
    }
    reference_file.close();
    reference_genome.erase(remove_if(reference_genome.begin(), reference_genome.end(), ::isspace), reference_genome.end());

    ifstream decompressed_file(decompressed_file_path);
    if (!decompressed_file.is_open()) {
        cerr << "Greska pri otvaranju datoteke: " << decompressed_file_path << "\n";
        exit(1);
    }

    string first_line, lowercase_indices_line, n_indices_line, encoded_genome_line;

    if (!getline(decompressed_file, first_line)) {
        cerr << "Greska pri citanju datoteke: " << decompressed_file_path << "\n";
        exit(1);
    }

    if (!first_line.empty() && first_line[0] == '>') {
        target_header = first_line;
        if (!getline(decompressed_file, lowercase_indices_line)) {
            cerr << "Greska pri citanju indeksa malih slova: " << decompressed_file_path << "\n";
            exit(1);
        }
        if (!getline(decompressed_file, n_indices_line)) {
            cerr << "Greska pri citanju indeksa nepoznatih nukleotida: " << decompressed_file_path << "\n";
            exit(1);
        }
        if (!getline(decompressed_file, encoded_genome_line)) {
            cerr << "Greska pri citanju kodiranog genoma: " << decompressed_file_path << "\n";
            exit(1);
        }
    } else {
        lowercase_indices_line = first_line;
        if (!getline(decompressed_file, n_indices_line)) {
            cerr << "Greska pri citanju indeksa nepoznatih nukleotida: " << decompressed_file_path << "\n";
            exit(1);
        }
        if (!getline(decompressed_file, encoded_genome_line)) {
            cerr << "Greska pri citanju kodiranog genoma: " << decompressed_file_path << "\n";
            exit(1);
        }
    }

    lowercase_indices = lowercase_indices_line;
    n_indices = n_indices_line;
    encoded_genome = encoded_genome_line;

    decompressed_file.close();

    if(n_indices_line == ",") {
        n_indices_line = "";
    } else {
        reference_genome.erase(remove(reference_genome.begin(), reference_genome.end(), 'N'), reference_genome.end());
    }
    transform(reference_genome.begin(), reference_genome.end(), reference_genome.begin(), ::toupper);

    cout << "Ucitan referentni genom (" << reference_genome.size() << " znakova)." << "\n";
    cout << "Ucitan ciljni genom (" << encoded_genome.size() << " znakova)." << "\n";
}

// Reconstruct the target genome using the reference genome and the compressed file.
string reconstruct_genome(const string& reference_genome, 
                          const string& encoded_genome,
                          const string& n_indices_str,
                          const string& lowercase_indices_str) {

    cout << "DEBUG: Starting genome reconstruction" << "\n";
    cout << "DEBUG: Reference genome size: " << reference_genome.size() << "\n";
    cout << "DEBUG: Encoded genome size: " << encoded_genome.size() << "\n";

    vector<int> lowercase_positions;
    int prev_lower = 0;
    size_t pos = 0;
    while (pos < lowercase_indices_str.size()) {
        if (lowercase_indices_str[pos] == '(') {
            size_t close = lowercase_indices_str.find(')', pos);
            string tuple_token = lowercase_indices_str.substr(pos + 1, close - pos - 1);
            size_t comma = tuple_token.find(',');
            int delta = stoi(tuple_token.substr(0, comma));
            int len = stoi(tuple_token.substr(comma + 1));
            int lower_start = prev_lower + delta;
            
            for (int j = 0; j < len; j++) {
                lowercase_positions.push_back(lower_start + j);
            }
            prev_lower = lower_start;
            pos = close + 1;
            if (pos < lowercase_indices_str.size() && lowercase_indices_str[pos] == ',')
                pos++;
        } else {
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

    vector<int> n_positions;
    int prev_n = 0;
    size_t pos_n = 0;
    while (pos_n < n_indices_str.size()) {
        if (n_indices_str[pos_n] == '(') {
            
            size_t close = n_indices_str.find(')', pos_n);
            string tuple_token = n_indices_str.substr(pos_n + 1, close - pos_n - 1);
            size_t comma = tuple_token.find(',');
            int delta = stoi(tuple_token.substr(0, comma));
            int len = stoi(tuple_token.substr(comma + 1));
            int n_start = prev_n + delta;
            
            for (int j = 0; j < len; j++) {
                n_positions.push_back(n_start + j);
            }
            prev_n = n_start;
            pos_n = close + 1;
            
            if (pos_n < n_indices_str.size() && n_indices_str[pos_n] == ',')
                pos_n++;
        } else {
            
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

    // Konacna rekonstrukcija
    string reconstructed_genome;
    size_t i = 0;
    int prev_abs_start = 0;
    while (i < encoded_genome.size()) {
        if (encoded_genome[i] == '(') {
            size_t end_pos = encoded_genome.find(')', i);
            size_t comma_pos = encoded_genome.find(',', i);
            
            int delta = stoi(encoded_genome.substr(i + 1, comma_pos - i - 1));
            int length = stoi(encoded_genome.substr(comma_pos + 1, end_pos - comma_pos - 1));
            int absolute_start = prev_abs_start + delta;
            
            prev_abs_start = absolute_start;
            if (absolute_start + length > (int)reference_genome.size()) {
                cerr << "ERROR: absolute_start + length (" 
                     << absolute_start + length 
                     << ") exceeds reference genome size (" 
                     << reference_genome.size() << ")" << "\n";
                exit(1);
            }
            reconstructed_genome += reference_genome.substr(absolute_start, length);
            i = end_pos + 1;
        } else {
            reconstructed_genome.push_back(encoded_genome[i]);
            ++i;
        }
    }
    cout << "DEBUG: Reconstructed genome length after token decoding: " 
         << reconstructed_genome.size() << "\n";

    // Umetni nepoznate nukleotide ('N') na odgovarajuce pozicije.
    string result;
    size_t temp_index = 0;
    size_t n_pos = 0;
    for (int i = 0; i < (int)(reconstructed_genome.size() + n_positions.size()); ++i) {
        if (n_pos < n_positions.size() && n_positions[n_pos] == i) {
            result += 'N';
            
            ++n_pos;
        } else {
            result += reconstructed_genome[temp_index++];
        }
    }

    // Pretvori mala slova u velika (gdje treba)
    for (int l_index : lowercase_positions) {
        if (l_index < (int)result.size()) {
            result[l_index] = tolower(result[l_index]);
        } else {
            cout << "WARNING: Lowercase index " << l_index << " out of bounds (" 
                 << result.size() << ")\n";
        }
    }
    
    cout << "DEBUG: Final reconstructed genome length: " << result.size() << "\n";

    string formatted;
    size_t chunk_size = 50;
    for (size_t pos = 0; pos < result.size(); pos += chunk_size) {
        formatted += result.substr(pos, chunk_size);
        if (pos + chunk_size < result.size()) {
            formatted += "\n";
        }
    }
    result = formatted + "\n";

    cout << "DEBUG: Reconstructed genome ready with " << result.size() << " characters.\n";

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
    string target_header;

    cout << "Starting decompression.";
    auto start = chrono::high_resolution_clock::now();
    
    try {
        decompress_genome(argv[1], argv[2], argv[3], reference_genome, encoded_genome, n_indices, lowercase_indices, target_header);
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
    auto end = chrono::high_resolution_clock::now();
    
    // Write the reconstructed genome to the output file.
    string output_file_path = string(argv[3]) + "/reconstructed_genome.fa";
    ofstream out_file(output_file_path);
    if (!out_file.is_open()) {
        cerr << "Error opening output file: " << argv[3] << "\n";
        return 1;
    }
    out_file << target_header << "\n" << reconstructed_genome;
    out_file.close();
    cout << "Reconstructed genome written (" << reconstructed_genome.size() << " characters) to: " << argv[3] << "\n";

    chrono::duration<double> time = end - start;
    cout << "Time taken to decompress: " << time.count() << " s\n";
    
    return 0;
}