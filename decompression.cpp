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

void decompress_genome(string compressed_file_path, string reference_genome_path, string output_folder) {
    filesystem::path compressedPath(compressed_file_path);
    filesystem::path outputPath(output_folder);

    // Stvori output datoteku ako ne postoji.
    if (!filesystem::exists(outputPath)) {
        filesystem::create_directory(outputPath);
    }

    std::string command = "7z e \"" + compressedPath.string() + "\" -o\"" + outputFolderPath.string() + "\" -y";
    int result = system(command.c_str());
    if (result != 0) {
        cerr << "Greska prilikom dekompresije datoteke 7-zipom: " << compressed_file_path << "\n";
        exit(1);
    } else {
        cout << "Datoteka uspjesno dekompriminirana: " << output_folder << "\n";
    }

    // Otvori dekomprimiranu datoteku.
    string decompressed_file_path = output_folder + "/" + compressedPath.stem().string();
    cout << "Pronadena dekomprimirana datoteka: " << decompressed_file_path << "\n";

    // Otvori referentni genom.
    ifstream reference_file(reference_genome_path);
    if (!reference_file.is_open()) {
        cerr << "Greska pri otvaranju referentne datoteke: " << reference_genome_path << "\n";
        exit(1);
    }
    string reference_genome;
    getline(reference_file, reference_genome);
    reference_file.close();

    // Otvori dekomprimiranu datoteku.
    ifstream decompressed_file(decompressed_file_path);
    if (!decompressed_file.is_open()) {
        cerr << "Greska pri otvaranju dekomprimirane datoteke: " << decompressed_file_path << "\n";
        exit(1);
    }
    string target_genome;
    getline(decompressed_file, target_genome);
    decompressed_file.close();

    
}

int main(int argc, char* argv[]) {
    cout << "Received " << argc - 1 << " arguments.\n";
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <compressed_file> <output_folder>\n";
        return 1;
    }
    try {
        // Ensure output folder exists.
        if (!filesystem::exists(argv[3])) {
            filesystem::create_directory(argv[3]);
        }
        cout << "Successfully created output folder: " << argv[2] << "\n";

        decompress_genome(argv[1], argv[2]);

    } catch (const std::exception& ex) {
        cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
    return 0;

    return 0;
}