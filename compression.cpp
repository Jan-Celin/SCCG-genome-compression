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

struct Position {
    string mismatch = "";
    int start_reference = -1; // p
    int length = 0;          // l
};

void read_genomes_from_files(const string& reference_file, const string& target_file, string& reference_genome, string& target_genome) {
    // Ucitavanje genoma iz datoteka.
    ifstream ref_stream(reference_file);
    if (!ref_stream.is_open()) {
        cerr << "Error opening reference file: " << reference_file << "\n";
        exit(1);
    }

    stringstream ref_buffer;
    ref_buffer << ref_stream.rdbuf();
    reference_genome = ref_buffer.str();
    ref_stream.close();
    
    // Obrisi znakove novog reda (i ostale praznine ako postoje).
    reference_genome.erase(remove_if(reference_genome.begin(), reference_genome.end(), ::isspace), reference_genome.end());

    ifstream target_stream(target_file);
    if (!target_stream.is_open()) {
        cerr << "Error opening target file: " << target_file << "\n";
        exit(1);
    }

    stringstream target_buffer;
    target_buffer << target_stream.rdbuf();
    target_genome = target_buffer.str();
    target_stream.close();

    // Obrisi znakove novog reda (i ostale praznine ako postoje).
    target_genome.erase(remove_if(target_genome.begin(), target_genome.end(), ::isspace), target_genome.end());
}

List<Position> matching(string Sr, string St, int k, int m, bool global){
    return List<Position>();
}

void delta_encode(string file_path) {
    // Izmijeni pozicije iz datoteke na nacin da se koristi delta enkodiranje.
    // Svaka pocetna pozicija (osim prve) je zabiljezena kao razlika od prethodne.
    ifstream file(file_path);
    // TODO: provjeriti radi li dobro.
    for (int i = 1; i < temp_file.size(); ++i) {
        if (file[i] == '(') {
            size_t start_pos = i + 1;
            size_t end_pos = temp_file.find(')', start_pos);
            if (end_pos != string::npos) {
                string position = file.substr(start_pos, end_pos - start_pos);
                size_t comma_pos = position.find(',');
                if (comma_pos != string::npos) {
                    int start_ref = stoi(position.substr(0, comma_pos));
                    int length = stoi(position.substr(comma_pos + 1));
                    // Delta enkodiranje.
                    if (i > 0) {
                        start_ref -= previous_start_ref;
                    }
                    previous_start_ref = start_ref;
                    file.replace(start_pos, end_pos - start_pos + 1, to_string(start_ref) + "," + to_string(length));
                }
            }
        }
    }
}

void compress_genome(const string& input_file, const string& output_file) {
    // Komprimiraj datoteku (PPMd algoritmom od 7-zipa).
    string command = "7z a -mx=9 " + output_file + ".7z " + output_file;
    int result = system(command.c_str());

    if (result != 0) {
        cerr << "Greska prilikom komprimiranja datoteke 7-zipom: " << result << " !\n";
        exit(1);
    } else {
        cout << "Datoteka uspjesno komprimirana: ": << output_file << ".7z !\n";
    }
}

void compress_genome(const string& reference_file, const string& target_file, const string& output_folder) {
    // Algoritam 2 iz rada.
    string reference_genome;
    string target_genome;

    read_genomes_from_files(reference_file, target_file, reference_genome, target_genome);
    
    string temp_file_path = output_folder + "temp_file.txt";
    string output_file_path = output_folder + "compressed_genome.txt";

    // Stvaranje privremene datoteke za spremanje rezultata.
    filesystem::create_directories(output_folder);
    ofstream temp_file(temp_file_path, ofstream::out | ofstream::trunc);

    // Prebaci sve znakove u velika slova i zapamti indekse malih slova u temp_file
    for (int i = 0; i < reference_genome.length(); ++i) {
        if (islower(reference_genome[i])) {
            temp_file << i << ",";
        }
    }
    temp_file << "\n";
    transform(reference_genome.begin(), reference_genome.end(), reference_genome.begin(), ::toupper);
    transform(target_genome.begin(), target_genome.end(), target_genome.begin(), ::toupper);

    // Pocetni parametri algoritma:
    int k = 14;  // primarna duljina k-mera
    int k2 = 10;  // sekundarna duljina k-mera
    int L = 1000;  // duljina segmenta
    int m = 100;  // udaljenost pretrazivanja preklapanja kod globalnog poravnanja
    float T1 = 0.5;  // prag za odredivanje kvalitete preklapanja;
    int T2 = 4;  // koliko segmenata zaredom morea biti lose kvalitete
    bool local = true;  // lokalno poravnanje

    // Dijeljenje reference_genome i target_genome na segmente duljine L.
    vector<string> reference_segments;
    vector<string> target_segments;

    for (size_t i = 0; i < reference_genome.length(); i += L) {
        reference_segments.push_back(reference_genome.substr(i, L));
    }
    for (size_t i = 0; i < target_genome.length(); i += L) {
        target_segments.push_back(target_genome.substr(i, L));
    }

    // Prodi kroz sve parove segmenata reference i targeta.
    int num_iterations = 0;
    if (reference_segments.size() < target_segments.size()) {
        num_iterations = reference_segments.size();
    } else {
        num_iterations = target_segments.size();
    }

    int mismatch = 0;
    for (int i = 0; i < num_iterations; ++i) {
        r_i = reference_segments[i];
        t_i = target_segments[i];

        // Lokalno poravnanje (prvi prolaz, s duljinom kmera k).
        List<Position> positions = r_i, t_i, k, 0, false
        if (positions.length() == 1 && positions[0].mismatch == "" || positions.length() > 1) {
            // Poravnanje je uspjesno.
            int count_mismatches = 0;
            int total_length = 0;
            for (Position pos : positions) {
                if (pos.mismatch == "") {
                    temp_file << " (" << pos.start_reference << ", " << pos.length << ") ";
                    total_length += pos.length;
                }
                else {
                    temp_file << " " << pos.mismatch " ";
                    count_mismatches += pos.mismatch.length();
                }
            }
            // Provjeri kvalitetu poravnanja
            float mismatch_ratio = (float)count_mismatches / total_length;
            if (mismatch_ratio > T1 && t_i.find_first_not_of('N') != string::npos) {
                mismatch++;
            }
            continue;
        }

        // Lokalno poravnanje (drugi prolaz, s duljinom k-mera k').
        List<Position> positions = r_i, t_i, k, 0, false
        if (positions.length() == 1 && positions[0].mismatch == "" || positions.length() > 1) {
            // Poravnanje je uspjesno.
            int count_mismatches = 0;
            int total_length = 0;
            for (Position pos : positions) {
                if (pos.mismatch == "") {
                    temp_file << " (" << pos.start_reference << ", " << pos.length << ") ";
                    total_length += pos.length;
                }
                else {
                    temp_file << " " << pos.mismatch " ";
                    count_mismatches += pos.mismatch.length();
                }
            }
            // Provjeri kvalitetu poravnanja
            float mismatch_ratio = (float)count_mismatches / total_length;
            if (mismatch_ratio > T1  && t_i.find_first_not_of('N') != string::npos) {
                mismatch++;
            }
            continue;
        }

        if (mismatch > T2) {
            // Broj neuspjeha je prevelik, zaustavlja se lokalno poravanje i pocinje globalno
            cout << "Lokalno poravnanje neuspjesno, zapocinje globalno!\n";

            // Obrisi prethodni sadrzaj temp_file.
            temp_file.close();
            temp_file.open("temp_file.txt", ofstream::out | ofstream::trunc);
            temp_file.close()
            
            local = false;
            break;
        }
    }

    // Ako lokalno poravnanje nije uspjelo
    if (!local) {
        temp_file.open(output_file, ofstream::out | ofstream::trunc);
        // Prebaci sve znakove u velika slova i zapamti indekse malih slova u temp_file
        temp_file << "lower: ";
        for (int i = 0; i < reference_genome.length(); ++i) {
            if (islower(reference_genome[i])) {
                temp_file << i << ",";
            }
        }
        temp_file << "\n";
        transform(reference_genome.begin(), reference_genome.end(), reference_genome.begin(), ::toupper);
        transform(target_genome.begin(), target_genome.end(), target_genome.begin(), ::toupper);

        // Obrisi sve N znakove u target_genome, zapamti njihove indekse u temp_file.
        temp_file << "N: ";
        for (size_t i = 0; i < target_genome.length(); ++i) {
            if (target_genome[i] == 'N') {
                temp_file << i << ",";
            }
        }
        temp_file << "\n";
        target_genome.erase(remove(target_genome.begin(), target_genome.end(), 'N'), target_genome.end());

        // Globalno poravnanje (s duljinom k-mera k).
        List<Position> positions = matching(reference_genome, target_genome, k, m, true);

        // Spremi rezultate u temp_file.
        for (Position pos : positions) {
            if (pos.mismatch == "") {
                temp_file << " (" << pos.start_reference << ", " << pos.length << ") ";
            } else {
                temp_file << " " << pos.mismatch << " ";
            }
        }
    }

    temp_file.close();

    // Delta kodiranje rezultata.
    delta_encode(temp_file_path);

    compress_genome(temp_file_path, output_file_path);
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
