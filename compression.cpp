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
    int start_reference = -1; // p
    int length = 0;          // l
    string mismatch = "";
};

// Function: extend_alignment
// Inputs:
//    Sr: Reference segment/sequence.
//    St: Target segment/sequence.
//    p:  Starting position in Sr.
//    index: Starting position in St.
//    k:  Length of the k-mer.
// Returns: Length of the extended alignment.
// This function extends the alignment from position p in Sr and index in St,
// checking for matches until a mismatch is found or the end of either sequence is reached.
int extend_alignment(const string& Sr, const string& St, int p, int index, int k) {
    int l = k; // starting with k characters that are known to match
    while (p + l < (int)Sr.size() && index + l < (int)St.size() &&
           Sr[p + l] == St[index + l]) {
        l++;
    }
    return l;
}

// Function: match_sequences
// Inputs:
//    Sr: Reference segment/sequence.
//    St: Target segment/sequence.
//    k: Length of the k-mer.
//    m: Limit of search range (for global matching distance)
//    global: Boolean flag, if true enables global matching.
vector<Position> match_sequences(const string& Sr, const string& St, int k, int m, bool global) {
    // 1. Get length L of target St.
    int L = St.size();
    
    // 2. Build hash table H from all k-mers in Sr.
    //    This maps each k-mer to a vector of its starting positions in Sr.
    unordered_map<string, vector<int>> H;
    for (int i = 0; i <= (int)Sr.size() - k; i++) {
        string kmer = Sr.substr(i, k);
        H[kmer].push_back(i);
    }
    
    // 3. Initialize index, previous match end, results,
    //    and an accumulator for mismatched characters.
    int index = 0;
    int prev_match_end = -1;
    vector<Position> results;
    Position currentMismatch;   // Will accumulate mismatched characters; start_reference remains -1
    
    // 4. Main loop: while index < L - k + 1.
    while (index < L - k + 1) {
        // 4.1 Extract the k-mer starting at index in St.
        string kmer_prime = St.substr(index, k);
        
        // 4.2 If no match for kmer_prime is found in H:
        if (H.find(kmer_prime) == H.end()) {
            // Accumulate the mismatched character.
            currentMismatch.mismatch += St.substr(index, 1);
            index++;  // Increment index by 1.
            continue;
        } else {
            if (global) {
                bool candidateInRange = false;
                // Check if there is at least one candidate p that is in range.
                for (int p : H[kmer_prime]) {
                    if (prev_match_end == -1 || abs(p - prev_match_end) <= m) {
                        candidateInRange = true;
                        break;
                    }
                }
                if (!candidateInRange) {
                    index++; // no candidate in range, continue searching
                    continue;
                }
                // Push current mismatches to clear them.
                if (!currentMismatch.mismatch.empty()) {
                    results.push_back(currentMismatch);
                    currentMismatch.mismatch = "";
                }
            } else {
                // For local matching, push mismatches only if there's something accumulated.
                if (!currentMismatch.mismatch.empty()) {
                    results.push_back(currentMismatch);
                    currentMismatch.mismatch = "";
                }
            }
            
            // 9. Initialize variables for best candidates.
            int lmax1 = 0, lmax2 = 0;
            int pn1 = 0, pn2 = 0;
            int ln1 = 0, ln2 = 0;
            
            // 10. For each occurrence of kmer_prime in H.
            for (int p : H[kmer_prime]) {
                // 11. Extend the alignment from position p in Sr and index in St.
                int l = extend_alignment(Sr, St, p, index, k);
                
                // 12. If global is enabled and the candidate's starting position p is 
                //     within m positions of the previous match's end:
                if (global && (prev_match_end == -1 || abs(p - prev_match_end) <= m)) {
                    // 13-20: Update global candidate.
                    if (l == lmax2) {
                        // 14. If this candidate's p is nearer (closer) to prev_match_end than the current candidate:
                        if (pn2 == 0 || abs(p - prev_match_end) < abs(pn2 - prev_match_end)) {
                            // 15. Update pn2 with this candidate.
                            pn2 = p;
                        }
                    } else if (l > lmax2) {
                        // 17-18. New longer match: update lmax2, pn2, and ln2.
                        lmax2 = l;
                        pn2 = p;
                        ln2 = l;
                    }
                }
                
                // 21. Otherwise (or regardless) update the best local candidate.
                if (l == lmax1) {
                    // 22. If p of the candidate is nearer to prev_match_end than the current local candidate:
                    if (pn1 == 0 || abs(p - prev_match_end) < abs(pn1 - prev_match_end)) {
                        // 23. Update pn1.
                        pn1 = p;
                    }
                } else if (l > lmax1) {
                    // 25-26. New longer match: update lmax1, pn1, and ln1.
                    lmax1 = l;
                    pn1 = p;
                    ln1 = l;
                }
            }
            
            // 29-31: Choose final candidate based on the global flag.
            int final_p, final_l;
            if (global && pn2 != 0) {
                final_p = pn2;
                final_l = ln2;
            } else {
                final_p = pn1;
                final_l = ln1;
            }
            
            // Update previous match end (for global evaluation).
            prev_match_end = final_p + final_l - 1;
            
            // 32. Record the match.
            Position matchRes;
            matchRes.start_reference = final_p;
            matchRes.length = final_l;
            // For a match, we leave mismatch as an empty string.
            matchRes.mismatch = "";
            results.push_back(matchRes);
            
            // 33. Update index.
            index = index + final_l + 1;
        }
    }
    
    // After the loop, if any mismatches remain, add them as a final entry.
    if (!currentMismatch.mismatch.empty()) {
        results.push_back(currentMismatch);
    }
    
    // 35. Return the results.
    return results;
}

void read_genomes_from_files(const string& reference_file, const string& target_file, string& reference_genome, string& target_genome) {
    // Ucitavanje genoma iz datoteka.
    ifstream ref_stream(reference_file);
    if (!ref_stream.is_open()) {
        cerr << "Error opening reference file: " << reference_file << "\n";
        exit(1);
    }

    // Obrisi linije koje pocinju s '>' (header linije).
    string line;
    while (getline(ref_stream, line)) {
        if (line.empty() || line[0] == '>') {
            continue; // Preskoci header linije.
        }
        reference_genome += line; // Dodaj liniju u referentni genom.
    }
    ref_stream.close();
    
    // Obrisi znakove novog reda (i ostale praznine ako postoje).
    reference_genome.erase(remove_if(reference_genome.begin(), reference_genome.end(), ::isspace), reference_genome.end());

    ifstream target_stream(target_file);
    if (!target_stream.is_open()) {
        cerr << "Error opening target file: " << target_file << "\n";
        exit(1);
    }

    while (getline(target_stream, line)) {
        if (line.empty() || line[0] == '>') {
            continue; // Preskoci header linije.
        }
        target_genome += line; // Dodaj liniju u target genom.
    }
    target_stream.close();

    // Obrisi znakove novog reda (i ostale praznine ako postoje).
    target_genome.erase(remove_if(target_genome.begin(), target_genome.end(), ::isspace), target_genome.end());
}

void delta_encode(string file_path) {
    // Izmijeni pozicije iz datoteke na nacin da se koristi delta enkodiranje.
    // Svaka pocetna pozicija (osim prve) je zabiljezena kao razlika od prethodne.
    ifstream file(file_path);
    if (!file.is_open()) {
        cerr << "Greska pri otvaranju datoteke: " << file_path << "\n";
        exit(1);
    }
    stringstream buffer;
    buffer << file.rdbuf();
    string temp_file = buffer.str();
    file.close();

    int previous_start_ref = 0;

    // TODO: provjeriti radi li dobro.
    for (int i = 1; i < temp_file.size(); ++i) {
        if (temp_file[i] == '(') {
            size_t start_pos = i + 1;
            size_t end_pos = temp_file.find(')', start_pos);
            string position = temp_file.substr(start_pos, end_pos - start_pos);
            size_t comma_pos = position.find(',');
            int start_ref = stoi(position.substr(0, comma_pos));
            int length = stoi(position.substr(comma_pos + 1));
            // Delta enkodiranje.
            if (i > 0) {
                start_ref -= previous_start_ref;
            }
            previous_start_ref = start_ref;
            temp_file.replace(start_pos, end_pos - start_pos + 1, to_string(start_ref) + "," + to_string(length));
        }
    }

    // Spremi izmijenjeni sadrzaj u datoteku.
    ofstream out_file(file_path, ofstream::out | ofstream::trunc);
    if (!out_file.is_open()) {
        cerr << "Greska pri otvaranju datoteke: " << file_path << "\n";
        exit(1);
    }
    out_file << temp_file;
    out_file.close();
}

void compress_genome_7z(const string& input_file, const string& output_file) {
    // Komprimiraj datoteku (7-zip-om, kao u radu).
    string command = "7z a -mx=9 " + output_file + ".7z " + input_file;
    int result = system(command.c_str());

    if (result != 0) {
        cerr << "Greska prilikom komprimiranja datoteke 7-zipom: " << result << " !\n";
        exit(1);
    } else {
        cout << "Datoteka uspjesno komprimirana: " << output_file << ".7z !\n";
    }
}

void compress_genome(const string& reference_file, const string& target_file, const string& output_folder) {
    // Algoritam 2 iz rada.
    string reference_genome;
    string target_genome;

    read_genomes_from_files(reference_file, target_file, reference_genome, target_genome);
    
    string temp_file_path = output_folder + "/compressed_genome.txt";
    string output_file_path = output_folder + "/compressed_genome.txt";

    // Stvaranje privremene datoteke za spremanje rezultata.
    filesystem::create_directories(output_folder);
    ofstream temp_file(temp_file_path, ofstream::out | ofstream::trunc);

    // Prebaci sve znakove u velika slova i zapamti indekse malih slova u temp_file
    for (int i = 0; i < target_genome.length(); ++i) {
        if (islower(target_genome[i])) {
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
        string r_i = reference_segments[i];
        string t_i = target_segments[i];

        // Lokalno poravnanje (prvi prolaz, s duljinom kmera k).
        vector<Position> positions = match_sequences(r_i, t_i, k, 0, false);
        if (positions.size() == 1 && positions[0].mismatch == "" || positions.size() > 1) {
            // Poravnanje je uspjesno.
            int count_mismatches = 0;
            int total_length = 0;
            for (Position pos : positions) {
                if (pos.mismatch == "") {
                    temp_file << "(" << pos.start_reference << "," << pos.length << ")";
                    total_length += pos.length;
                }
                else {
                    temp_file << pos.mismatch;
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
        positions = match_sequences(r_i, t_i, k, 0, false);
        if (positions.size() == 1 && positions[0].mismatch == "" || positions.size() > 1) {
            // Poravnanje je uspjesno.
            int count_mismatches = 0;
            int total_length = 0;
            for (Position pos : positions) {
                if (pos.mismatch == "") {
                    temp_file << "(" << pos.start_reference << "," << pos.length << ")";
                    total_length += pos.length;
                }
                else {
                    temp_file << pos.mismatch;
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
            temp_file.close();
            
            local = false;
            break;
        }
    }

    // Ako lokalno poravnanje nije uspjelo
    if (!local) {
        reference_genome.clear();
        target_genome.clear();
        read_genomes_from_files(reference_file, target_file, reference_genome, target_genome);

        temp_file.open(temp_file_path, ofstream::out | ofstream::trunc);
        // Prebaci sve znakove u velika slova i zapamti indekse malih slova u temp_file
        temp_file << "lower: ";
        for (int i = 0; i < target_genome.length(); ++i) {
            if (islower(target_genome[i])) {
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
        vector<Position> positions = match_sequences(reference_genome, target_genome, k, m, true);

        // Spremi rezultate u temp_file.
        for (Position pos : positions) {
            if (pos.mismatch == "") {
                temp_file << "(" << pos.start_reference << "," << pos.length << ")";
            } else {
                temp_file<<  pos.mismatch;
            }
        }
    }

    temp_file.close();

    // Delta kodiranje rezultata.
    //delta_encode(temp_file_path);

    compress_genome_7z(temp_file_path, output_file_path);
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
