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

struct Position {
    int start_reference = -1; // p
    int length = 0;           // l
    string mismatch = "";
};

// Inline to encourage inlining in performance-critical code
inline int extend_alignment(const string& Sr, const string& St, int p, int index, int k) {
    int l = k; // starting with k characters that are known to match
    while (p + l < (int)Sr.size() && index + l < (int)St.size() &&
           Sr[p + l] == St[index + l]) {
        ++l;
    }
    return l;
}

vector<Position> match_sequences(const string& Sr, const string& St, int k, int m, bool global, int offset = 0) {
    // 1. Get length of target.
    int L = St.size();
    
    // 2. Build hash table H from all k-mers in Sr.
    unordered_map<string_view, vector<int>> H;
    // Reserve capacity to avoid rehashing (number of k-mers = Sr.size()-k+1)
    H.reserve(Sr.size() - k + 1);
    for (int i = 0; i <= (int)Sr.size() - k; i++) {
        string_view kmer(&Sr[i], k);
        H[kmer].push_back(i);
    }
    
    // 3. Initialize variables.
    int index = 0;
    int prev_match_end = -1;
    vector<Position> results;
    // Heuristic: reserve results capacity based on target length.
    results.reserve((L / k) * 2);
    
    Position currentMismatch;
    // Preallocate to reduce dynamic allocation overhead.
    currentMismatch.mismatch.reserve(L);
    
    // For progress indicator if global alignment is enabled.
    int lastPrintedProgress = 0;
    
    // 4. Main loop.
    while (index < L - k + 1) {
        // Print progress every 2% if global.
        if (global) {
            int progress = (index * 100) / (L - k + 1);
            if (progress >= lastPrintedProgress + 2) {
                cout << "Progress1: " << progress << "%\n";
                lastPrintedProgress = progress;
            }
        }
        
        // Use string_view to extract k-mer from target (avoiding allocation).
        string_view kmer_prime(&St[index], k);
        
        if (H.find(kmer_prime) == H.end()) {
            // Instead of using substr to extract one character, use push_back.
            currentMismatch.mismatch.push_back(St[index]);
            index++;  // Increment index by 1.
            continue;
        } else {
            if (global) {
                bool candidateInRange = false;
                // Check if any candidate is within range.
                for (int p : H[kmer_prime]) {
                    if (prev_match_end == -1 || abs(p - prev_match_end) <= m) {
                        candidateInRange = true;
                        break;
                    }
                }
                if (!candidateInRange) {
                    currentMismatch.mismatch.push_back(St[index]);
                    index++;
                    continue;
                }
                if (!currentMismatch.mismatch.empty()) {
                    results.push_back(currentMismatch);
                    currentMismatch.mismatch.clear();
                    currentMismatch.mismatch.reserve(L - index);
                }
            } else {
                if (!currentMismatch.mismatch.empty()) {
                    results.push_back(currentMismatch);
                    currentMismatch.mismatch.clear();
                    currentMismatch.mismatch.reserve(L - index);
                }
            }
            
            // 9. Candidate selection.
            int lmax1 = 0, lmax2 = 0;
            int pn1 = 0, pn2 = 0;
            int ln1 = 0, ln2 = 0;
            for (int p : H[kmer_prime]) {
                int l = extend_alignment(Sr, St, p, index, k);
                if (global && (prev_match_end == -1 || abs(p - prev_match_end) <= m)) {
                    if (l == lmax2) {
                        if (pn2 == 0 || abs(p - prev_match_end) < abs(pn2 - prev_match_end))
                            pn2 = p;
                    } else if (l > lmax2) {
                        lmax2 = l; pn2 = p; ln2 = l;
                    }
                }
                if (l == lmax1) {
                    if (pn1 == 0 || abs(p - prev_match_end) < abs(pn1 - prev_match_end))
                        pn1 = p;
                } else if (l > lmax1) {
                    lmax1 = l; pn1 = p; ln1 = l;
                }
            }
            
            // 29-31: Choose final candidate.
            int final_p, final_l;
            if (global && pn2 != 0) {
                final_p = pn2; final_l = ln2;
            } else {
                final_p = pn1; final_l = ln1;
            }

            if (global) {
                int progress = ((index+final_l) * 100) / (L);
                if (progress >= lastPrintedProgress + 2) {
                    cout << "Progress2: " << progress << "%\n";
                    lastPrintedProgress = progress;
                }
            }
            
            // Update previous match end.
            prev_match_end = final_p + final_l - 1;
            
            // 32. Record the match.
            Position matchRes;
            matchRes.start_reference = final_p + offset; // Adjust for local subset offset.
            matchRes.length = final_l;
            matchRes.mismatch = "";
            results.push_back(matchRes);
            
            // 33. Update index.
            index += final_l;
        }
    }
    cout << index << "/" << L << " characters processed.\n";
    
    if (index < L)
        currentMismatch.mismatch.append(St, index, string::npos);
    if (!currentMismatch.mismatch.empty())
        results.push_back(currentMismatch);
    
    if (global) {
        int progress = ((index + currentMismatch.mismatch.size()) * 100) / (L);
        if (progress >= lastPrintedProgress + 2) {
            cout << "Progress3: " << progress << "%\n";
            lastPrintedProgress = progress;
        }
    }
    
    cout << "DEBUG: match_sequences() finished with " << results.size() << " results.\n";
    return results;
}

void read_genomes_from_files(const string& reference_file, 
                             const string& target_file, 
                             string& reference_genome, 
                             string& target_genome, 
                             string& target_header) {
    cout << "DEBUG: read_genomes_from_files() started.\n";
    ifstream ref_stream(reference_file);
    if (!ref_stream.is_open()) {
        cerr << "Error opening reference file: " << reference_file << "\n";
        exit(1);
    }
    string line;
    while (getline(ref_stream, line)) {
        if (line.empty() || line[0] == '>') {
            continue; // Preskoci header linije.
        }
        reference_genome += line; // Dodaj liniju u referentni genom.
    }
    ref_stream.close();
    reference_genome.erase(remove_if(reference_genome.begin(), reference_genome.end(), ::isspace), reference_genome.end());

    ifstream target_stream(target_file);
    if (!target_stream.is_open()) {
        cerr << "Error opening target file: " << target_file << "\n";
        exit(1);
    }
    bool headerFound = false;
    while (getline(target_stream, line)) {
        if (line.empty()) continue;
        if (!headerFound && line[0] == '>') {
            target_header = line;
            headerFound = true;
            continue; // do not add header line to target_genome.
        }
        target_genome += line;
    }
    target_stream.close();
    target_genome.erase(remove_if(target_genome.begin(), target_genome.end(), ::isspace), target_genome.end());
    cout << "DEBUG: Completed reading genomes.\n";
}

void delta_encode(string file_path) {
    cout << "DEBUG: delta_encode() started for file: " << file_path << "\n";

    ifstream file(file_path);
    if (!file.is_open()) {
        cerr << "Greska pri otvaranju datoteke: " << file_path << "\n";
        exit(1);
    }
    stringstream buffer;
    buffer << file.rdbuf();
    string temp_file = buffer.str();
    file.close();
    cout << "DEBUG: Read file of size " << temp_file.size() << " characters.\n";
    
    size_t search_start = 0;
    if (!temp_file.empty() && temp_file[0] == '>') {
        size_t header_end = temp_file.find('\n');
        if (header_end != string::npos) {
            size_t second_newline = temp_file.find('\n', header_end + 1);
            if (second_newline != string::npos) {
                size_t third_newline = temp_file.find('\n', second_newline + 1);
                if (third_newline != string::npos) {
                    search_start = third_newline + 1;
                }
            }
        }
    } else {
        size_t first_newline = temp_file.find('\n');
        if (first_newline != string::npos) {
            size_t second_newline = temp_file.find('\n', first_newline + 1);
            if (second_newline != string::npos) {
                search_start = second_newline + 1;
            }
        }
    }
    
    int previous_start_ref = 0;
    int token_replacements = 0;

    // Pronadi sve parove (poc,duljina) i zamijeni ih s delta kodiranim.
    while (true) {
        size_t open_pos = temp_file.find('(', search_start);
        if (open_pos == string::npos)
            break; 

        size_t start_pos = open_pos + 1;
        size_t end_pos = temp_file.find(')', start_pos);
        if (end_pos == string::npos)  // Provjera u slucaju da nema zatvorene zagrade.
            break;

        string token = temp_file.substr(start_pos, end_pos - start_pos);
        size_t comma_pos = token.find(',');
        if (comma_pos == string::npos) {
            search_start = end_pos + 1;  // Provjera ako nema zareza
            continue;
        }

        int start_ref = stoi(token.substr(0, comma_pos));
        int delta = start_ref - previous_start_ref;

        previous_start_ref = start_ref;

        string new_token = to_string(delta) + token.substr(comma_pos);

        cout << "DEBUG " << token_replacements << ": Replacing token (" << token 
             << ") with (" << new_token << ") at position " << start_pos << "\n";

        temp_file.replace(start_pos, token.size(), new_token);
        token_replacements++;

        search_start = start_pos + new_token.size();
    }
    cout << "DEBUG: Total tokens replaced: " << token_replacements << "\n";

    ofstream out_file(file_path, ofstream::out | ofstream::trunc);
    if (!out_file.is_open()) {
        cerr << "Greska pri otvaranju datoteke: " << file_path << "\n";
        exit(1);
    }
    out_file << temp_file;
    out_file.close();
    cout << "DEBUG: delta_encode() finished. Final file size: " << temp_file.size() << "\n";
}

void compress_genome_7z(const string& input_file, const string& output_file) {
    cout << "DEBUG: compress_genome() using 7z started.\n";
    string command = "7z a -mx=9 \"" + output_file + ".7z\" \"" + input_file + "\"";
    int result = system(command.c_str());

    if (result != 0) {
        cerr << "Greska prilikom komprimiranja datoteke 7-zipom: " << result << " !\n";
        exit(1);
    } else {
        cout << "Datoteka uspjesno komprimirana: " << output_file << ".7z !\n";
    }
    cout << "DEBUG: compress_genome() (reference/target) finished.\n";
}

void compress_genome(const string& reference_file, const string& target_file, const string& output_folder) {
    cout << "DEBUG: compress_genome() (reference/target) started.\n";
    // Algoritam 2 iz rada.
    string reference_genome;
    string target_genome;
    string target_header;

    read_genomes_from_files(reference_file, target_file, reference_genome, target_genome, target_header);

    string temp_file_path = output_folder + "/compressed_genome.txt";
    string output_file_path = output_folder + "/compressed_genome.txt";

    // Stvaranje privremene datoteke za spremanje rezultata.
    cout << "DEBUG: Creating output folder if necessary.\n";
    filesystem::create_directories(output_folder);
    ofstream temp_file(temp_file_path, ofstream::out | ofstream::trunc);

    if (!target_header.empty()) {
        temp_file << target_header << "\n";
    }
    
    int previous_lower_start = 0;
    int lower_start = -1;
    int lower_len = 0;
    for (int i = 0; i < target_genome.length(); ++i) {
        if (islower(target_genome[i])) {
            if (lower_len == 0)
                lower_start = i;
            ++lower_len;
        } else {
            if (lower_len != 0) {
                int delta = lower_start - previous_lower_start;
                if (lower_len == 1)
                    temp_file << delta << ",";
                else
                    temp_file << "(" << delta << "," << lower_len << ")";
                previous_lower_start = lower_start;
                lower_len = 0;
            }
        }
    }
    if (lower_len != 0) {
        int delta = lower_start - previous_lower_start;
        if (lower_len == 1)
            temp_file << delta;
        else
            temp_file << "(" << delta << "," << lower_len << ")";
    }
    temp_file << "\n,\n";
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

    int num_iterations = min(reference_segments.size(), target_segments.size());
    cout << "DEBUG: Starting segment comparisons for " << num_iterations << " iterations.\n";
    int mismatch = 0;
    for (int i = 0; i < num_iterations; ++i) {
        cout << "DEBUG: Processing segment " << i << "\n";
        string r_i = reference_segments[i];
        string t_i = target_segments[i];

        // Lokalno poravnanje (prvi prolaz, s duljinom kmera k).
        vector<Position> positions = match_sequences(r_i, t_i, k, 0, false, i*L);
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
            float mismatch_ratio = (float)count_mismatches / t_i.length();
            cout << "DEBUG: Local alignment mismatch ratio: " << mismatch_ratio << "\n";
            if (mismatch_ratio > T1 && t_i.find_first_not_of('N') != string::npos) {
                mismatch++;
                continue;
            }
            mismatch = 0;
            continue;
        }

        // Lokalno poravnanje (drugi prolaz, s duljinom k-mera k').
        positions = match_sequences(r_i, t_i, k2, 0, false, i*L);
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
            float mismatch_ratio = (float)count_mismatches / t_i.length();
            cout << "DEBUG: Second pass mismatch ratio: " << mismatch_ratio << "\n";            
            if (mismatch_ratio > T1  && t_i.find_first_not_of('N') != string::npos) {
                mismatch++;
                continue;
            }
            mismatch = 0;
            continue;
        }

        if (t_i.find_first_not_of('N') != string::npos) {
            // Ako je poravnanje neuspjesno, povecaj broj neuspjesnih pokusaja.
            mismatch++;
        } else {
            // Ako je cijeli segment N, resetiraj broj neuspjesnih pokusaja.
            mismatch = 0;
        }

        if (mismatch > T2) {
            // Broj neuspjeha je prevelik, zaustavlja se lokalno poravanje i pocinje globalno
            cout << "DEBUG: Mismatch threshold exceeded. Switching to global alignment.\n";

            // Obrisi prethodni sadrzaj temp_file.
            temp_file.close();
            temp_file.open(temp_file_path, ofstream::out | ofstream::trunc);
            temp_file.close();
            
            local = false;
            break;
        }
    }

    if (local && target_segments.size() > num_iterations) {
        cout << "DEBUG: Adding remaining target segments as mismatches.\n";
        for (size_t j = num_iterations; j < target_segments.size(); ++j) {
            temp_file << target_segments[j];
        }
    }

    // Ako lokalno poravnanje nije uspjelo
    if (!local) {
        reference_genome.clear();
        target_genome.clear();
        read_genomes_from_files(reference_file, target_file, reference_genome, target_genome, target_header);

        temp_file.open(temp_file_path, ofstream::out | ofstream::trunc);

        if (!target_header.empty()) {
            temp_file << target_header << "\n";
        }

        int previous_lower_start = 0;
        int lower_start = -1;
        int lower_len = 0;
        for (int i = 0; i < target_genome.length(); ++i) {
            if (islower(target_genome[i])) {
                if (lower_len == 0)
                    lower_start = i;
                ++lower_len;
            } else {
                if (lower_len != 0) {
                    int delta = lower_start - previous_lower_start;
                    if (lower_len == 1)
                        temp_file << delta << ",";
                    else
                        temp_file << "(" << delta << "," << lower_len << ")";
                    previous_lower_start = lower_start;
                    lower_len = 0;
                }
            }
        }
        if (lower_len != 0) {
            int delta = lower_start - previous_lower_start;
            if (lower_len == 1)
                temp_file << delta;
            else
                temp_file << "(" << delta << "," << lower_len << ")";
        }
        temp_file << "\n";
        transform(target_genome.begin(), target_genome.end(), target_genome.begin(), ::toupper);
        transform(reference_genome.begin(), reference_genome.end(), reference_genome.begin(), ::toupper);
        cout << "DEBUG: Recorded lowercase indices in temp_file.\n";

        int previous_n_start = 0;
        int n_start = -1;
        int n_len = 0;
        for (int i = 0; i < target_genome.length(); ++i) {
            if (target_genome[i] == 'N') {
                if (n_start == -1)
                    n_start = i;
                ++n_len;
            } else {
                if (n_len != 0) {
                    int delta = n_start - previous_n_start;
                    if (n_len == 1)
                        temp_file << delta << ",";
                    else
                        temp_file << "(" << delta << "," << n_len << ")";
                    previous_n_start = n_start;
                    n_start = -1;
                    n_len = 0;
                }
            }
        }
        if (n_len != 0) {
            int delta = n_start - previous_n_start;
            if (n_len == 1)
                temp_file << delta;
            else
                temp_file << "(" << delta << "," << n_len << ")";
        }
        temp_file << "\n";
        target_genome.erase(remove(target_genome.begin(), target_genome.end(), 'N'), target_genome.end());
        reference_genome.erase(remove(reference_genome.begin(), reference_genome.end(), 'N'), reference_genome.end());
        cout << "DEBUG: Recorded 'N' characters in temp_file.\n";

        // Globalno poravnanje (s duljinom k-mera k).
        vector<Position> positions = match_sequences(reference_genome, target_genome, k, m, true);

        // Spremi rezultate u temp_file.
        for (Position pos : positions) {
            //cout << " (" << pos.start_reference << ", " << pos.length << ")\n";
            //cout << pos.mismatch;
            //cout << "\n\n";
            if (pos.mismatch == "") {
                temp_file << "(" << pos.start_reference << "," << pos.length << ")";
            } else {
                temp_file <<  pos.mismatch;
            }
        }
    }

    temp_file.close();

    // Delta kodiranje rezultata.
    delta_encode(temp_file_path);

    compress_genome_7z(temp_file_path, output_file_path);
} 

int main(int argc, char* argv[]) {
    cout << "DEBUG: main() started.\n";
    cout << "Received " << argc << " arguments.\n";
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <reference_file> <target_file> <output_folder>\n";
        return 1;
    }
    try {
        if (!filesystem::exists(argv[3])) {
            filesystem::create_directory(argv[3]);
        }
        cout << "Successfully created output folder: " << argv[3] << "\n";
        cout << "Starting compression.";
        auto start = chrono::high_resolution_clock::now();
        compress_genome(argv[1], argv[2], argv[3]);
        auto end = chrono::high_resolution_clock::now();

        chrono::duration<double> time = end - start;
        cout << "Time taken to compress: " << time.count() << " s\n";

    } catch (const std::exception& ex) {
        cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
    cout << "DEBUG: main() finished successfully.\n";
    return 0;
}
