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

struct Position {
    int start_reference = -1; // p
    int length = 0;           // l
    string mismatch = "";
};

// Inline to encourage inlining in performance-critical code
inline int extend_alignment(const string& Sr, const string& St, int p, int index, int k) {
    int l = k; // starting with k characters that are known to match
    while (p + l + 1 < (int)Sr.size() && index + l + 1 < (int)St.size() &&
           Sr[p + l + 1] == St[index + l + 1]) {
        ++l;
    }
    return l;
}

vector<Position> match_sequences(const string& Sr, const string& St, int k, int m, bool global, int ind = 0) {
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
    int index = ind;
    int prev_match_end = -1;
    vector<Position> results;
    // Heuristic: reserve results capacity based on target length.
    results.reserve(L / k);
    
    Position currentMismatch;
    // Preallocate to reduce dynamic allocation overhead.
    currentMismatch.mismatch.reserve(L);
    
    // For progress indicator if global alignment is enabled.
    int lastPrintedProgress = 0;
    
    // 4. Main loop.
    while (index < L - k + 1) {
        // Print progress every 2% if global.
        if (global) {
            int progress = (index * 100) / L;
            if (progress >= lastPrintedProgress + 2) {
                cout << "Progress: " << progress << "%\n";
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
            
            // Update previous match end.
            prev_match_end = final_p + final_l - 1;
            
            // 32. Record the match.
            Position matchRes;
            matchRes.start_reference = final_p;
            matchRes.length = final_l;
            matchRes.mismatch = "";
            results.push_back(matchRes);
            
            // 33. Update index.
            index += final_l;
        }
    }
    
    // Append any remaining characters (using append overload avoids temporary substring).
    if (index < L)
        currentMismatch.mismatch.append(St, index, string::npos);
    if (!currentMismatch.mismatch.empty())
        results.push_back(currentMismatch);
    
    cout << "DEBUG: match_sequences() finished with " << results.size() << " results.\n";
    return results;
}

void read_genomes_from_files(const string& reference_file, const string& target_file, string& reference_genome, string& target_genome) {
    // Ucitavanje genoma iz datoteka.
    cout << "DEBUG: read_genomes_from_files() started.\n";
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
    cout << "DEBUG: Completed reading genomes.\n";
}

void delta_encode(string file_path) {
    // Read the entire file into a string.
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
    size_t search_start = 0;

    // Find all tokens in the format (start_reference,length) and replace start_reference with delta.
    while (true) {
        size_t open_pos = temp_file.find('(', search_start);
        if (open_pos == string::npos)
            break;  // No more tokens.

        size_t start_pos = open_pos + 1;
        size_t end_pos = temp_file.find(')', start_pos);
        if (end_pos == string::npos)
            break;  // Malformed token; no matching ')'.

        // Extract the token between the parentheses.
        string token = temp_file.substr(start_pos, end_pos - start_pos);
        size_t comma_pos = token.find(',');
        if (comma_pos == string::npos) {
            // If no comma, skip this token.
            search_start = end_pos + 1;
            continue;
        }
        // Parse the original start reference.
        int start_ref = stoi(token.substr(0, comma_pos));
        int delta = start_ref - previous_start_ref;
        // Update previous_start_ref with the absolute value.
        previous_start_ref = start_ref;

        // Build the new token content: replace the original number with the delta,
        // and keep the rest (from the comma onward) unchanged.
        string new_token = to_string(delta) + token.substr(comma_pos);
        // Replace the existing token content with the new token.
        temp_file.replace(start_pos, token.size(), new_token);

        // Move search_start to after the closing parenthesis.
        search_start = end_pos;
    }

    // Write the modified content back into the file.
    ofstream out_file(file_path, ofstream::out | ofstream::trunc);
    if (!out_file.is_open()) {
        cerr << "Greska pri otvaranju datoteke: " << file_path << "\n";
        exit(1);
    }
    out_file << temp_file;
    out_file.close();
    cout << "DEBUG: delta_encode() finished.\n";
}

void compress_genome_7z(const string& input_file, const string& output_file) {
    cout << "DEBUG: compress_genome() using 7z started.\n";
    // Komprimiraj datoteku (7-zip-om, kao u radu).
    string command = "7z a -mx=9 " + output_file + ".7z " + input_file;
    int result = system(command.c_str());

    if (result != 0) {
        cerr << "Greska prilikom komprimiranja datoteke 7-zipom: " << result << " !\n";
        exit(1);
    } else {
        cout << "Datoteka uspjesno komprimirana: " << output_file << ".7z !\n";
    }
    cout << "DEBUG: compress_genome() using 7z finished.\n";
}

void compress_genome(const string& reference_file, const string& target_file, const string& output_folder) {
    cout << "DEBUG: compress_genome() (reference/target) started.\n";
    // Algoritam 2 iz rada.
    string reference_genome;
    string target_genome;

    read_genomes_from_files(reference_file, target_file, reference_genome, target_genome);
    
    string temp_file_path = output_folder + "/compressed_genome.txt";
    string output_file_path = output_folder + "/compressed_genome.txt";

    // Stvaranje privremene datoteke za spremanje rezultata.
    cout << "DEBUG: Creating output folder if necessary.\n";
    filesystem::create_directories(output_folder);
    ofstream temp_file(temp_file_path, ofstream::out | ofstream::trunc);

    // Prebaci sve znakove u velika slova i zapamti indekse malih slova u temp_file
    cout << "DEBUG: Writing lowercase indices from reference_genome.\n";
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

    // After the main for loop of segment comparisons:
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
    // delta_encode(temp_file_path);
    cout << "DEBUG: delta_encode() started.\n";

    compress_genome_7z(temp_file_path, output_file_path);
    cout << "DEBUG: compress_genome() (reference/target) finished.\n";
} 

// Main function.
int main(int argc, char* argv[]) {
    cout << "DEBUG: main() started.\n";
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
    cout << "DEBUG: main() finished successfully.\n";
    return 0;
}
