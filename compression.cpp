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
            // Before processing a match, if we've accumulated mismatches, add them to results.
            if (!currentMismatch.mismatch.empty()) {
                results.push_back(currentMismatch);
                currentMismatch.mismatch = "";
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
