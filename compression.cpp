#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <climits>
#include <fstream>
#include <sstream>
#include <functional>
#include <cmath>
#include <set>

using namespace std;


class Encoder {
public:
    std::string reference, target;
    int kmer_length;
    double T1;
    double T2;

    class Kmer {
    public:
        std::string kmer;
        int kmerstart;
    };
    
    std::unordered_map<int, std::vector<Kmer>> kmer_map;

    class Position {
    public:
        int start_reference;
        int end_reference;
        int start_target;
        int end_target;
    };

    void create_hash_map_L(
        const std::string& reference, 
        int kmer_length
    ) {
        for (int i=0; i < reference.length() - kmer_length + 1; i++) {
            std::string kmer = reference.substr(i, kmer_length);
            int hash_value = std::hash<std::string>{}(kmer);
            Kmer element;
            element.kmer = kmer;
            element.kmerstart = i;
            
            if (kmer_map.find(hash_value) == kmer_map.end()) {
                kmer_map[hash_value] = std::vector<Kmer>();
                kmer_map[hash_value].push_back(element);
            }
            else {
                kmer_map[hash_value].push_back(element);
            }

            // TODO-JC: Preskoci n znakove ( provjeriti je li tocno )
            if (reference[i] == 'N' || reference[i] == 'n') {
                while (i < reference.length() - kmer_length + 1 && (reference[i] == 'N' || reference[i] == 'n')) {
                    i++;
                }
            }
        }
    }

    std::vector<Encoder::Position> Lmatch(
        const std::string& reference, 
        const std::string& target, 
        int kmer_length, 
        double T1, 
        double T2,
        int i = 0
    ) {
        create_hash_map_L(reference, kmer_length);
        std::vector<Position> pos_list;
        std::set<int> used_reference_positions;  // Track used reference indices

        while (i <= (int)target.length() - kmer_length) {
            std::string kmer = target.substr(i, kmer_length);
            int hash_value = std::hash<std::string>{}(kmer);

            if (kmer_map.find(hash_value) == kmer_map.end()) {
                ++i;
                continue;
            }

            const std::vector<Kmer>& kmer_list = kmer_map[hash_value];
            int max_increment = -1;
            int best_start_ref = -1;

            for (const auto& new_kmer : kmer_list) {
                if (new_kmer.kmer != kmer) continue;

                int kmer_start = new_kmer.kmerstart;
                int end_reference = kmer_start + kmer_length - 1;
                int end_target = i + kmer_length - 1;

                int max_extend = std::min((int)reference.size() - 1 - end_reference, (int)target.size() - 1 - end_target);
                int increment = 0;
                for (; increment < max_extend; ++increment) {
                    if (reference[end_reference + 1 + increment] != target[end_target + 1 + increment]) {
                        break;
                    }
                }

                // Check if this reference region overlaps any already-used positions
                bool overlaps = false;
                for (int p = kmer_start; p <= kmer_start + kmer_length + increment - 1; ++p) {
                    if (used_reference_positions.count(p)) {
                        overlaps = true;
                        break;
                    }
                }

                if (overlaps) continue;

                if (increment > max_increment || (increment == max_increment && kmer_start < best_start_ref)) {
                    max_increment = increment;
                    best_start_ref = kmer_start;
                }
            }

            if (best_start_ref == -1) {
                ++i;
                continue;
            }

            Position pos;
            pos.start_reference = best_start_ref;
            pos.end_reference = best_start_ref + kmer_length + max_increment - 1;
            pos.start_target = i;
            pos.end_target = i + kmer_length + max_increment - 1;
            pos_list.push_back(pos);

            // Mark used reference positions
            for (int p = pos.start_reference; p <= pos.end_reference; ++p) {
                used_reference_positions.insert(p);
            }

            i = pos.end_target + 1;  // Advance past the matched region in target
        }

        return pos_list;
    }

    std::vector<Position> Gmatch(
        const std::string& reference,
        const std::string& target,
        int kmer_length
    ) {
        constexpr int LIMIT = 100;
        constexpr int MAXSEQ = 1073741824;

        std::vector<int> kmer_location(MAXSEQ, -1);
        std::vector<int> next_kmer(reference.size(), -1);
        std::vector<Position> matches;
        std::set<int> used_reference_positions; // Track used reference indices

        // Build global hash table
        for (size_t i = 0; i + kmer_length <= reference.size(); ++i) {
            std::string kmer = reference.substr(i, kmer_length);
            long long key = std::abs(static_cast<long long>(std::hash<std::string>{}(kmer)));
            if (key == static_cast<long long>(INT_MIN)) key = 0;
            while (key > MAXSEQ - 1) key = key / 2;
            next_kmer[i] = kmer_location[(int)key];
            kmer_location[(int)key] = static_cast<int>(i);
        }

        int index = 0;
        int lastEIR = 0;
        while (index + kmer_length <= static_cast<int>(target.size())) {
            int increment = 0, most_incre = 0;
            std::string kmer = target.substr(index, kmer_length);
            long long key = std::abs(static_cast<long long>(std::hash<std::string>{}(kmer)));
            if (key == static_cast<long long>(INT_MIN)) key = 0;
            while (key > MAXSEQ - 1) key = key / 2;

            if (kmer_location[(int)key] == -1) {
                index = index + 1;
                continue;
            }

            int startIndex = INT_MAX;
            most_incre = 0;
            bool match = false;

            // Try to find matches in a limit range
            for (int k = kmer_location[(int)key]; k != -1; k = next_kmer[k]) {
                increment = 0;
                std::string Rkmer = reference.substr(k, kmer_length);
                if (kmer != Rkmer) continue;

                if (!matches.empty()) {
                    lastEIR = matches.back().end_reference;
                } else {
                    lastEIR = 0;
                }
                if (k - lastEIR > LIMIT || k - lastEIR < -LIMIT) continue;

                // Extend match
                int ref_idx = k + kmer_length;
                int tar_idx = index + kmer_length;
                while (ref_idx < (int)reference.size() && tar_idx < (int)target.size() &&
                       reference[ref_idx] == target[tar_idx]) {
                    ref_idx++;
                    tar_idx++;
                    increment++;
                }

                // Check for overlap in reference
                bool overlaps = false;
                for (int p = k; p <= k + kmer_length + increment - 1; ++p) {
                    if (used_reference_positions.count(p)) {
                        overlaps = true;
                        break;
                    }
                }
                if (overlaps) continue;

                match = true;
                if (increment == most_incre) {
                    if (matches.size() > 1) {
                        if (k - lastEIR < startIndex - lastEIR)
                            startIndex = k;
                    }
                } else if (increment > most_incre) {
                    most_incre = increment;
                    startIndex = k;
                }
            }

            // If no match in limit range, search whole sequence
            if (!match) {
                for (int k = kmer_location[(int)key]; k != -1; k = next_kmer[k]) {
                    increment = 0;
                    std::string Rkmer = reference.substr(k, kmer_length);
                    if (kmer != Rkmer) continue;

                    int ref_idx = k + kmer_length;
                    int tar_idx = index + kmer_length;
                    while (ref_idx < (int)reference.size() && tar_idx < (int)target.size() &&
                           reference[ref_idx] == target[tar_idx]) {
                        ref_idx++;
                        tar_idx++;
                        increment++;
                    }

                    // Check for overlap in reference
                    bool overlaps = false;
                    for (int p = k; p <= k + kmer_length + increment - 1; ++p) {
                        if (used_reference_positions.count(p)) {
                            overlaps = true;
                            break;
                        }
                    }
                    if (overlaps) continue;

                    if (increment == most_incre) {
                        if (matches.size() > 1) {
                            if (k - lastEIR < startIndex - lastEIR)
                                startIndex = k;
                        }
                    } else if (increment > most_incre) {
                        most_incre = increment;
                        startIndex = k;
                    }
                }
            }

            if (startIndex == INT_MAX) {
                index = index + 1;
                continue;
            }

            Position pos;
            pos.start_target = index;
            pos.end_target = index + kmer_length + most_incre - 1;
            pos.start_reference = startIndex;
            pos.end_reference = startIndex + kmer_length + most_incre - 1;
            matches.push_back(pos);

            // Mark used reference positions
            for (int p = pos.start_reference; p <= pos.end_reference; ++p) {
                used_reference_positions.insert(p);
            }

            index = pos.end_target + 1;
        }
        return matches;
    }
};

namespace FileUtils {

// Writes mismatches + (start,end) pairs with respect to the target sequence
void formatMatchesWithTarget(const std::vector<Encoder::Position>& positions, 
                             const std::string& target, 
                             const std::string& filename) {
    std::ostringstream oss;
    int prevEndTar = -1;

    for (size_t i = 0; i < positions.size(); ++i) {
        const auto& pos = positions[i];
        if (prevEndTar >= 0 && pos.start_target > prevEndTar + 1) {
            oss << target.substr(prevEndTar + 1, pos.start_target - prevEndTar - 1) << "\n";
        } else if (i == 0 && pos.start_target > 0) {
            oss << target.substr(0, pos.start_target) << "\n";
        }

        oss << pos.start_reference << "," << pos.end_reference << "\n";
        prevEndTar = pos.end_target;
    }

    if (!positions.empty() && prevEndTar < static_cast<int>(target.size()) - 1) {
        oss << target.substr(prevEndTar + 1) << "\n";
    }

    std::ofstream file(filename, std::ios::app);
    file << oss.str();
}

// Writes only (start,end) pairs, no target sequence output
void formatMatchesSimple(const std::vector<Encoder::Position>& positions, const std::string& filename) {
    std::ofstream file(filename, std::ios::app);
    for (const auto& pos : positions) {
        file << pos.start_reference << "," << pos.end_reference << "\n";
    }
}

// Write arbitrary text to a file, optionally appending
void writeTextToFile(const std::string& filename, const std::string& text, bool append) {
    std::ofstream file(filename, append ? std::ios::app : std::ios::trunc);
    file << text;
}

// Write positions delta-encoded with an auxiliary string at the start
void writePositionsDeltaEncoded(const std::string& filename, const std::vector<Encoder::Position>& positions, bool append, const std::string& auxiliary) {
    std::ofstream file(filename, append ? std::ios::app : std::ios::trunc);
    file << auxiliary;

    int prevEnd = 0;
    for (const auto& pos : positions) {
        int start = pos.start_target;
        int end = pos.end_target;
        file << (start - prevEnd) << " " << (end - start) << " ";
        prevEnd = end;
    }
    file << "\n";
}

// Postprocess input file by merging consecutive position ranges, then delta encode and write to output file
void postprocessPositions(const std::string& infile, const std::string& outfile) {
    std::ifstream inputFile(infile);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening input file: " << infile << "\n";
        return;
    }

    std::string line;
    std::stringstream merged;
    std::vector<int> range;

    // Merge consecutive ranges from the input file
    while (std::getline(inputFile, line)) {
        if (line.find(',') != std::string::npos) {
            auto commaPos = line.find(',');
            int begin = std::stoi(line.substr(0, commaPos));
            int end = std::stoi(line.substr(commaPos + 1));

            if (!range.empty() && range.back() != begin - 1) {
                merged << range.front() << "," << range.back() << "\n";
                range.clear();
            }
            range.push_back(begin);
            range.push_back(end);
        } else if (!line.empty()) {
            if (!range.empty()) {
                merged << range.front() << "," << range.back() << "\n";
                range.clear();
            }
            merged << line << "\n";
        }
    }
    if (!range.empty()) {
        merged << range.front() << "," << range.back() << "\n";
    }
    inputFile.close();

    // Delta encode merged ranges and write to output file
    std::istringstream mergedStream(merged.str());
    std::ostringstream deltaEncoded;
    int prev = 0;
    bool firstLine = true;

    while (std::getline(mergedStream, line)) {
        if (line.find(',') != std::string::npos) {
            int commaPos = line.find(',');
            int begin = std::stoi(line.substr(0, commaPos));
            int end = std::stoi(line.substr(commaPos + 1));

            if (firstLine) {
                deltaEncoded << begin << "," << (end - begin) << "\n";
                firstLine = false;
            } else {
                deltaEncoded << (begin - prev) << "," << (end - begin) << "\n";
            }
            prev = end;
        } else if (!line.empty()) {
            deltaEncoded << line << "\n";
        }
    }

    std::ofstream outputFile(outfile, std::ios::app);
    outputFile << deltaEncoded.str();
}

}  // namespace FileUtils

int main() {
    using namespace std;

    // Example genomes
    string reference = "ATCGGATTACGATCGTTAGGCTAGGCTAATCGGCTAGGATCGAAGCTTAGGCTTAGGCTAA";
    string target =    "TTAGGCTAGGCTAACCGGATCGTTAGGCTAGGCTAACGTTAGGCTAGGATCGTTAAGCTTAG";

    // Create Encoder and find matches with kmer_length = 5
    Encoder encoder;
    vector<Encoder::Position> matches = encoder.Gmatch(reference, target, 5);

    // Format matches with target (including unmatched regions)
    FileUtils::formatMatchesWithTarget(matches, target, "matches_with_target.txt");

    // Write simple matches (start, end pairs) to a file
    FileUtils::formatMatchesSimple(matches, "matches_simple.txt");

    // Delta encode matches and write to a file
    FileUtils::writePositionsDeltaEncoded("matches_delta.txt", matches, false, "Delta Encoded Matches:\n");

    // Postprocess matches (merge consecutive ranges and delta encode)
    FileUtils::postprocessPositions("matches_simple.txt", "matches_postprocessed.txt");

    // Output results to console
    cout << "Matches with target written to 'matches_with_target.txt'.\n";
    cout << "Simple matches written to 'matches_simple.txt'.\n";
    cout << "Delta encoded matches written to 'matches_delta.txt'.\n";
    cout << "Postprocessed matches written to 'matches_postprocessed.txt'.\n";

    return 0;
}
