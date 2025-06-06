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

class Encoder {
public:
    // Hyperparameters (same as in Java)
    int kmer_length = 21; //21
    double T1 = 0.5;
    int T2 = 4;
    const int sub_length = 30000;
    bool local = true;

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

    void create_hash_map_L(const std::string& reference, int kmer_length) {
        std::string Nkmer = "";
        for (int j = 0; j < kmer_length; ++j) {
            Nkmer += 'N';
        }
        int Nkmer_hash_value = std::hash<std::string>{}(Nkmer);
        for (int i = 0; i < (int)reference.length() - kmer_length + 1; ++i) {
            std::string kmer = reference.substr(i, kmer_length);
            int hash_value = std::hash<std::string>{}(kmer);
            Kmer new_kmer;
            new_kmer.kmer = kmer;
            new_kmer.kmerstart = i;
            if (kmer_map.find(hash_value) == kmer_map.end()) {
                kmer_map[hash_value] = std::vector<Kmer>();
            }
            kmer_map[hash_value].push_back(new_kmer);

            // Skip N characters if hash_value corresponds to Nkmer.
            if (hash_value == Nkmer_hash_value) {
                while ((reference[i] == 'N' || reference[i] == 'n') && 
                       i < (int)reference.length() - kmer_length + 1) {
                    ++i;
                }
            }
        }
    }

    std::vector<Position> Lmatch(const std::string& reference, const std::string& target, 
                                 int kmer_length, double T1, double T2, int i = 0) {
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
                if (new_kmer.kmer != kmer)
                    continue;
                int kmer_start = new_kmer.kmerstart;
                int end_reference = kmer_start + kmer_length - 1;
                int end_target = i + kmer_length - 1;
                int max_extend = std::min((int)reference.size() - 1 - end_reference, 
                                          (int)target.size() - 1 - end_target);
                int increment = 0;
                for (; increment < max_extend; ++increment) {
                    if (reference[end_reference + 1 + increment] != target[end_target + 1 + increment])
                        break;
                }
                // Check if this region overlaps any already-used positions.
                bool overlaps = false;
                for (int p = kmer_start; p <= kmer_start + kmer_length + increment - 1; ++p) {
                    if (used_reference_positions.count(p)) {
                        overlaps = true;
                        break;
                    }
                }
                if (overlaps)
                    continue;
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
            // Mark used positions.
            for (int p = pos.start_reference; p <= pos.end_reference; ++p) {
                used_reference_positions.insert(p);
            }
            i = pos.end_target + 1;  // Advance past the matched region.
        }
        return pos_list;
    }

    std::vector<Position> Gmatch(const std::string& reference, const std::string& target, int kmer_length) {
        constexpr int LIMIT = 100;  // same as Java's limit
        constexpr int MAXSEQ = 1073741824;  // same as Java's maxseq
        std::vector<int> kmer_location(MAXSEQ, -1);
        std::vector<int> next_kmer(reference.size(), -1);
        std::vector<Position> matches;
        std::set<int> used_reference_positions; // Track used reference indices

        // Build global hash table.
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
        while (index + kmer_length <= (int)target.size()) {
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
            // Try to find matches in a limit range.
            for (int k = kmer_location[(int)key]; k != -1; k = next_kmer[k]) {
                increment = 0;
                std::string Rkmer = reference.substr(k, kmer_length);
                if (kmer != Rkmer) continue;
                if (!matches.empty())
                    lastEIR = matches.back().end_reference;
                else
                    lastEIR = 0;
                if (k - lastEIR > LIMIT || k - lastEIR < -LIMIT)
                    continue;
                // Extend match.
                int ref_idx = k + kmer_length;
                int tar_idx = index + kmer_length;
                while (ref_idx < (int)reference.size() && tar_idx < (int)target.size() &&
                       reference[ref_idx] == target[tar_idx]) {
                    ref_idx++;
                    tar_idx++;
                    increment++;
                }
                // Check for overlap.
                bool overlaps = false;
                for (int p = k; p <= k + kmer_length + increment - 1; ++p) {
                    if (used_reference_positions.count(p)) {
                        overlaps = true;
                        break;
                    }
                }
                if (overlaps)
                    continue;
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
            // If no match in limit, search whole sequence.
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
                    bool overlaps = false;
                    for (int p = k; p <= k + kmer_length + increment - 1; ++p) {
                        if (used_reference_positions.count(p)) {
                            overlaps = true;
                            break;
                        }
                    }
                    if (overlaps)
                        continue;
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
            for (int p = pos.start_reference; p <= pos.end_reference; ++p) {
                used_reference_positions.insert(p);
            }
            index = pos.end_target + 1;
        }
        return matches;
    }

    std::string readFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open file: " + filename);
        }
        std::ostringstream oss;
        oss << file.rdbuf();
        return oss.str();
    }

    void writeTextToFile(const std::string& filename, const std::string& text, bool append) {
        std::ofstream file(filename, append ? std::ios::app : std::ios::trunc);
        if (!file.is_open()) {
            throw std::runtime_error("Could not write to file: " + filename);
        }
        file << text;
    }

    void postprocessPositions(const std::string& infile, const std::string& outfile) {
        std::ifstream inputFile(infile);
        if (!inputFile.is_open()) {
            std::cerr << "Error opening input file: " << infile << "\n";
            return;
        }
        std::string line;
        std::stringstream merged;
        std::vector<int> range;
        // Merge consecutive ranges.
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
        std::ofstream outputFile(outfile, std::ios::trunc);
        outputFile << deltaEncoded.str();
    }

    void run7zip(const std::string& filename, const std::string& outputFolder) {
        std::string compressed = outputFolder + "/compressed.7z";
        std::string command = "\"C:\\7-Zip\\7z.exe\" a -t7z \"" + compressed + "\" \"" + filename + "\"";
        if (std::system(command.c_str()) != 0) {
            std::cerr << "7zip compression failed for " << filename << "\n";
        } else {
            std::cout << "File compressed successfully: " << compressed << "\n";
        }
    }

    void compress_genome(const std::string& referenceFile, const std::string& targetFile,
                         const std::string& outputFolder) {
        if (!std::filesystem::exists(outputFolder)) {
            std::filesystem::create_directory(outputFolder);
        }
        std::string refSeq = readFile(referenceFile);
        std::string tarSeq = readFile(targetFile);

        auto posRef = refSeq.find('\n');
        if (posRef != std::string::npos) {
            refSeq = refSeq.substr(posRef + 1);
        }
        auto posTar = tarSeq.find('\n');
        std::string meta_data;
        if (posTar != std::string::npos) {
            meta_data = tarSeq.substr(0, posTar);
            tarSeq = tarSeq.substr(posTar + 1);
        }

        if (tarSeq.size() < sub_length * 5) {
            local = false;
            std::cout << "Target too short for local matching. Will use global matching.\n";
        } else if (tarSeq.size() < sub_length * 1333) {
            T1 = 0.1;
            T2 = 0;
            std::cout << "Adjusted thresholds for short target.\n";
        }

        std::string finalFile = outputFolder + "/final.txt";
        std::string interimFile = outputFolder + "/interim.txt";
        {
            std::ofstream temp(interimFile, std::ios::trunc);
            if (!temp.is_open()) {
                throw std::runtime_error("Could not create interim file: " + interimFile);
            }
        }
        writeTextToFile(finalFile, meta_data + "\n" + std::to_string(tarSeq.size()) + "\n", false);

        bool fallbackGlobal = false;
        std::string accumulatedText;
        if (local) {
            // Local Matching Phase.
            int sot = 0;                     // Start offset in target.
            int eot = sub_length;            // End offset for target segment.
            int sor = 0;                     // Start offset in reference.
            int eor = sub_length;            // End offset for reference segment.
            int mismatch = 0;
            while (sot < (int)tarSeq.size()) {
                // If the segment exceeds target or reference bounds, append the rest.
                if (eot >= (int)tarSeq.size() || eor >= (int)refSeq.size()) {
                    accumulatedText += tarSeq.substr(sot);
                    std::cout << "Reached end of sequence, appending remainder.\n";
                    break;
                }
                std::string refSegment = refSeq.substr(sor, eor - sor);
                std::string tarSegment = tarSeq.substr(sot, eot - sot);
                auto matches = Lmatch(refSegment, tarSegment, kmer_length, T1, T2);
                if (matches.empty()) {
                    // Try with a lower kmer length.
                    int temp_kmer = 11;
                    matches = Lmatch(refSegment, tarSegment, temp_kmer, T1, T2);
                }
                if (matches.empty()) {
                    // No match found; count a mismatch and output the raw segment.
                    mismatch++;
                    accumulatedText += tarSegment + "\n";
                    std::cout << "No local match found at target pos " << sot 
                              << ", mismatch count = " << mismatch << "\n";
                    sot += sub_length;
                    eot = sot + sub_length;
                    sor += sub_length;
                    eor = sor + sub_length;
                    // Fallback if mismatches exceed threshold.
                    if (mismatch > T2) {
                        fallbackGlobal = true;
                        std::cout << "Too many mismatches (" << mismatch 
                                  << "), falling back to global matching.\n";
                        break;
                    }
                    continue;
                }
                // Merge found positions.
                Position mergedPos = format_matches(matches);
                accumulatedText += std::to_string(mergedPos.start_reference) + "," +
                                   std::to_string(mergedPos.end_reference) + "\n";
                std::cout << "Local match: target [" << sot << "," << mergedPos.end_target 
                          << "] mapped to reference [" << mergedPos.start_reference 
                          << "," << mergedPos.end_reference << "]\n";
                sot += (mergedPos.end_target + 1);
                eot = sot + sub_length;
                sor += (mergedPos.end_reference + 1);
                eor = sor + sub_length;
            }
            writeTextToFile(interimFile, accumulatedText, true);
        }
        
        // Global Matching Phase fallback.
        if (fallbackGlobal || accumulatedText.empty()) {
            std::cout << "Using global matching.\n";
            // Clear or recreate the interim file.
            {
                std::ofstream temp(interimFile, std::ios::trunc);
                if (!temp.is_open()) {
                    throw std::runtime_error("Could not recreate interim file: " + interimFile);
                }
            }
            std::string globalText;
            auto globalMatches = Gmatch(refSeq, tarSeq, kmer_length);
            if (globalMatches.empty()) {
                std::cout << "Global matching did not find any matches. Output raw target.\n";
                globalText = tarSeq;
            } else {
                for (const auto &pos : globalMatches) {
                    globalText += std::to_string(pos.start_reference) + "," +
                                  std::to_string(pos.end_reference) + "\n";
                    std::cout << "Global match: reference [" << pos.start_reference << "," 
                              << pos.end_reference << "]\n";
                }
            }
            writeTextToFile(interimFile, globalText, true);
        }

        // Postprocess interim file into final file.
        postprocessPositions(interimFile, finalFile);
        std::cout << "Compression complete. Final file is located at " << finalFile << "\n";
        // Clean up temporary file.
        std::remove(interimFile.c_str());
    }

    // Merges a vector of Position objects.
    Position format_matches(const std::vector<Position>& matches) {
        Position merged;
        if (matches.empty()) return merged;
        merged.start_reference = matches[0].start_reference;
        merged.end_reference = matches[0].end_reference;
        merged.start_target = matches[0].start_target;
        merged.end_target = matches[0].end_target;
        for (size_t i = 1; i < matches.size(); ++i) {
            merged.start_reference = std::min(merged.start_reference, matches[i].start_reference);
            merged.end_reference = std::max(merged.end_reference, matches[i].end_reference);
            merged.start_target = std::min(merged.start_target, matches[i].start_target);
            merged.end_target = std::max(merged.end_target, matches[i].end_target);
        }
        return merged;
    }
};
  
// Main function.
int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <reference_file> <target_file> <output_folder>\n";
        return 1;
    }
    try {
        // Ensure output folder exists.
        if (!std::filesystem::exists(argv[3])) {
            std::filesystem::create_directory(argv[3]);
        }
        Encoder encoder;
        encoder.compress_genome(argv[1], argv[2], argv[3]);
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
