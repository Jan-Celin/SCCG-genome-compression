#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <climits>
#include <fstream>
#include <sstream>
#include <functional>
#include <cmath>

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

        while (i <= (int)target.length() - kmer_length) {
            std::string kmer = target.substr(i, kmer_length);
            int hash_value = std::hash<std::string>{}(kmer);

            if (kmer_map.find(hash_value) == kmer_map.end()) {
                // No match, move forward
                i++;
                continue;
            }

            const std::vector<Kmer>& kmer_list = kmer_map[hash_value];
            int max_increment = -1;
            int best_start_ref = -1;

            // Find the best extension of the kmer match
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

                if (increment > max_increment || (increment == max_increment && kmer_start < best_start_ref)) {
                    max_increment = increment;
                    best_start_ref = kmer_start;
                }
            }

            if (best_start_ref == -1) {
                // No valid match found, move forward
                i++;
                continue;
            }

            Position pos;
            pos.start_reference = best_start_ref;
            pos.end_reference = best_start_ref + kmer_length + max_increment - 1;
            pos.start_target = i;
            pos.end_target = i + kmer_length + max_increment - 1;
            pos_list.push_back(pos);

            // Advance i to the end of the matched region + 1 to avoid overlapping matches
            i = pos.end_target + 1;
        }

        return pos_list;
    }


    std::vector<Position> Gmatch(
        const std::string& reference,
        const std::string& target,
        int kmerLength
    ) {
        constexpr int LIMIT = 100;
        constexpr int MAXSEQ = 1073741824;

        std::vector<int> kmerLocation(MAXSEQ, -1);
        std::vector<int> nextKmer(reference.size(), -1);
        std::vector<Position> matches;

        for (size_t i = 0; i + kmerLength <= reference.size(); ++i) {
            std::string kmer = reference.substr(i, kmerLength);
            size_t key = std::hash<std::string>{}(kmer) % MAXSEQ;
            nextKmer[i] = kmerLocation[key];
            kmerLocation[key] = static_cast<int>(i);
        }

        int index = 0;
        while (index + kmerLength <= static_cast<int>(target.size())) {
            std::string kmer = target.substr(index, kmerLength);
            size_t key = std::hash<std::string>{}(kmer) % MAXSEQ;

            if (kmerLocation[key] == -1) {
                ++index;
                continue;
            }

            int bestStart = INT_MAX;
            int maxIncrement = 0;
            int lastEndInRef = matches.empty() ? 0 : matches.back().end_reference;

            for (int pos = kmerLocation[key]; pos != -1; pos = nextKmer[pos]) {
                if (reference.compare(pos, kmerLength, kmer) != 0) continue;
                if (std::abs(pos - lastEndInRef) > LIMIT) continue;

                int refIdx = pos + kmerLength;
                int tarIdx = index + kmerLength;
                int increment = 0;

                while (refIdx < (int)reference.size() && tarIdx < (int)target.size() &&
                    reference[refIdx] == target[tarIdx]) {
                    ++refIdx; ++tarIdx; ++increment;
                }

                if (increment > maxIncrement ||
                    (increment == maxIncrement && std::abs(pos - lastEndInRef) < std::abs(bestStart - lastEndInRef))) {
                    bestStart = pos;
                    maxIncrement = increment;
                }
            }

            if (bestStart == INT_MAX) {
                ++index;
                continue;
            }

            matches.push_back(Position{
                bestStart,
                bestStart + kmerLength + maxIncrement - 1,
                index,
                index + kmerLength + maxIncrement - 1
            });

            index += kmerLength + maxIncrement + 1;
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
    string reference = "ATCGATCGATCG";
    string target =    "ATCGTTTCGATCG";

    // Create Encoder and find matches with kmer_length = 3, T1 and T2 unused in Lmatch here
    Encoder encoder;
    vector<Encoder::Position> matches = encoder.Lmatch(reference, target, 3, 0, 0);

    // In-memory version of formatMatchesWithTarget (returns string)
    auto formatMatchesWithTargetToString = [&](const vector<Encoder::Position>& positions, const string& target) {
        ostringstream oss;
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
        return oss.str();
    };

    // In-memory version of formatMatchesSimple (returns string)
    auto formatMatchesSimpleToString = [&](const vector<Encoder::Position>& positions) {
        ostringstream oss;
        for (const auto& pos : positions) {
            oss << pos.start_reference << "," << pos.end_reference << "\n";
        }
        return oss.str();
    };

    cout << "Matches with target (including unmatched target sequences):\n";
    cout << formatMatchesWithTargetToString(matches, target) << "\n";

    cout << "Simple matches (only position pairs):\n";
    cout << formatMatchesSimpleToString(matches) << "\n";

    return 0;
}
