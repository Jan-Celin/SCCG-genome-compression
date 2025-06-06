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
#include <fstream>
#include <stdexcept>
#include <ranges>

using namespace std;

class Encoder {
public:
    std::string reference, target;
    int kmer_length = 21;
    double T1 = 0.5;
    int T2 = 4;
    int maxchar = 268435456;
    int maxseq = 1073741824;
    int limit = 100;
    SUB_LENGTH = 30000;

    std::string meta_data = "";
    std::string text = "";

    int length = 0;
    int mismatch = 0;
    int endref = SUB_LENGTH - 1;
    int sot = 0;
    int eot = 0;
    int sor = 0;
    int eor = 0;

    class Kmer {
    public:
        std::string kmer;
        int kmerstart;
    };
    
    std::unordered_map<int, std::vector<Kmer>> kmer_map;
    std::vector<int> next_kmer;
    std::vector<int> kmer_location;

    bool local = true;

    class Position {
    public:
        int start_reference;
        int end_reference;
        int start_target;
        int end_target;

        Position() : start_reference(0), end_reference(0), start_target(0), end_target(0) {}
    };

    void create_hash_map_L(
        const std::string& reference, 
        int kmer_length
    ) {
        kmer_map.clear();
        std::string Nkmer = "";
        for (int j=0; j<kmer_length; ++j) {
            Nkmer += 'N';
        }

        int Nkmer_hash_value = std::hash<std::string>{}(Nkmer);

        for (int i=0; i < reference.length() - kmer_length + 1;) {
            std::string kmer = reference.substr(i, kmer_length);
            int hash_value = std::hash<std::string>{}(kmer);

            Kmer new_kmer;
            new_kmer.kmer = kmer;
            new_kmer.kmerstart = i;

            
            if (kmer_map.find(hash_value) == kmer_map.end()) {
                kmer_map[hash_value] = std::vector<Kmer>();
            }
            kmer_map[hash_value].push_back(new_kmer);

            ++i;

            // Preskoci n znakove
            if (hash_value == Nkmer_hash_value) {
                while ((reference[i] == 'N' || reference[i] == 'n') && i < reference.length() - kmer_length + 1) {
                    ++i;
                }
            }
        }
    }

    std::vector<Position> Lmatch(
        const std::string& reference, 
        const std::string& target, 
        int kmer_length
    ) {
        create_hash_map_L(reference, kmer_length);
        std::vector<Position> pos_list;
        int i = 0;

        while (true) {
            int increment = 0;
            int max_extend = 0;
            int start_index = INT_MAX;

            if (i + kmer_length > target.size()) {
                break;  // Vise nema dovoljno znakova za izraditi kmer
            }

            std::string kmer = target.substr(i, kmer_length);
            int kmer_hash = std::hash<std::string>{}(kmer);

            if (kmer_map.find(kmer_hash) == kmer_map.end()) {
                ++i;  // Kmer nije pronaden, preskoci jedan znak
                continue;
            }

            const std::vector<Kmer>& kmer_list = kmer_map[kmer_hash];

            for (const auto& new_kmer : kmer_list) {
                if (new_kmer.kmer != kmer) continue;

                int kmer_start = new_kmer.kmerstart;
                int end_reference = kmer_start + kmer_length - 1;
                int end_target = i + kmer_length - 1;

                // Pronadi najveci broj znakova koji se podudaraju
                int max_extend = std::min((int)reference.size() - 1 - end_reference, (int)target.size() - 1 - end_target);
                for (increment; increment < max_extend; ++increment) {
                    if (reference[end_reference + 1 + increment] != target[end_target + 1 + increment]) {
                        break;
                    }
                }

                if (kmer_list.size() > 1) {
                    if (increment == max_extend) {
                        // Ako su svi znakovi podudarni, uzmi prvi kmer
                        if (pos_list.size() > 1) {
                            int last_end_reference = pos_list.back().end_reference;
                            if (kmer_start < start_index) {
                                start_index = kmer_start;
                            }
                        }
                    }  
                    else if (increment > max_extend) {
                        // Ako je pronaden veci broj podudarnih znakova, azuriraj
                        start_index = kmer_start;
                        max_extend = increment;
                    }
                } else {
                    max_extend = increment;
                    start_index = kmer_start;
                    break;
                }
            }

            if (start_index == INT_MAX) {
                ++i;
                continue;
            }

            Position pos;
            pos.start_reference = start_index;
            pos.end_reference = start_index + kmer_length + max_extend - 1;
            pos.start_target = i;
            pos.end_target = i + kmer_length + max_extend - 1;
            pos_list.push_back(pos);

            i += kmer_length + max_extend + 1;
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

    std::string read_sequence_l_match(const std::string& file_path) {
        std::ifstream file(file_path);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open file: " + file_path);
        }

        std::string sequence;
        std::string line;
        std::string metadata;

        std::getline(file, metadata);
        while (std::getline(file, line)) {
            sequence += line;
        }

        file.close();

        return sequence;
    }

    // Writes mismatches + (start,end) pairs with respect to the target sequence
    void formatMatchesWithTarget(const std::vector<Position>& positions, 
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
    void formatMatchesSimple(const std::vector<Position>& positions, const std::string& filename) {
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
    void writePositionsDeltaEncoded(const std::string& filename, const std::vector<Position>& positions, bool append, const std::string& auxiliary) {
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

    std::vector<Position> lowercase_position(const std::string& sequence) {
        std::vector<Position> positions;
        
        bool successive = false;
        int start = 0;
        int end = 0;

        for (int i=0; i<sequence.length(); ++i) {
            if (std::islower(sequence[i])) {
                if (successive) {
                    end += 1;
                } else {
                    successive = true;
                    start = i;
                    end += 1;
                }
            } else {
                if (successive) {
                    Position pos;
                    pos.start_target = start;
                    pos.end_target = end - 1;
                    positions.push_back(pos);
                }
                successive = false;
                start = 0;
                end = i + 1;
            }
        }

        if (successive) {
            Position pos;
            pos.start_target = start;
            pos.end_target = end;
            positions.push_back(pos);
        }

        return positions;
    }

    Position format_matches(List<Position> matches) {
        int start_reference;
        int end_reference;
        int start_target;
        int trouble = 0;

        for (int i = 0; i < list.size(); ++i) {
            if (i == 0) {
                start_reference = matches[i].start_reference;
                end_reference = matches[i].end_reference;
                start_target = matches[i].start_target;
                if (end_reference >= 0)
            }
        }
    }

    void save_7zip(const std::string& file_path, const std::string& output_dir) {
        std::string command = "7z a -t7z \"" + output_dir + "/compressed.7z\" \"" + file_path + "\"";
        int result = system(command.c_str());
        if (result != 0) {
            std::cerr << "Error compressing file with 7zip: " << file_path << "\n";
        } else {
            std::cout << "File compressed successfully: " << file_path << "\n";
        }
    }

    void compress_genome(
        const std::string& reference_genome_path,
        const std::string& target_genome_path,
        const std::string& output_dir,
        int kmer_length = 21
    ) {
        int controuble = 1;
        bool is_con = false;

        std::cout << "Sazimanje genoma: " << target_genome_path << std::endl;
        int mismatch = 0;

        std::remove(output_file.c_str());

        std::string chr_file_name = std::filesystem::path(target_genome_path).filename().string();
        std::string output_file = (std::filesystem::path(output_dir) / chr_file_name).string();
        std::string temp_file = (std::filesystem::path(output_dir) / (chr_file_name + ".temp")).string();

        std::remove(output_file.c_str());

        // TODO: Izmjeriti vremensko i prostorno zauzece algoritma.

        while (local) {
            mismatch = 0;
            std::string reference_sequence = read_sequence_l_match(reference_genome_path);
            std::string target_sequence = read_sequence_l_match(target_genome_path);

            int target_length = target_sequence.length();

            if (target_sequence.length() < SUB_LENGTH*5) {
                local = false;
                break;
            } else if (target_sequence.length() < SUB_LENGTH*1333) {
                T1 = 0.1;
                T2 = 0;
            }

            std::vector<Position> L_list =lowercase_position(target_sequence);

            // std::string meta_data = "neki_metapodaci\n" + std::to_string(target_length) + "\n";

            writePositionsDeltaEncoded(output_file, L_list, false, "Lmatch positions for " + chr_file_name + "\n");
            writeTextToFile(output_file, "\n", true);

            std::transform(reference_sequence.begin(), reference_sequence.end(), reference_sequence.begin(), ::toupper);
            std::transform(target_sequence.begin(), target_sequence.end(), target_sequence.begin(), ::toupper);

            std::remove(temp_file.c_str());

            int sot = 0, eot = SUB_LENGTH, sor = 0, eor = SUB_LENGTH;

            Encoder encoder;
            Position pos;

            while (true) {
                if (eor > reference_sequence.length() || eot > target_sequence.length()) {
                    std::string text = target_sequence.substr(sot);
                    if (text.length() <= 0) {
                        break;
                    } else {
                        FileUtils::writeTextToFile(temp_file, text, true);
                        break;
                    }
                }

                std::string reference = reference_sequence.substr(sor, eor - sor);
                std::string target = target_sequence.substr(sot, eot - sot);

                std::vector<Position> matches = encoder.Lmatch(reference, target, kmer_length);

                if (matches.size() <= 0) {
                    kmer_length = 11;
                    matches = encoder.Lmatch(reference, target, kmer_length);
                }

                if (matches.size() <= 0) {
                    mismatch++;

                    if (eot >= target_sequence.length() - 1) {
                        std::string text = target_sequence.substr(sot);
                        writeTextToFile(temp_file, text, true);
                        break;
                    }

                    if (is_con) {
                        controuble++;
                    }
                    is_con = true;

                    writeTextToFile(temp_file, target + "\n", true);
                    sot += SUB_LENGTH;
                    eot = sot + SUB_LENGTH;
                    sor += SUB_LENGTH;
                    eor = sor + SUB_LENGTH;

                    int difference = target_sequence.length() - sot;
                    if (difference <= kmer_length) {
                        text = target_sequence.substr(sot);
                        writeTextToFile(temp_file, text, true);
                        break;
                    } else if (difference < SUB_LENGTH) {
                        eot = target_sequence.length() - 1;
                    }

                    int difference_ref = reference_sequence.length() - sor;

                    if (difference_ref < SUB_LENGTH) {
                        eor = reference_sequence.length() - 1;
                    }
                    if (eot >= target_sequence.length()) {
                        break;
                    }
                    if (difference_ref <= kmer_length) {
                        text = target_sequence.substr(sot);
                        if (text.length() <= 0) {
                            break;
                        } else {
                            if (text.length() >= SUB_LENGTH * T1) {
                                mismatch++;
                            }
                            if (mismatch > T2) {
                                local = false;
                                break;
                            }
                            writeTextToFile(temp_file, text, true);
                            break;
                        }
                    }
                    continue;
                }

                is_con = false;
                if (controuble > 2) {
                    mismatch -= controuble;
                }
                controuble = 1;

                position = format_matches(matches)

                if (mismatch > T2) {
                    local = false;
                    break;
                }

                sot += position.end_target + 1;
                eot = sot + SUB_LENGTH;
                sor += position.end_reference + 1;
                eor = sor + SUB_LENGTH;

                writeTextToFile(temp_file, text, true);

                int difference = target_sequence.length() - sot;

                if (difference <= kmer_length) {
                    text = target_sequence.substr(sot);
                    if (text.length() <= 0) {
                        break;
                    } else {
                        writeTextToFile(temp_file, text, true);
                        break;
                    }
                } else if (difference < SUB_LENGTH) {
                    eot = target_sequence.length() - 1;
                }

                int difference_ref = reference_sequence.length() - sor;

                if (difference_ref < SUB_LENGTH) {
                    eor = reference_sequence.length() - 1;
                }
                if (difference_ref <= kmer_length) {
                    text = target_sequence.substr(sot);
                    if (text.length() <= 0) {
                        break;
                    } else {
                        if (text.length() > SUB_LENGTH * T1) {
                            mismatch++;
                        }
                        if (mismatch > T2) {
                            local = false;
                            break;
                        }
                        writeTextToFile(temp_file, text, true);
                        break;
                    }
                }
            }

            if (!local) {
                break;
            }

            postprocessPositions(temp_file, output_file);

            // TODO: use7zip
            break;
        }
        
        if (!local) {
            st::string reference_sequence = 
        }
    }
};


int main(int argc, char* argv[]) {
    Encoder encoder;
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <reference_genome_file> <target_genome_file> <output_dir>\n";
        return 1;
    }

    std::string reference_genome_path = argv[1];
    std::string target_genome_path = argv[2];
    std::string output_path = argv[3];

    encoder.compress_genome(reference_genome_path, target_genome_path, output_dir);
    
    return 0;
}
