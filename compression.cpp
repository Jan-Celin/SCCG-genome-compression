#include <iostream>
#include <unordered_map>
#include <vector>
#include <climits>
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

    std::vector<Position> Lmatch(
        const std::string& reference, 
        const std::string& target, 
        int kmer_length, 
        double T1, 
        double T2,
        int i = 0
    ) {
        create_hash_map_L(reference, kmer_length);
        std::vector<Position> pos_list;

        while (true) {
            std::string kmer = target.substr(i, kmer_length);
            int hash_value = std::hash<std::string>{}(kmer);

            if (kmer_map.find(hash_value) == kmer_map.end()) {
                i++;
                if (i > target.length() - kmer_length) {
                    break;
                }
            } else {
                int lmax1 = 0;
                int lmax2 = 0;
                int pn1 = 0;
                int pn2 = 0;
                int ln1 = 0;
                int ln2 = 0;
                int increment = 0;
                int max_increment = 0;
                int start_i = 0;

                std::vector<Kmer> kmer_list = kmer_map[hash_value];

                for (Kmer new_kmer : kmer_list) {
                    int kmer_start = new_kmer.kmerstart;
                    int end_reference = kmer_start + kmer_length - 1;
                    int end_target = i + kmer_length - 1;
                    
                    // Prosiri duljinu l pocevsi od zajednickog kmera, dok ne naides
                    // na znak koji nije jednak

                    if (new_kmer.kmer == kmer) {
                        kmer_start = new_kmer.kmerstart;
                        end_reference = kmer_start + kmer_length - 1;
                        end_target = i + kmer_length - 1;
                                int increment = 0;
                        int end_i;
                        if (reference.length() - 1 - end_reference < target.length() - 1 - end_target) {
                            end_i = reference.length() - end_reference;
                        } else {
                            end_i = target.length() - end_target;
                        }

                        for (int i=1; i<end_i; i++) {
                            if (reference[end_reference + i] == target[end_target + i]) {
                                increment++;
                            } else {
                                break;
                            }
                        }

                        if (kmer_list.size() == 1) {
                            max_increment = increment;
                            start_i = kmer_start;
                        } else {
                            if (increment > max_increment) {
                                max_increment = increment;
                                start_i = kmer_start;
                            } else if (increment == max_increment) {
                                if (kmer_start < start_i) {
                                    start_i = kmer_start;
                                }
                            }
                        }
                    }

                    if (start_i == INT_MAX) {
                        i++;
                        if (i > target.length() - kmer_length) {
                            break;
                        }
                        continue;
                    }

                    Position pos;
                    pos.start_reference = start_i;
                    pos.end_reference = start_i + kmer_length + max_increment - 1;
                    pos.start_target = i;
                    pos.end_target = i + kmer_length + max_increment - 1;
                    pos_list.push_back(pos);
                    if (i + kmer_length + max_increment - 1 > target.length() - kmer_length) {
                        break;
                    }
                }
            }

            return pos_list;
        }
    }
};

int main() {
    cout << "Hello World";

    return 0;
}