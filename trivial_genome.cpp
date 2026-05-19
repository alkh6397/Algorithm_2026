#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <unordered_map>
using namespace std;

// Trivial Mapping: slide read along reference, count mismatches
// if mismatches <= threshold, it is a match
int trivialMap(const string& ref, const string& read, int threshold) {
    int n = ref.size(), m = read.size();
    int bestPos = -1;
    int bestMismatch = threshold + 1;

    for (int i = 0; i <= n - m; i++) {
        int mismatch = 0;
        for (int j = 0; j < m; j++) {
            if (ref[i+j] != read[j]) mismatch++;
            if (mismatch > threshold) break; // early exit
        }
        if (mismatch <= threshold && mismatch < bestMismatch) {
            bestMismatch = mismatch;
            bestPos = i;
        }
    }
    return bestPos; // -1 if no match found
}

// read a single line from file
string readSingleLine(const string& path) {
    ifstream f(path);
    if (!f) { cerr << "Failed to open file: " << path << endl; exit(1); }
    string line;
    getline(f, line);
    return line;
}

// read all reads from file (one read per line)
vector<string> readReads(const string& path) {
    ifstream f(path);
    if (!f) { cerr << "Failed to open file: " << path << endl; exit(1); }
    vector<string> reads;
    string line;
    while (getline(f, line))
        if (!line.empty()) reads.push_back(line);
    return reads;
}

// snp_log.txt format: "index original_base variant_base"
struct SNP {
    int index;
    char original;
    char variant;
};

vector<SNP> readSNPLog(const string& path) {
    ifstream f(path);
    if (!f) { cerr << "Failed to open file: " << path << endl; exit(1); }
    vector<SNP> snps;
    int idx; char orig, var;
    while (f >> idx >> orig >> var)
        snps.push_back({idx, orig, var});
    return snps;
}

// map each read to reference, overwrite result at best matched position
// duplicate reads are cached to avoid recomputation
string reconstruct(const string& ref, const vector<string>& reads, int threshold,
                   int& mapped, int& unmapped) {
    int n = ref.size();
    string result(n, 'N');
    unordered_map<string, int> posCache;
    mapped = 0; unmapped = 0;

    int total = reads.size();
    int interval = max(1, total / 10);

    for (int idx = 0; idx < total; idx++) {
        const string& read = reads[idx];

        auto [it, inserted] = posCache.emplace(read, -1);
        if (inserted)
            it->second = trivialMap(ref, read, threshold);

        int bestPos = it->second;
        int L = read.size();
        if (bestPos >= 0 && bestPos + L <= n) {
            mapped++;
            for (int i = 0; i < L; i++)
                result[bestPos+i] = read[i];
        } else {
            unmapped++;
        }

        if ((idx + 1) % interval == 0 || idx + 1 == total)
            cout << "Progress: " << (idx + 1) << " / " << total
                 << " (" << (idx + 1) * 100 / total << "%)" << endl;
    }
    return result;
}

// overall accuracy: reconstructed vs variant (ground truth)
// 'N' (not covered) is counted as wrong
double calcAccuracy(const string& variant, const string& reconstructed) {
    int n = variant.size(), match = 0;
    for (int i = 0; i < n; i++)
        if (reconstructed[i] == variant[i]) match++;
    return (double)match / n * 100.0;
}

// SNP accuracy: check only SNP positions
void calcSNPAccuracy(const vector<SNP>& snps, const string& reconstructed) {
    int total = snps.size(), correct = 0, uncovered = 0;
    for (const SNP& snp : snps) {
        if (snp.index >= (int)reconstructed.size()) continue;
        char restored = reconstructed[snp.index];
        if (restored == 'N') {
            uncovered++;
        } else if (restored == snp.variant) {
            correct++;
        }
    }
    cout << "SNP accuracy: " << correct << " / " << total
         << " (" << (total ? (double)correct/total*100.0 : 0.0) << "%)" << endl;
    if (uncovered > 0)
        cout << "  Uncovered SNP positions: " << uncovered << endl;
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        cout << "Usage: " << argv[0]
             << " <reference.txt> <reads_shuffled.txt> <variant.txt> <snp_log.txt>" << endl;
        return 1;
    }

    string reference     = readSingleLine(argv[1]);
    vector<string> reads = readReads(argv[2]);
    string variant       = readSingleLine(argv[3]);
    vector<SNP> snps     = readSNPLog(argv[4]);
    const int threshold  = 2;

    cout << "=== Trivial Mapping Algorithm ===" << endl;
    cout << "Reference length: " << reference.size() << endl;
    cout << "Variant length:   " << variant.size() << endl;
    cout << "Number of reads:  " << reads.size() << endl;
    cout << "Read length:      " << (reads.empty() ? 0 : (int)reads[0].size()) << endl;
    cout << "Number of SNPs:   " << snps.size() << endl;
    cout << "Threshold:        " << threshold << endl << endl;

    int mapped = 0, unmapped = 0;

    auto start = chrono::high_resolution_clock::now();
    string reconstructed = reconstruct(reference, reads, threshold, mapped, unmapped);
    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double, milli>(end - start).count();

    int covered = 0;
    for (char c : reconstructed) if (c != 'N') covered++;

    double accuracy = calcAccuracy(variant, reconstructed);

    cout << "Mapped reads:      " << mapped << endl;
    cout << "Unmapped reads:    " << unmapped << endl;
    cout << "Covered positions: " << covered << " / " << reference.size() << endl;
    cout << "Overall accuracy:  " << accuracy << "%" << endl;
    calcSNPAccuracy(snps, reconstructed);
    cout << "Elapsed time:      " << elapsed << "ms" << endl;

    ofstream out("restored_genome_trivial.txt");
    if (!out) { cerr << "Failed to save restored_genome_trivial.txt" << endl; return 1; }
    out << reconstructed;
    out.close();
    cout << "Restored genome saved to restored_genome_trivial.txt" << endl;

    return 0;
}
