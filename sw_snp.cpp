#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>
#include <unordered_map>
using namespace std;

// Smith-Waterman: full DP table
// MATCH=2, MISMATCH=-1, GAP=-1
// returns best score, sets bestPos to estimated alignment start (0-indexed)
int SmithWaterman(const string& ref, const string& pattern, int& bestPos) {
    int n = ref.size(), m = pattern.size();
    const int MATCH = 2, MISMATCH = -1, GAP = -1;

    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));

    int bestScore = 0;
    int bestCol = 0, bestRow = 0;

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int diag = dp[i-1][j-1] + (pattern[i-1] == ref[j-1] ? MATCH : MISMATCH);
            int up   = dp[i-1][j] + GAP;
            int left = dp[i][j-1] + GAP;
            dp[i][j] = max({0, diag, up, left});

            if (dp[i][j] > bestScore) {
                bestScore = dp[i][j];
                bestCol   = j;
                bestRow   = i;
            }
        }
    }

    bestPos = bestCol - bestRow;
    if (bestPos < 0) bestPos = 0;
    if (bestPos + m > n) bestPos = n - m;

    return bestScore;
}

// read a single line from file (for reference and variant)
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

// map each read to reference using SW, overwrite result at mapped position
string reconstruct(const string& ref, const vector<string>& reads) {
    int n = ref.size();
    string result(n, 'N');
    unordered_map<string, int> posCache;

    int total = reads.size();
    int interval = max(1, total / 10);

    for (int idx = 0; idx < total; idx++) {
        const string& read = reads[idx];

        auto [it, inserted] = posCache.emplace(read, 0);
        if (inserted)
            SmithWaterman(ref, read, it->second);

        int bestPos = it->second;
        int L = read.size();
        if (bestPos >= 0 && bestPos + L <= n)
            for (int i = 0; i < L; i++)
                result[bestPos+i] = read[i];

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
    if (argc < 6) {
        cout << "Usage: " << argv[0]
             << " <reference.txt> <reads_shuffled.txt> <variant.txt> <snp_log.txt> <snp_rate>" << endl;
        return 1;
    }

    string reference     = readSingleLine(argv[1]);
    vector<string> reads = readReads(argv[2]);
    string variant       = readSingleLine(argv[3]);
    vector<SNP> snps     = readSNPLog(argv[4]);
    double snpRate       = stod(argv[5]);

    int readLen = reads.empty() ? 0 : (int)reads[0].size();

    cout << "=== Smith-Waterman (SNP Rate Experiment) ===" << endl;
    cout << "Reference length: " << reference.size() << endl;
    cout << "Variant length:   " << variant.size() << endl;
    cout << "Number of reads:  " << reads.size() << endl;
    cout << "Read length:      " << readLen << endl;
    cout << "Number of SNPs:   " << snps.size() << endl;
    cout << "SNP rate:         " << snpRate * 100 << "%" << endl << endl;

    auto start = chrono::high_resolution_clock::now();
    string reconstructed = reconstruct(reference, reads);
    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double, milli>(end - start).count();

    int covered = 0;
    for (char c : reconstructed) if (c != 'N') covered++;

    double accuracy = calcAccuracy(variant, reconstructed);
    int snpCorrect = 0, snpTotal = snps.size();
    for (const SNP& snp : snps)
        if (snp.index < (int)reconstructed.size() && reconstructed[snp.index] == snp.variant)
            snpCorrect++;

    cout << "Covered positions: " << covered << " / " << reference.size() << endl;
    cout << "Overall accuracy:  " << accuracy << "%" << endl;
    calcSNPAccuracy(snps, reconstructed);
    cout << "Elapsed time:      " << elapsed << "ms" << endl;

    // save restored genome to file
    ofstream out("restored_genome_snp.txt");
    if (!out) { cerr << "Failed to save restored_genome_snp.txt" << endl; return 1; }
    out << reconstructed;
    out.close();
    cout << "Restored genome saved to restored_genome_snp.txt" << endl;

    // count existing trials for this SNP rate to assign trial number
    int trialNum = 1;
    ifstream logIn("sw_snp_log.txt");
    if (logIn) {
        string line;
        string snpRateStr = "SNP rate:         " + to_string((int)round(snpRate * 100)) + "%";
        while (getline(logIn, line))
            if (line.find(snpRateStr) != string::npos)
                trialNum++;
        logIn.close();
    }

    ofstream log("sw_snp_log.txt", ios::app);
    log << "=== Smith-Waterman ===" << endl;
    log << "SNP rate:         " << snpRate * 100 << "%" << endl;
    log << "Trial:            " << trialNum << endl;
    log << "Reference length: " << reference.size() << endl;
    log << "Read length:      " << readLen << endl;
    log << "Number of reads:  " << reads.size() << endl;
    log << "Number of SNPs:   " << snpTotal << endl;
    log << "Covered positions: " << covered << " / " << reference.size() << endl;
    log << "Overall accuracy:  " << accuracy << "%" << endl;
    log << "SNP accuracy:      " << snpCorrect << " / " << snpTotal
        << " (" << (snpTotal ? (double)snpCorrect/snpTotal*100.0 : 0.0) << "%)" << endl;
    log << "Elapsed time:      " << elapsed << "ms" << endl;
    log << endl;
    log.close();
    cout << "Log saved to sw_snp_log.txt (Trial " << trialNum << ")" << endl;

    return 0;
}
