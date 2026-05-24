#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

using namespace std;

// 두 read 간 overlap 길이 계산
int overlap(const string& a, const string& b) {
    int maxOverlap = 0;
    int minLen = min(a.size(), b.size());
    for (int len = 1; len <= minLen; len++) {
        if (a.substr(a.size() - len) == b.substr(0, len)) {
            maxOverlap = len;
        }
    }
    return maxOverlap;
}

// Greedy하게 가장 큰 overlap을 가진 read를 이어붙임
string assemble(vector<string> reads) {
    while (reads.size() > 1) {
        int bestA = 0, bestB = 0, bestOverlap = -1;
        string merged;

        // 모든 쌍 비교
        for (int i = 0; i < reads.size(); i++) {
            for (int j = 0; j < reads.size(); j++) {
                if (i == j) continue;
                int ov = overlap(reads[i], reads[j]);
                if (ov > bestOverlap) {
                    bestOverlap = ov;
                    bestA = i;
                    bestB = j;
                    merged = reads[i] + reads[j].substr(ov);
                }
            }
        }

        // 가장 좋은 쌍을 병합
        vector<string> newReads;
        for (int i = 0; i < reads.size(); i++) {
            if (i != bestA && i != bestB) newReads.push_back(reads[i]);
        }
        newReads.push_back(merged);
        reads = newReads;
    }
    return reads[0];
}

int main() {
    vector<string> reads;
    ifstream infile("reads_shuffled.txt");
    string line;
    while (getline(infile, line)) {
        if (!line.empty()) {
            reads.push_back(line);
        }
    }
    infile.close();
    string genome = assemble(reads);
    cout << "===== Assembled Genome (String Graph style =====" << endl;
    cout << genome << endl;

    return 0;
}
