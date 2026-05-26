#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <ctime>

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
    clock_t start, finish;
    double duration;
    start = clock();

    vector<string> reads;
    ifstream infile("reads_shuffled.txt");
    string line;
    while (getline(infile, line)) {
        if (!line.empty()) {
            reads.push_back(line);
        }
    }
    infile.close();
    // 예측한 정답
    string genome = assemble(reads);
    cout << "===== Assembled Genome (String Graph style =====" << endl;
    cout << genome << endl;

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    cout << duration << "초" << endl;

    // assembled_genome.txt 남겨놓기
    ofstream outfile("assembled_genome.txt");
    if (outfile.is_open()) {
        outfile << genome;
        outfile.close();
    }
    else {
        cout << "파일을 열 수 없습니다." << endl;
    }

    // 실제 정답
    string answer_genome;
    ifstream infile_answer("variant.txt");
    getline(infile_answer, answer_genome);
    infile_answer.close();

    // 정확도 계산
    int len1 = genome.length();
    int len2 = answer_genome.length();
    int min_len = min(len1, len2);
    int mismatch = 0;

    for (int i = 0; i < min_len; i++) {
        if (genome[i] != answer_genome[i])
            mismatch++;
    }
    int total = min_len;
    double accuracy = ((double)(total - mismatch) / total) * 100.0;
    cout << "Mismatch count: " << mismatch << endl;
    cout << "Accuracy: " << accuracy << "%" << endl;
    

    return 0;
}