#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <iomanip>

using namespace std;

// 최소 overlap 길이
const int MIN_OVERLAP = 5;

// overlap 계산
int overlap(const string& a, const string& b) {

    int maxOverlap = 0;

    int minLen = min(a.size(), b.size());

    for (int len = MIN_OVERLAP; len <= minLen; len++) {

        if (a.substr(a.size() - len) == b.substr(0, len)) {
            maxOverlap = len;
        }
    }

    return maxOverlap;
}

// 현재 read 뒤로 얼마나 더 연결 가능한지 탐색
int futureScore(const vector<string>& reads, int current) {

    int best = 0;

    for (int k = 0; k < reads.size(); k++) {

        if (k == current)
            continue;

        int ov = overlap(reads[current], reads[k]);

        best = max(best, ov);
    }

    return best;
}

// 개선된 Greedy Assembly
string assemble(vector<string> reads) {

    int initialSize = reads.size();

    while (reads.size() > 1) {

        int bestA = -1;
        int bestB = -1;

        int bestScore = -1;

        string merged;

        // 모든 쌍 비교
        for (int i = 0; i < reads.size(); i++) {

            for (int j = 0; j < reads.size(); j++) {

                if (i == j)
                    continue;

                // overlap 계산
                int ov = overlap(reads[i], reads[j]);

                // overlap 너무 작으면 무시
                if (ov < MIN_OVERLAP)
                    continue;

                // 미래 연결 가능성
                int future = futureScore(reads, j);

                // 총 점수
                int totalScore = ov + future;

                // 가장 좋은 경로 선택
                if (totalScore > bestScore) {

                    bestScore = totalScore;

                    bestA = i;
                    bestB = j;

                    merged =
                        reads[i] +
                        reads[j].substr(ov);
                }
            }
        }

        // 더 이상 연결 불가능
        if (bestA == -1 || bestB == -1) {

            cout << "\nNo more valid overlaps found.\n";

            break;
        }

        // merge 수행
        vector<string> newReads;

        for (int i = 0; i < reads.size(); i++) {

            if (i != bestA && i != bestB) {
                newReads.push_back(reads[i]);
            }
        }

        newReads.push_back(merged);

        reads = newReads;

        // 진행률 계산
        double progress =
            ((double)(initialSize - reads.size()) /
                (initialSize - 1)) * 100.0;

        // 진행률 출력
        cout << "\rProgress: "
            << fixed << setprecision(2)
            << progress << "%   "
            << flush;
    }

    cout << endl;

    return reads[0];
}

int main() {

    clock_t start, finish;

    double duration;

    start = clock();

    vector<string> reads;

    ifstream infile("C:\\Users\\alkh6\\source\\repos\\SeedCreator\\x64\\Debug\\Seed_Data\\reads_shuffled.txt");

    string line;

    while (getline(infile, line)) {

        if (!line.empty()) {
            reads.push_back(line);
        }
    }

    infile.close();

    cout << "Read count: "
        << reads.size() << endl;

    cout << "Minimum overlap: "
        << MIN_OVERLAP << endl;

    cout << "Assembling genome..." << endl;

    // 조립 수행
    string genome = assemble(reads);

    cout << endl;

    cout << "===== Assembled Genome =====" << endl;

    cout << genome << endl;

    finish = clock();

    duration =
        (double)(finish - start) / CLOCKS_PER_SEC;

    cout << endl;

    cout << "Running time: "
        << fixed << setprecision(3)
        << duration << " sec" << endl;

    // assembled genome 저장
    ofstream outfile("assembled_genome.txt");

    if (outfile.is_open()) {

        outfile << genome;

        outfile.close();
    }

    // 실제 정답 genome 불러오기
    string answer_genome;

    ifstream infile_answer("C:\\Users\\alkh6\\source\\repos\\SeedCreator\\x64\\Debug\\Seed_Data\\variant.txt");

    getline(infile_answer, answer_genome);

    infile_answer.close();

    // 정확도 계산
    int len1 = genome.length();

    int len2 = answer_genome.length();

    int min_len = min(len1, len2);

    int mismatch = 0;

    for (int i = 0; i < min_len; i++) {

        if (genome[i] != answer_genome[i]) {
            mismatch++;
        }
    }

    // 길이 차이 반영
    mismatch += abs(len1 - len2);

    int total = max(len1, len2);

    double accuracy =
        ((double)(total - mismatch) / total) * 100.0;

    cout << "Mismatch count: "
        << mismatch << endl;

    cout << "Accuracy: "
        << fixed << setprecision(2)
        << accuracy << "%" << endl;

    // 결과 저장
    ofstream resultfile("result.txt");

    if (resultfile.is_open()) {

        resultfile << "Mismatch count: "
            << mismatch << endl;

        resultfile << "Accuracy: "
            << fixed << setprecision(2)
            << accuracy << "%" << endl;

        resultfile << "Running time: "
            << fixed << setprecision(3)
            << duration << " sec" << endl;

        resultfile.close();
    }

    return 0;
}