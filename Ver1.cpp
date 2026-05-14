#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

struct MatchResult {
    int bestScore;
    int bestPos;
};

// Smith-Waterman으로
// read가 genome 어디와 가장 유사한지 탐색
MatchResult smithWaterman(const string& genome,
    const string& read)
{
    int n = genome.size();
    int m = read.size();

    vector<vector<int>> dp(
        m + 1,
        vector<int>(n + 1, 0)
    );

    int match = 2;
    int mismatch = -1;
    int gap = -1;

    int bestScore = 0;
    int bestEndPos = 0;

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {

            int diag = dp[i - 1][j - 1] +
                (read[i - 1] == genome[j - 1]
                    ? match
                    : mismatch);

            int up = dp[i - 1][j] + gap;
            int left = dp[i][j - 1] + gap;

            dp[i][j] = max({
                0,
                diag,
                up,
                left
                });

            if (dp[i][j] > bestScore) {
                bestScore = dp[i][j];
                bestEndPos = j - 1;
            }
        }
    }

    // read 길이 기준으로 시작 위치 추정
    int estimatedStart =
        bestEndPos - m + 1;

    if (estimatedStart < 0)
        estimatedStart = 0;

    return {
        bestScore,
        estimatedStart
    };
}

int main() {

    // =========================
    // 1. reference 읽기
    // =========================

    ifstream refFile("reference.txt");

    if (!refFile.is_open()) {
        cout << "reference.txt 열기 실패" << endl;
        return 1;
    }

    string genome;
    getline(refFile, genome);

    refFile.close();

    // =========================
    // 2. reads 읽기
    // =========================

    ifstream readFile("reads_shuffled.txt");

    if (!readFile.is_open()) {
        cout << "reads_shuffled.txt 열기 실패" << endl;
        return 1;
    }

    vector<string> reads;

    string line;

    while (getline(readFile, line)) {

        if (!line.empty()) {
            reads.push_back(line);
        }
    }

    readFile.close();

    // =========================
    // 3. 복원용 문자열 생성
    // =========================

    // 처음에는 전부 '?'
    string restored(
        genome.size(),
        '?'
    );

    // =========================
    // 4. 각 read 정렬
    // =========================

    cout << "===== Matching Start =====" << endl;

    for (int r = 0; r < reads.size(); r++) {

        string read = reads[r];

        MatchResult result =
            smithWaterman(genome, read);

        cout << "Read " << r
            << " -> Position "
            << result.bestPos
            << " | Score "
            << result.bestScore
            << endl;

        // =========================
        // 5. 복원 문자열에 삽입
        // =========================

        for (int i = 0; i < read.size(); i++) {

            int genomePos =
                result.bestPos + i;

            if (genomePos >= genome.size())
                break;

            // 비어있으면 바로 삽입
            if (restored[genomePos] == '?') {
                restored[genomePos] = read[i];
            }
            else {
                // 이미 문자가 존재하면
                // reference와 더 잘 맞는 쪽 유지

                if (restored[genomePos]
                    != genome[genomePos])
                {
                    restored[genomePos] =
                        read[i];
                }
            }
        }
    }

    // =========================
    // 6. 빈 공간 보정
    // =========================

    // 현재 프로젝트 조건상
    // overlap이 거의 없으므로
    // 남은 '?'는 reference 기반 보정

    for (int i = 0; i < restored.size(); i++) {

        if (restored[i] == '?') {
            restored[i] = genome[i];
        }
    }

    // =========================
    // 7. 결과 출력
    // =========================

    ofstream outFile("restored_result.txt");

    if (!outFile.is_open()) {
        cout << "출력 파일 생성 실패" << endl;
        return 1;
    }

    outFile << restored;

    outFile.close();

    cout << endl;
    cout << "===== Restoration Complete =====" << endl;
    cout << "Result saved to restored_result.txt" << endl;

    cout << endl;
    cout << "Restored Genome:" << endl;
    cout << restored << endl;

    return 0;
}