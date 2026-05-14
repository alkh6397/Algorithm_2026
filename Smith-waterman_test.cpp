#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

struct Result {
    int bestScore;
    int endPos; // genome 기준 끝 위치
};

Result smithWaterman(const string& genome, const string& read) {
    int n = genome.size();
    int m = read.size();

    // DP 테이블
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));

    // 점수 설정
    int match = 2;
    int mismatch = -1;
    int gap = -1;

    int bestScore = 0;
    int bestCol = 0;

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {

            int diag =
                dp[i - 1][j - 1] +
                (read[i - 1] == genome[j - 1] ? match : mismatch);

            int up = dp[i - 1][j] + gap;
            int left = dp[i][j - 1] + gap;

            dp[i][j] = max({
                0,
                diag,
                up,
                left
                });

            // 최고 점수 갱신
            if (dp[i][j] > bestScore) {
                bestScore = dp[i][j];
                bestCol = j;
            }
        }
    }

    return { bestScore, bestCol - 1 };
}

int main() {

    string genome =
        "ATCGTACGATCGGATCGTTAACCG";

    string read =
        "GATAGTT";

    Result result = smithWaterman(genome, read);

    cout << "Best Score : " << result.bestScore << endl;
    cout << "Best Match End Position : "
        << result.endPos << endl;

    // 대략적인 시작 위치 추정
    int startPos =
        max(0, result.endPos - (int)read.size() + 1);

    cout << "Approx Start Position : "
        << startPos << endl;

    cout << "Matched Region : ";

    for (int i = startPos;
        i <= result.endPos && i < genome.size();
        i++) {
        cout << genome[i];
    }

    cout << endl;

    return 0;
}