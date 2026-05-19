#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <random>
#include <ctime>

using namespace std;

struct Individual {
    vector<int> order; // read 순서
    int fitness;       // overlap 점수
};

vector<string> reads = {
    "ATCGTAAGGCCT",
    "TACCTGAACGTT",
    "CTTACGATTCGA",
    "TTAGCCATGGAT",
    "CATGGATCGTAA",//
    "GAACGTTAGCCA",
    "AAGGCCTTACGA",
    "ATGCGTACCTGA"
};

mt19937 rng(time(0));


// overlap 길이 계산
const int MIN_OVERLAP = 4;

int overlap(const string& a, const string& b) {

    int maxOverlap = 0;

    int minLen = min(a.size(), b.size());

    for (int len = MIN_OVERLAP; len <= minLen; len++) {

        string suffix = a.substr(a.size() - len);
        string prefix = b.substr(0, len);

        if (suffix == prefix)
            maxOverlap = len;
    }

    return maxOverlap;
}


// fitness 계산
int calculateFitness(const vector<int>& order) {

    int score = 0;

    for (int i = 0; i < order.size() - 1; i++) {

        string a = reads[order[i]];
        string b = reads[order[i + 1]];

        int ov = overlap(a, b);

        if (ov == 0)
            score -= 100;
        else
            score += ov;
    }

    return score;
}


// genome 조립
string assembleGenome(const vector<int>& order) {

    string result = reads[order[0]];

    for (int i = 1; i < order.size(); i++) {

        string prev = reads[order[i - 1]];
        string curr = reads[order[i]];

        int ov = overlap(prev, curr);

        // overlap 제외 부분만 추가
        result += curr.substr(ov);
    }

    return result;
}


// 랜덤 개체 생성
Individual createRandomIndividual(int n) {

    Individual ind;

    ind.order.resize(n);

    for (int i = 0; i < n; i++)
        ind.order[i] = i;

    shuffle(ind.order.begin(), ind.order.end(), rng);

    ind.fitness = calculateFitness(ind.order);

    return ind;
}


// 돌연변이
void mutate(Individual& ind) {

    int a = rng() % ind.order.size();
    int b = rng() % ind.order.size();

    swap(ind.order[a], ind.order[b]);

    ind.fitness = calculateFitness(ind.order);
}


int main() {

    const int POP_SIZE = 100;
    const int GENERATIONS = 1000;

    vector<Individual> population;

    // 초기 개체 생성
    for (int i = 0; i < POP_SIZE; i++) {
        population.push_back(createRandomIndividual(reads.size()));
    }

    // 진화 시작
    for (int gen = 0; gen < GENERATIONS; gen++) {

        // fitness 기준 정렬
        sort(population.begin(), population.end(),
            [](const Individual& a, const Individual& b) {
                return a.fitness > b.fitness;
            });

        // 상위 개체 출력
        if (gen % 100 == 0) {
            cout << "Generation " << gen
                << " Best Fitness = "
                << population[0].fitness << "\n";
        }

        // 상위 절반 생존
        population.resize(POP_SIZE / 2);

        // 복제 + 돌연변이
        while (population.size() < POP_SIZE) {

            Individual child =
                population[rng() % (POP_SIZE / 2)];

            mutate(child);

            population.push_back(child);
        }
    }

    // 최종 정렬
    sort(population.begin(), population.end(),
        [](const Individual& a, const Individual& b) {
            return a.fitness > b.fitness;
        });

    // 최고 개체
    Individual best = population[0];

    cout << "\n===== BEST ORDER =====\n";

    for (int idx : best.order)
        cout << idx << " ";

    cout << "\n";

    // 복원 결과
    string genome = assembleGenome(best.order);

    cout << "\n===== ASSEMBLED GENOME =====\n";
    cout << genome << "\n";

    return 0;
}