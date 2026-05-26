#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <ctime>

using namespace std;

class GenomeSeedGenerator {
private:
    mt19937 rng;
    const string BASES = "ATGC";

public:
    GenomeSeedGenerator() {
        rng.seed(time(nullptr));
    }

    // 랜덤 염기 생성
    char randomBase() {
        uniform_int_distribution<int> dist(0, 3);
        return BASES[dist(rng)];
    }

    // 길이 N의 랜덤 게놈 생성
    string generateReferenceGenome(int N) {
        string genome;
        genome.reserve(N);

        for (int i = 0; i < N; i++) {
            genome += randomBase();
        }

        return genome;
    }

    // SNP 적용
    string generateVariantGenome(
        const string& reference,
        double snpRate,
        vector<string>& snpLogs
    ) {
        string variant = reference;

        uniform_real_distribution<double> prob(0.0, 1.0);

        for (int i = 0; i < variant.size(); i++) {

            // SNP 발생
            if (prob(rng) < snpRate) {

                char original = variant[i];
                char mutated;

                // 기존 문자 제외하고 변경
                do {
                    mutated = randomBase();
                } while (mutated == original);

                variant[i] = mutated;

                // 로그 기록
                string log =
                    to_string(i) + " " +
                    original + " " +
                    mutated;

                snpLogs.push_back(log);
            }
        }

        return variant;
    }

    // Read 생성
    vector<string> generateReads(
        const string& genome,
        int L,
        int minOverlap,
        int maxOverlap
    ) {
        vector<string> reads;

        uniform_int_distribution<int> overlapDist(
            minOverlap,
            maxOverlap
        );

        int start = 0;

        while (true) {

            // 마지막 read 처리
            if (start + L >= genome.size()) {

                start = genome.size() - L;

                if (start < 0)
                    start = 0;

                string read = genome.substr(start, L);

                reads.push_back(read);

                break;
            }

            // 일반 read 생성
            string read = genome.substr(start, L);

            reads.push_back(read);

            int overlap = overlapDist(rng);

            // 다음 시작 위치
            start += (L - overlap);
        }

        return reads;
    }

    // 파일 저장
    void saveSingleLineFile(
        const string& data,
        const string& filename
    ) {
        ofstream fout(filename);

        fout << data;

        fout.close();
    }

    // 여러 줄 저장
    void saveVectorFile(
        const vector<string>& data,
        const string& filename
    ) {
        ofstream fout(filename);

        for (const string& line : data) {
            fout << line << '\n';
        }

        fout.close();
    }
};

int main() {

    GenomeSeedGenerator generator;

    // =========================
    // 입력
    // =========================

    int N;
    int L;

    double snpRate;

    int minOverlap;
    int maxOverlap;

    cout << "Genome Length (N): ";
    cin >> N;

    cout << "Read Length (L): ";
    cin >> L;

    cout << "SNP Rate (0.0 ~ 1.0): ";
    cin >> snpRate;

    cout << "Min Overlap: ";
    cin >> minOverlap;

    cout << "Max Overlap: ";
    cin >> maxOverlap;

    // =========================
    // 예외 처리
    // =========================

    if (L > N) {
        cout << "Error: L cannot be greater than N.\n";
        return 1;
    }

    if (minOverlap > maxOverlap) {
        cout << "Error: minOverlap > maxOverlap\n";
        return 1;
    }

    if (maxOverlap >= L) {
        cout << "Error: maxOverlap must be smaller than L.\n";
        return 1;
    }

    // =========================
    // STEP 1
    // Reference Genome 생성
    // =========================

    string referenceGenome =
        generator.generateReferenceGenome(N);

    // =========================
    // STEP 2
    // Variant Genome 생성
    // =========================

    vector<string> snpLogs;

    string variantGenome =
        generator.generateVariantGenome(
            referenceGenome,
            snpRate,
            snpLogs
        );

    // =========================
    // STEP 3
    // Read 생성
    // =========================

    vector<string> orderedReads =
        generator.generateReads(
            variantGenome,
            L,
            minOverlap,
            maxOverlap
        );

    // shuffle용 복사
    vector<string> shuffledReads = orderedReads;

    shuffle(
        shuffledReads.begin(),
        shuffledReads.end(),
        mt19937(random_device()())
    );

    // =========================
    // 파일 저장
    // =========================

    generator.saveSingleLineFile(
        referenceGenome,
        "Seed_Data\\reference.txt"
    );

    generator.saveSingleLineFile(
        variantGenome,
        "Seed_Data\\variant.txt"
    );

    generator.saveVectorFile(
        orderedReads,
        "Seed_Data\\reads_ordered.txt"
    );

    generator.saveVectorFile(
        shuffledReads,
        "Seed_Data\\reads_shuffled.txt"
    );

    generator.saveVectorFile(
        snpLogs,
        "Seed_Data\\snp_log.txt"
    );

    // =========================
    // 완료 메시지
    // =========================

    cout << "\n===== Generation Complete =====\n";

    cout << "reference.txt generated\n";
    cout << "variant.txt generated\n";
    cout << "reads_ordered.txt generated\n";
    cout << "reads_shuffled.txt generated\n";
    cout << "snp_log.txt generated\n";

    cout << "\nTotal Reads: "
        << orderedReads.size()
        << '\n';

    cout << "Total SNPs: "
        << snpLogs.size()
        << '\n';

    return 0;
}