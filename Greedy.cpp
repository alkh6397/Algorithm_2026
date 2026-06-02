#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;

const int MIN_OVERLAP = 5;

//두 문자열의 최대 overlap 길이 계산 함수
int get_overlap(const string& a, const string& b) {

	int max_overlap = 0;

	int minlen = min(a.size(), b.size());
	for (int len = MIN_OVERLAP; len <= minlen; len++) {

		//a의 suffix == b의 prefix ??
		if (a.substr(a.size() - len) == b.substr(0, len)) {
			max_overlap = len;
		}
	}

	return max_overlap;
}

int main() {

	//reads 저장용 vector
	vector <string> reads;

	//파일 열기
	ifstream fin;
	fin.open("C:\\Users\\유성\\Desktop\\Seed_Creater\\Seed_Data\\reads_shuffled.txt");

	//파일 열기 실패시
	if (!fin.is_open()) {
		cout << "파일 열기 실패\n";
		return 1;
	}

	string line;

	//reads vector에 reads 값들을 저장
	while (getline(fin, line)) {

		if (line.empty()) continue;

		reads.push_back(line);
	}

	fin.close();

	//초기 Read 길이 측정용 String
	string L = reads[0];

	cout << "총 읽어온 reads 갯수: " << reads.size() << endl;
	cout << "탐색을 시작합니다. 잠시만 기다려주세요 ...\n" << endl;

	int total_reads = reads.size();

	auto start_time = high_resolution_clock::now();

	while (reads.size() > 1) {

		int bestI = -1;
		int bestJ = -1;
		int best_overlap = -1;

		for (int i = 0; i < reads.size(); i++) {
			for (int j = 0; j < reads.size(); j++) {

				if (i == j) continue;

				int overlap = get_overlap(reads[i], reads[j]);

				if (overlap < MIN_OVERLAP) {
					continue;
				}

				if (overlap > best_overlap) {
					best_overlap = overlap;
					bestI = i;
					bestJ = j;
				}
			}
		}

		if (bestI == -1 || bestJ == -1) {
			cout << "\n[알림] 더 이상 " << MIN_OVERLAP << "bp 이상 겹치는 조각이 없어 조립을 조기 종료합니다." << endl;
			break;
		}

		string merged = reads[bestI] + reads[bestJ].substr(best_overlap);

		if (reads.size() % 10 == 0) {
			cout << "\r[진행 상황] 남은 조각: " << reads.size() << " / " << total_reads
				<< " | 현재 결합된 최대 길이: " << merged.size() << " bp    " << flush;
		}

		//더 큰 값을 뒤쪽 인덱스로 옮기기
		if (bestI > bestJ)
			swap(bestI, bestJ);

		//인덱스 혼란을 막기 위해 뒤쪽 인덱스부터 삭제
		reads.erase(reads.begin() + bestJ);
		reads.erase(reads.begin() + bestI);

		//결합한 read 추가하기
		reads.push_back(merged);

	}

	auto end_time = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(end_time - start_time);

	cout << "\r[진행 상황] 남은 조각: 1 / " << total_reads
		<< " | 현재 결합된 최대 길이: " << reads[0].size() << " bp          " << flush;

	cout << "\n\n조립 완료!\n" << endl;

	/*
	cout << "\nFinal Genome: \n";
	cout << reads[0] << "\n\n" << endl;
	*/


	//조합된 reads를 string형식으로 저장
	string Merged_reads = reads[0];

	//조합해야할 목표 Genome 파일 open
	ifstream Goal_Gene;
	Goal_Gene.open("C:\\Users\\유성\\Desktop\\Seed_Creater\\Seed_Data\\variant.txt");

	if (!Goal_Gene.is_open()) {
		cout << "파일 열기 실패\n";
		return 1;
	}

	//목표 Genome을 저장할 GoalGene string
	string GoalGene, line2;

	while (getline(Goal_Gene, line2)) {

		if (line2.empty()) continue;

		GoalGene += line2;
	}

	Goal_Gene.close();

	/*
	cout << "Target Genome: \n" << GoalGene << "\n\n" << endl;
	*/

	int diff = 0;	int match_count = 0;
	int min_len = min(GoalGene.size(), Merged_reads.size());

	//목표 Genome과 조합한 Reads 비교용 반복문. 다르다면 다른 인덱스와 원소 출력
	for (int i = 0; i < min_len; i++) {
		if (Merged_reads[i] != GoalGene[i]) {
			diff++;
		}
		else {
			match_count++;
		}
	}

	int length_diff = abs((int)GoalGene.size() - (int)Merged_reads.size());
	diff += length_diff;

	double accuracy = (double)match_count / GoalGene.size() * 100.0;

	if (diff == 0 && length_diff == 0) {
		cout << "원본과 일치합니다.\n" << endl;
	}
	else {
		cout << "총 " << diff << "개의 염기가 다르게 조합되었습니다.(길이 차이 포함)\n" << endl;
	}

	cout << "전체 Genome의 길이(N): " << GoalGene.size() << endl;
	cout << "Reads의 개수: " << total_reads << endl;
	cout << "Read의 길이(L): " << L.size() << endl;
	cout << "SNP 비율: 0.05\n" << endl;

	cout << "최종 조합 정확도: " << accuracy << "%" << endl;
	cout << "**알고리즘 실행 시간: " << duration.count() << " ms (" << duration.count() / 1000.0 << " 초)**" << endl;

	return 0;
}
