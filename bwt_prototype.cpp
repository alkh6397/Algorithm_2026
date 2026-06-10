#include<iostream>
#include<string>
#include<vector>
#include<algorithm>
#include<map>
#include<fstream>
#include<chrono>
using namespace std;
using namespace std::chrono;

map<char, int> C;
vector<vector<int>> occ;

int base_to_idx(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 4; // $ 처리
    }
}

//bwt 추출 함수
pair<string, vector<int>> generateBWT(string text) {
    text += "$";
    int n = text.length();

    vector<int> indices(n);
    for(int i = 0; i < n; i++) {
        indices[i] = i;
    }

    sort(indices.begin(), indices.end(), [&](int a, int b) {
        for (int k = 0; k < n; k++) {
            char char_a = text[(a + k) % n];
            char char_b = text[(b + k) % n];
            if (char_a != char_b) {
                return char_a < char_b;
            }
        }
        return false;
    });

    string L = "";
    vector<int> SA;
    L.reserve(n); 
    SA.reserve(n);

    for(int i = 0; i < n; i++) {
        int start_idx = indices[i];
        L += text[(start_idx + n - 1) % n]; 
        SA.push_back(start_idx); 
    }

    return pair(L, SA);
}

//재귀 탐색 진행 함수
void search_bwt(const string& read, int read_idx, int top, int bottom, int mismatches_left, 
                vector<int>& matched_rows) {
    
    //끝까지 탐색 성공
    if (read_idx < 0) {
        for (int i = top; i <= bottom; i++) {
            matched_rows.push_back(i);
        }
        return;
    }

    // A, C, G, T 4방향으로 탐
    char bases[] = {'A', 'C', 'G', 'T'};
    for (char c : bases) {
        int c_idx = base_to_idx(c);
        
        int new_top = C[c] + occ[top][c_idx];
        int new_bottom = C[c] + occ[bottom + 1][c_idx] - 1;

        if (new_top <= new_bottom) {
            int cost = (read[read_idx] == c) ? 0 : 1;
            
            //잔여 허용 오차가 남아있으면 재귀로 탐색 진행
            if (mismatches_left >= cost) {
                search_bwt(read, read_idx - 1, new_top, new_bottom, mismatches_left - cost, matched_rows);
            }
        }
    }
}


int main(void) {
    // 원본 염기서열 입력부
    string reference = "";
    ifstream fin_ref("Seed_Data/reference.txt");
    if (!fin_ref.is_open()) {
        cout << "에러: 파일을 찾을 수 없습니다!" << '\n';
    }
    fin_ref >> reference; 
    fin_ref.close();

    // read들 입력
    vector<string> reads;
    ifstream fin_rd("Seed_Data/reads_shuffled.txt");
    if (!fin_rd.is_open()) {
        cout << "에러: 파일을 찾을 수 없습니다!" << '\n';
    }

    string line;
    while (getline(fin_rd, line)) {
        if (!line.empty()) {
            reads.push_back(line);
        }
    }
    fin_rd.close();

    // 허용하는 최대 오차 수

    // 선형으로 설정
    int MAX_MISMATCH = (int) ((reads[0].length() * 0.05) + 1);

    // 상수 설정
    //int MAX_MISMATCH = 6;
    
    //====================================================================================//
    auto start = high_resolution_clock::now();

    pair<string, vector<int>> bwtResult = generateBWT(reference);
    string L = bwtResult.first, F = bwtResult.first; sort(F.begin(), F.end());
    vector<int> SA = bwtResult.second;

    C.clear(); // 전역변수 초기화
    for(int i = 0; i < F.length(); i++) {
        if(C.find(F[i]) == C.end()) {
            C[F[i]] = i;
        }
    }

    //occ배열 값 추출
    int l_len = L.length();
    occ.assign(l_len + 1, vector<int>(5, 0)); // 전역변수 할당

    //BWT 문자열을 0번 칸부터 한 글자씩 읽어가며 누적합(Prefix Sum) 계산
    for (int i = 0; i < l_len; i++) {
        //이전 칸의 누적 데이터들을 그대로 복사
        for (int j = 0; j < 5; j++) {
            occ[i + 1][j] = occ[i][j];
        }
        //현재 인덱스(i)에 등장한 문자의 누적 카운트를 1 증가시키기
        int current_char_idx = base_to_idx(L[i]);
        occ[i + 1][current_char_idx]++;
    }

    //투표함 만들기
    vector<map<char, int>> pileup(reference.length());

    //조각(reads)들을 하나씩 꺼내서 탐색 수행
    int time_count = 0;
    for (string read : reads) {
        vector<int> matched_rows; 

        search_bwt(read, read.length() - 1, 0, L.length() - 1, MAX_MISMATCH, matched_rows);

        time_count++;
        cout << "\r" << " (" << time_count << "/" << reads.size() << ")" << "개 조각 탐색 완료..." << flush;

        //찾은 줄 번호를 진짜 위치(SA)로 변환하고 투표하기
        for (int row : matched_rows) {
            int original_pos = SA[row]; 

            for (int i = 0; i < read.length(); i++) {
                if (original_pos + i < reference.length()) {
                    pileup[original_pos + i][read[i]]++;
                }
            }
        }
    }
    cout << endl;

    //투표 결과 확인
    string final_genome = "";
    for (int i = 0; i < reference.length(); i++) {
        char best_base = reference[i]; 
        int max_votes = 0;

        for (auto const& [base, votes] : pileup[i]) {
            if (votes > max_votes) {
                max_votes = votes;
                best_base = base; 
            }
        }
        final_genome += best_base; 
    }

    //분석부
    string variant = "";
    ifstream fin_var("Seed_Data/variant.txt");
    if (!fin_var.is_open()) {
        cout << "에러: 파일을 찾을 수 없습니다!" << '\n';
    }
    fin_var >> variant; 
    fin_var.close();

    int cnt = 0;
    for(int i = 0; i < variant.length(); i++) {
        if(final_genome[i] == variant[i]) {
            cnt++;
        }
    }
    double accuracy = (double)cnt / variant.length();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    //결과 출력
    cout << "\n=====================================" << endl;
    cout << "적용된 Max_Mismatch: " << MAX_MISMATCH << endl;
    cout << "정확도: " << accuracy << endl;
    cout << "걸린 시간(초): " << duration.count() / 1000.0 << "초\n"; 
    cout << "=====================================" << endl;

    return 0;
}
