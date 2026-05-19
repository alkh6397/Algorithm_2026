#include<iostream>
#include<string>
#include<vector>
#include<algorithm>
#include<map>
#include<fstream>
#include<chrono>
using namespace std;
using namespace std::chrono;

//LF매핑, SA값 반환
pair<string, vector<int>> generateBWT(string text) {
    text += "$";
    int n = text.length();

    vector<pair<string, int>> rotations; // (pair로 선언하여 int값에는 각각의 BWT문자열이 몇번 째 글자(index)부터 읽은 문자열이였는지를 저장)
    for(int i = 0;i<n;i++) {
        string rotated = text.substr(n-1-i,n) + text.substr(0,n-1-i);
        rotations.push_back({rotated,n-1-i});
    }

    sort(rotations.begin(),rotations.end());

    // for(int i = 0;i<n;i++) {
    //     cout << rotations[i].first << '\n';
    // }

    string L; // BWT행렬의 LAST만 모아놓은 문자열
    for(int i = 0;i<rotations.size();i++) {
        L += rotations[i].first[n-1];
    }
    vector<int> SA;
    for(int i = 0;i<rotations.size();i++) {
        SA.push_back(rotations[i].second); // 표지판
    }

    // cout << L;
    return pair(L, SA);
}

//각 문자의 (bwt=L 문자열에서의) 개수 반한 함수
int get_occ(string L, char c, int idx) {
    int cnt = 0;
    for(int i = 0;i<=idx;i++) {
        if(L[i] == c) cnt++;
    }
    return cnt;
}


// 매개변수 설명:
// read: 우리가 찾을 조각 (예: "ATCCA")
// read_idx: 현재 찾고 있는 조각의 글자 인덱스 (맨 뒤부터 시작하므로 처음엔 read.length() - 1)
// top, bottom: 현재 좁혀진 탐색 범위
// mismatches_left: 남은 돌연변이(에러) 허용 횟수
// matched_rows: 최종적으로 조각을 찾은 줄 번호들을 담아갈 바구니 (참조자 & 사용)
void search_bwt(string read, int read_idx, int top, int bottom, 
                int mismatches_left, string bwt, map<char, int>& C, vector<int>& matched_rows) {
    
    // [종료 조건 1] 조각의 첫 글자(인덱스 0)까지 무사히 다 찾은 경우 (성공!)
    if (read_idx < 0) {
        // 현재 살아남은 top부터 bottom까지의 행 번호를 matched_rows에 모두 넣고 함수 종료
        for (int i = top; i <= bottom; i++) {
            matched_rows.push_back(i);
        }
        return;
    }

    // [탐색 진행] A, C, G, T 4방향 길을 모두 찔러보기A
    char bases[] = {'A', 'C', 'G', 'T'};
    for (char c : bases) {
        
        // 1. LF 매핑 공식으로 다음 스텝의 새로운 범위(new_top, new_bottom) 계산
        // (가이드 2에서 만든 get_occ 함수를 여기서 씁니다!)
        int new_top = C[c] + get_occ(bwt, c, top - 1);
        int new_bottom = C[c] + get_occ(bwt, c, bottom) - 1;

        // 2. 해당 글자로 가는 길이 존재하는지 확인 (범위가 뒤집히지 않았는지)
        if (new_top <= new_bottom) {
            
            // 3. 에러 티켓(비용) 계산
            // 내가 지금 찾으려는 조각의 글자(read[read_idx])와 
            // 현재 탐색하려는 길의 글자(c)가 다르면 에러 티켓 1장 소모!
            int cost;
            if (read[read_idx] == c) {
                cost = 0; // 일치하면 비용 0
            } else {
                cost = 1; // 다르면(돌연변이) 비용 1
            }
            
            // 4. 에러 티켓이 충분히 남아있다면 다음 글자로 계속 탐색 진행! (재귀 호출)
            if (mismatches_left >= cost) {
                // read_idx를 1 줄여서 앞 글자로 넘어가고, 범위는 new_top/new_bottom으로 갱신
                search_bwt(read, read_idx - 1, new_top, new_bottom, mismatches_left - cost, bwt, C, matched_rows);
            }
        }
    }
}


int main(void) {
    // 원본 염기서열 입력부 (후에 txt파일 추출한다던지 하는 방식으로 변경)
    string reference = "";
    
    ifstream fin_ref("Seed_Data/reference.txt");
    if (!fin_ref.is_open()) {
        cout << "에러: " << " 파일을 찾을 수 없습니다!" << '\n';
    }
    fin_ref >> reference; 
    fin_ref.close();

    // read들 입력
    vector<string> reads;
    
    ifstream fin_rd("Seed_Data/reads_shuffled.txt");
    
    if (!fin_rd.is_open()) {
        cout << "에러: " << " 파일을 찾을 수 없습니다!" << '\n';
    }

    string line;
    // 2. 파일의 끝(EOF)에 도달할 때까지 한 줄씩 계속 읽어오기
    while (getline(fin_rd, line)) {
        // 빈 줄이 아닐 경우에만 벡터에 추가
        if (!line.empty()) {
            reads.push_back(line);
        }
    }
    fin_rd.close();


    // 허용하는 최대 오차 수
    int MAX_MISMATCH = 1;
    
    //====================================================================================//
    auto start = high_resolution_clock::now();


    pair<string, vector<int>> bwtResult = generateBWT(reference);
    string L = bwtResult.first, F = bwtResult.first; sort(F.begin(), F.end());
    vector<int> SA = bwtResult.second;

    map<char, int> C; //각 문자의 시작 인덱스를 map으로 저장
    for(int i = 0;i<F.length();i++) {
        if(C.find(F[i]) == C.end()) {
            C[F[i]] = i;
        }
    }

    // 1. 투표함(Pile-up) 만들기
    // reference 길이만큼의 배열을 만들고, 각 칸마다 염기별 득표수를 담을 map을 둡니다.
    vector<map<char, int>> pileup(reference.length());

    // 2. 조각(reads)들을 하나씩 꺼내서 BWT 행렬(L열)에서 위치 찾기
    int time_count=0;
    for (string read : reads) {
        vector<int> matched_rows; // 이 조각이 매칭된 줄 번호들을 담을 벡터

        // 검색 시작! (가이드 3에서 만든 search_bwt 함수 호출)
        // 파라미터: 조각, 현재 인덱스(맨 뒤), top(0), bottom(L길이-1), 에러허용치, L열, C배열, 결과저장벡터
        search_bwt(read, read.length() - 1, 0, L.length() - 1, MAX_MISMATCH, L, C, matched_rows);

        time_count++;
        if(time_count%100 == 0) cout << time_count << "개 조각 탐색 완료...\n";

        // 3. 찾은 줄 번호(matched_rows)를 진짜 위치(SA)로 변환하고 투표하기
        for (int row : matched_rows) {
            int original_pos = SA[row]; // 표지판을 보고 원래 유전자의 인덱스 획득

            // 조각의 길이만큼 투표함에 득표수 더해주기
            for (int i = 0; i < read.length(); i++) {
                // 원본 서열의 범위를 벗어나지 않도록 안전 검사 (Out of bounds 방지)
                if (original_pos + i < reference.length()) {
                    pileup[original_pos + i][read[i]]++;
                }
            }
        }
    }

    // 4. 투표 결과(Consensus) 확인 및 최종 유전자 조립
    string final_genome = "";
    for (int i = 0; i < reference.length(); i++) {
        char best_base = reference[i]; // 기본값은 원본(reference)의 글자로 셋팅
        int max_votes = 0;

        // 투표함(pileup[i])을 열어보고 가장 표가 많은 염기 찾기
        for (auto const& [base, votes] : pileup[i]) {
            if (votes > max_votes) {
                max_votes = votes;
                best_base = base; // 다수결로 이긴 염기로 덮어쓰기 (SNP 변이 반영)
            }
        }
        final_genome += best_base; // 찾아낸 염기를 최종 유전자에 이어붙이기
    }


    //분석
    string variant = "";
    
    ifstream fin_var("Seed_Data/variant.txt");
    if (!fin_var.is_open()) {
        cout << "에러: " << " 파일을 찾을 수 없습니다!" << '\n';
    }
    fin_var >> variant; 
    fin_var.close();

    int cnt = 0;
    for(int i = 0;i<variant.length();i++) {
        if(final_genome[i] == variant[i]) {
            cnt++;
        }
    }
    double accuracy = (double)cnt / variant.length();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    // 5. 결과 출력
    cout << "\n=====================================" << endl;
    cout << "원본 레퍼런스: " << reference << endl;
    cout << "나의 최종 유전자: " << final_genome << endl;
    cout << "정확도: " << accuracy << endl;
    cout << "걸린 시간(초): " << duration.count() / 1000.0 << "초\n"; 
    cout << "=====================================" << endl;

    return 0;
}