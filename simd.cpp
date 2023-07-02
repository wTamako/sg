#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <immintrin.h> 

using namespace std;

void readMatrixData(const string& filePath, int N, vector<vector<int>>& matrix) {
    ifstream file;
    file.open(filePath, ios::binary);

    for (int i = 0; i < N; i++) {
        matrix.push_back(vector<int>(N, 0));
        int col;
        char ch = ' ';
        file >> col;
        int r = N - 1 - col;
        matrix[i][r >> 5] = 1 << (31 - (r & 31));
        file.get(ch);
        while (file.peek() != '\r') {
            file >> col;
            int diff = N - 1 - col;
            matrix[i][diff >> 5] += 1 << (31 - (diff & 31));
            file.get(ch);
        }
    }

    file.close();
}

void processMatrix(vector<vector<int>>& EE, vector<vector<int>>& ER, vector<int>& flag) {
    int N = EE[0].size();
    for (int i = 0; i < EE.size(); i++) {
        int byte = 0;
        int bit = 0;
        while (true) {
            while (byte < N && EE[i][byte] == 0) {
                byte++;
                bit = 0;
            }
            if (byte >= N) {
                break;
            }
            int temp = EE[i][byte] << bit;
            while (temp >= 0) {
                bit++;
                temp <<= 1;
            }
            int& isExist = flag[N - 1 - (byte << 5) - bit];
            if (!isExist == 0) {
                vector<int>& er = isExist > 0 ? ER[isExist - 1] : EE[~isExist];
                for (int j = 0; j + 8 <= N; j += 8) {
                    __m256i vaki = _mm256_loadu_si256((__m256i*)&EE[i][j]);
                    __m256i vakj = _mm256_loadu_si256((__m256i*)&er[j]);
                    vaki = _mm256_xor_si256(vaki, vakj);
                    _mm256_storeu_si256((__m256i*)&EE[i][j], vaki);
                    if (j + 16 > N) { // 处理串行末尾
                        while (j < N) {
                            EE[i][j] ^= er[j];
                            j++;
                        }
                        break;
                    }
                }
            } else {
                isExist = ~i;
                break;
            }
        }
    }
}

void writeResult(const string& filePath, int N, vector<vector<int>>& EE) {
    ofstream resFile(filePath, ios::trunc);
    for (int i = 0; i < EE.size(); i++) {
        int count = N - 1;
        for (int j = 0; j < EE[i].size(); j++) {
            int dense = EE[i][j];
            while (dense != 0) {
                int bit = __builtin_ctz(dense);
                dense ^= 1 << bit;
                resFile << count - bit << " ";
            }
            count -= 32;
        }
        resFile << "\n";
    }
    resFile.close();
}

int main() {
    int N = 0; // 矩阵的宽度和高度
    string inputPath = "input.txt"; // 输入文件路径
    string outputPath = "output.txt"; // 输出文件路径

    cout << "请输入矩阵的宽度和高度N：";
    cin >> N;

    vector<vector<int>> EE;
    vector<vector<int>> ER(N, vector<int>(N, 0));
    vector<int> flag(N, 0);

    readMatrixData(inputPath, N, EE);
    processMatrix(EE, ER, flag);
    writeResult(outputPath, N, EE);

    cout << "处理完成，结果已保存至output.txt文件中。" << endl;

    return 0;
}
