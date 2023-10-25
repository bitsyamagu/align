#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <pthread.h>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <zlib.h> // For gzipped file reading

static void* computeScore(void* arg);

using namespace std;

class Fasta {
public:
    std::string name;
    std::string seq;
    void setName(std::string name) {
        this->name = name;
    }
    void setSequence(std::string seq) {
        this->seq = seq;
    }
    std::string getName() {
        return name;
    }
    std::string getSequence() {
        return seq;
    }
};

class Align {
private:
    double gap_open_penalty = 2.0;
    double gap_extend_penalty = 0.4;
    double end_gap_open_penalty = 0.2;
    double end_gap_extend_penalty = 0.1;
    double mismatch_penalty = -1.0;

public:
    Align() {}
    Align(double gap_open_penalty, double gap_extend_penalty, double end_gap_open_penalty,
          double end_gap_extend_penalty, double mismatch_penalty) :
        gap_open_penalty(gap_open_penalty), gap_extend_penalty(gap_extend_penalty),
        end_gap_open_penalty(end_gap_open_penalty), end_gap_extend_penalty(end_gap_extend_penalty),
        mismatch_penalty(mismatch_penalty) {}
    std::vector<std::string> align(Fasta& a, Fasta& b);
    double calculateGapPenalty(int position, int length, bool isOpen, bool isEnd) {
        if (isEnd) {
            return isOpen ? end_gap_open_penalty : end_gap_extend_penalty;
        } else {
            return isOpen ? gap_open_penalty : gap_extend_penalty;
        }
    }
    double calculateGapExtensionPenalty(int gapLength) {
        return gap_extend_penalty * log(1 + gapLength);
    }
    double get_mismatch_penalty() {
        return mismatch_penalty;
    }
};

struct ThreadData {
    int i;
    Fasta* a;
    Fasta* b;
    std::vector<std::vector<double>>* score;
    std::vector<std::vector<double>>* gap_a;
    std::vector<std::vector<double>>* gap_b;
    std::vector<std::vector<std::vector<int>>>* traceback;
    std::vector<std::vector<int>>* existingGapLengthA;
    std::vector<std::vector<int>>* existingGapLengthB;
    Align* alignInstance;
};

static void* initializeMatrix(void* arg) {
    ThreadData* data = static_cast<ThreadData*>(arg);
    int i = data->i;

    for (int j = 0; j <= data->b->seq.length(); j++) {
        (*data->score)[i][j] = 0.0;
        (*data->gap_a)[i][j] = std::numeric_limits<double>::lowest();
        (*data->gap_b)[i][j] = std::numeric_limits<double>::lowest();
        (*data->traceback)[i][j][0] = 0;
        (*data->traceback)[i][j][1] = 0;
        (*data->existingGapLengthA)[i][j] = 0;
        (*data->existingGapLengthB)[i][j] = 0;
    }
    return nullptr;
}

std::vector<std::string> Align::align(Fasta& a, Fasta& b) {
    // Declare matrices here, e.g.
    std::vector<std::vector<double>> score(a.seq.length() + 1, std::vector<double>(b.seq.length() + 1));
    std::vector<std::vector<double>> gap_a(a.seq.length() + 1, std::vector<double>(b.seq.length() + 1, std::numeric_limits<double>::lowest()));
    std::vector<std::vector<double>> gap_b(a.seq.length() + 1, std::vector<double>(b.seq.length() + 1, std::numeric_limits<double>::lowest()));
    std::vector<std::vector<std::vector<int>>> traceback(a.seq.length() + 1, std::vector<std::vector<int>>(b.seq.length() + 1, std::vector<int>(2)));
    std::vector<std::vector<int>> existingGapLengthA(a.seq.length() + 1, std::vector<int>(b.seq.length() + 1));
    std::vector<std::vector<int>> existingGapLengthB(a.seq.length() + 1, std::vector<int>(b.seq.length() + 1));

    // For each thread creation:
    // 1. Create a ThreadData instance and set its members
    // 2. Create the thread and pass the ThreadData instance to the thread function
    // 3. Join the threads after they finish execution
    // Score computation
    std::vector<pthread_t> init_threads(a.seq.length() + 1);
    for (int cur_i = 0; cur_i <= a.seq.length(); cur_i++) {
        ThreadData* data = new ThreadData();
        data->i = cur_i;
        data->a = &a;
        data->b = &b;
        data->score = &score;
        data->gap_a = &gap_a;
        data->gap_b = &gap_b;
        data->traceback = &traceback;
        data->existingGapLengthA = &existingGapLengthA;
        data->existingGapLengthB = &existingGapLengthB;
        data->alignInstance = this;

        pthread_t thread;
        pthread_create(&thread, NULL, initializeMatrix, data);
        init_threads[cur_i] = thread;
    }

    // Joining the initialization threads
    for (pthread_t thread : init_threads) {
        pthread_join(thread, NULL);
    }

    std::vector<pthread_t> score_threads(a.seq.length());
    for (int currentI = 1; currentI <= a.seq.length(); currentI++) {
        ThreadData* data = new ThreadData();
        data->i = currentI;
        data->a = &a;
        data->b = &b;
        data->score = &score;
        data->gap_a = &gap_a;
        data->gap_b = &gap_b;
        data->traceback = &traceback;
        data->existingGapLengthA = &existingGapLengthA;
        data->existingGapLengthB = &existingGapLengthB;
        data->alignInstance = this;

        pthread_t thread;
        pthread_create(&thread, NULL, computeScore, data);
        score_threads[currentI - 1] = thread;  // Adjust the index since the loop starts from 1
    }

    // Joining the score computation threads
    for (pthread_t thread : score_threads) {
        pthread_join(thread, NULL);
    }
    std::string alignedSeqA;
    std::string alignedSeqB;

    int i = a.seq.length();
    int j = b.seq.length();

    while (i > 0 && j > 0) {
        if (traceback[i][j][0] == 0) { // match or mismatch
            alignedSeqA = a.seq[i - 1] + alignedSeqA;
            alignedSeqB = b.seq[j - 1] + alignedSeqB;
            i--;
            j--;
        } else if (traceback[i][j][0] == 1) { // gap in sequence b
            alignedSeqA = a.seq[i - 1] + alignedSeqA;
            alignedSeqB = '-' + alignedSeqB;
            i--;
        } else { // gap in sequence a
            alignedSeqA = '-' + alignedSeqA;
            alignedSeqB = b.seq[j - 1] + alignedSeqB;
            j--;
        }
    }

    while (i > 0) {
        alignedSeqA = a.seq[i - 1] + alignedSeqA;
        alignedSeqB = '-' + alignedSeqB;
        i--;
    }

    while (j > 0) {
        alignedSeqA = '-' + alignedSeqA;
        alignedSeqB = b.seq[j - 1] + alignedSeqB;
        j--;
    }

    return {alignedSeqA, alignedSeqB};
}

static void* computeScore(void* arg) {
    ThreadData* data = static_cast<ThreadData*>(arg);
    int i = data->i;
    Fasta* a = data->a;
    Fasta* b = data->b;

    for (int j = 1; j <= b->seq.length(); j++) {
        (*data->existingGapLengthA)[i][j] = ((*data->traceback)[i - 1][j][0] == 1) ? (*data->existingGapLengthA)[i - 1][j] + 1 : 0;
        (*data->existingGapLengthB)[i][j] = ((*data->traceback)[i][j - 1][0] == 2) ? (*data->existingGapLengthB)[i][j - 1] + 1 : 0;

        bool isEndA = (i == 1 || i == a->seq.length());
        bool isEndB = (j == 1 || j == b->seq.length());

        double penalty_gap_open_a = data->alignInstance->calculateGapPenalty(i, a->seq.length(), true, isEndA);
        double penalty_gap_open_b = data->alignInstance->calculateGapPenalty(j, b->seq.length(), true, isEndB);
        double penalty_gap_extend_a = data->alignInstance->calculateGapPenalty(i, a->seq.length(), false, isEndA) + data->alignInstance->calculateGapExtensionPenalty((*data->existingGapLengthA)[i][j]);
        double penalty_gap_extend_b = data->alignInstance->calculateGapPenalty(j, b->seq.length(), false, isEndB) + data->alignInstance->calculateGapExtensionPenalty((*data->existingGapLengthB)[i][j]);

        (*data->gap_a)[i][j] = std::max((*data->score)[i - 1][j] - penalty_gap_open_a, (*data->gap_a)[i - 1][j] - penalty_gap_extend_a);
        (*data->gap_b)[i][j] = std::max((*data->score)[i][j - 1] - penalty_gap_open_b, (*data->gap_b)[i][j - 1] - penalty_gap_extend_b);

        double match = (*data->score)[i - 1][j - 1] + (a->seq[i - 1] == b->seq[j - 1] ? 1.0 : data->alignInstance->get_mismatch_penalty());

        double max_score = match;
        (*data->traceback)[i][j][0] = 0;
        (*data->traceback)[i][j][1] = 0;

        if (max_score < (*data->gap_a)[i][j]) {
            max_score = (*data->gap_a)[i][j];
            (*data->traceback)[i][j][0] = 1;
        }
        if (max_score < (*data->gap_b)[i][j]) {
            max_score = (*data->gap_b)[i][j];
            (*data->traceback)[i][j][0] = 2;
        }

        (*data->score)[i][j] = max_score;
    }

    return nullptr;
}

// Function to load fasta sequences from a gzipped file
std::unordered_map<std::string, std::string> loadFasta(const std::string& path) {
    gzFile file = gzopen(path.c_str(), "r");
    if (!file) {
        std::cerr << "Cannot open file: " << path << std::endl;
        exit(1);
    }

    char buffer[1024];
    std::string line;
    std::string name;
    std::unordered_map<std::string, std::string> map;
    std::string sequence;
    std::cerr << "Loading " << path << std::endl;

    while (gzgets(file, buffer, sizeof(buffer))) {
        line = buffer;
        // ここで行末の改行文字を削除します
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());

        if (line[0] == '>') {
            if (!name.empty()) {
                map[name] = sequence;
                cerr << "loaded: " << name << endl;
            }
            std::istringstream iss(line);
            iss >> name;
            name = name.substr(1); // Remove '>'
            sequence.clear();
        } else {
            sequence += line;
        }
    }
    if (!name.empty()) {
        map[name] = sequence;
    }
    gzclose(file);
    return map;
}


int main(int argc, char** argv) {
    try {
        std::unordered_map<std::string, std::string> human_mfasta = loadFasta(argv[1]);
        std::unordered_map<std::string, std::string> mouse_mfasta = loadFasta(argv[2]);
        Align align;

        std::ifstream file(argv[3]);
        std::string raw;
        Fasta human;
        Fasta mouse;

        while (std::getline(file, raw)) {
            std::istringstream iss(raw);
            std::vector<std::string> line;
            std::string temp;
            while (iss >> temp) {
                line.push_back(temp);
            }

            if ( line.size() < 2) {
                continue;
            }

            std::string mouseID = line[0];
            std::string humanID = line[1];
            std::string humanSeq = human_mfasta[humanID];
            std::string mouseSeq = mouse_mfasta[mouseID];

            if (humanSeq.empty() || mouseSeq.empty()) {
                std::cerr << "Skip " << raw << std::endl;
                std::cerr << "Reason: " << (humanSeq.empty() ? humanID : mouseID) << std::endl;
                continue;
            }

            std::cerr << "Processing " << line[0] << " " << line[1] << std::endl;
            human.setSequence(humanSeq);
            mouse.setSequence(mouseSeq);
            human.setName(humanID);
            mouse.setName(mouseID);

            std::vector<std::string> result = align.align(human, mouse);
            std::string match;
            // check uninitialized result
            if (result.size() != 2) {
                std::cerr << "Skip " << raw << std::endl;
                std::cerr << "Reason: " << "Uninitialized result" << std::endl;
                continue;
            }
            for (size_t i = 0; i < result[0].length(); i++) {
                match += (result[0][i] == result[1][i]) ? '|' : ' ';
            }

            std::cout << human.getName() << "\t" << result[0] << std::endl;
            std::cout << "                \t" << match << std::endl;
            std::cout << mouse.getName() << "\t" << result[1] << std::endl;
            std::cout << "----" << std::endl;
        }
        file.close();
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
