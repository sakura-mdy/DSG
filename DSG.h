
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>
#include <map>
#include <string>
#include <tuple>

using namespace std;

extern vector<pair<int, int>> shell_new_edges;
extern vector<pair<int, int>> new_edges;

class EfficientCMAlgorithm {
public:
    EfficientCMAlgorithm() {
        cannot_insert = false;
    }
    void output_data(const string &out_file, const double &total_t, const string &algorithm, int &edges, int &num_followers, const string &graph_name);
    tuple<int, int, vector<pair<int, int>>> greedy_solution_selection(int b);
    void partition_shell();
    void partition_shell_FC(const string &graph_name);
    void complete_conversion();
    tuple<int, int, vector<pair<int, int>>> complete_conversion_FC(int b, const string &graph_name);
    tuple<int, int, vector<pair<int, int>>> EfficientCM(int b);
    tuple<int, int, vector<pair<int, int>>> EfficientCM_FG(int b, const string &graph_name);
    bool cannot_insert;
private:
    map<int, int> is_collapse;
    unordered_map<int, vector<vector<pair<int, int>>>> solution_condidates; 
    unordered_map<int, vector<int>> followers_number; 
    unordered_map<int, vector<int>> biedge_number; 
    unordered_map<int, vector<vector<int>>> followers_items; 
    unordered_map<int, vector<vector<int>>> biedge_nodes;
    map<int, unordered_map<int, int>> k_degree; 
    unordered_map<int, int> degree_FC;
    vector<int> single_component;
    unordered_map<int, vector<vector<int>>> cc; 
    unordered_map<int, pair<int, int>> node_partition;
    unordered_map<int, vector<double>> cc_num;
};