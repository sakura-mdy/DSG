#include "graph.h"
#include "DSG.h"
#include <iostream>
#include <chrono>
#include <memory>
#include <fstream>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif
using namespace std;

void create_directory(const string& path) {
#ifdef _WIN32
    _mkdir(path.c_str());
#else
    mkdir(path.c_str(), 0777);
#endif
}

string extract_filename(const string& path) {
    size_t pos = path.find_last_of("/\\");
    if (pos != string::npos) {
        return path.substr(pos + 1);
    }
    return path;
}

int K, b, depth;
int lambda = 1;

int main(int argc, char* argv[]){
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <graph_name> <K> <b>" << endl;
        return 1;
    }

    string graph_name = argv[1];
    string graph_file = extract_filename(graph_name);
    K = atoi(argv[2]);
    b = atoi(argv[3]);
    string out_file = "results.txt";
    string output_folder = "output";

    // Create output folder if it does not exist
    create_directory(output_folder);

    Readin_Graph(graph_name.c_str());
    unique_ptr<EfficientCMAlgorithm> algo = make_unique<EfficientCMAlgorithm>();
    core_decompostion();

    // *****  DSG CG procedure:
    auto start = std::chrono::high_resolution_clock::now();
    tuple<int, int, vector<pair<int, int>>> res = algo->EfficientCM(b);
    auto end = std::chrono::high_resolution_clock::now();
    double total_t = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
    int num_edges = get<0>(res), num_nodes = get<1>(res);
    vector<pair<int, int>> edges = get<2>(res);
    int count = 0;
    while(num_edges <= 0.98*b){
        depth++;
        count++;
        res = algo->EfficientCM(b);
        num_nodes = get<1>(res), num_edges = get<0>(res);
        vector<pair<int, int>> edges = get<2>(res);
    }
    string edge_file_name = output_folder + "/added_edges_" + graph_file;
    ofstream edge_file(edge_file_name);
    if (!edge_file) {
        cerr << "Error opening file: " << edge_file_name << endl;
        return 1;
    }
    for (const auto &edge : edges) {
        edge_file << edge.first << " " << edge.second << endl;
    }
    edge_file.close();
    algo->output_data(out_file, total_t, string("DSG CG"), num_edges, num_nodes, graph_name);

    // *****  DSG PG procedure:
    depth -= count;
    auto start_FG = std::chrono::high_resolution_clock::now();
    tuple<int, int, vector<pair<int, int>>> res_FG = algo->EfficientCM_FG(b, graph_file);
    auto end_FG = std::chrono::high_resolution_clock::now();
    double total_t_FG = double(std::chrono::duration_cast<std::chrono::milliseconds>(end_FG - start_FG).count()) / 1000;
    int num_nodes_FG = get<1>(res_FG), num_edges_FG = get<0>(res_FG);
    vector<pair<int, int>> edges_FG = get<2>(res_FG);
    // while(num_edges_FG <= 0.98*b){
    //     depth++;
    //     res_FG = algo->EfficientCM_FG(b, graph_file);
    //     num_nodes_FG = get<1>(res_FG), num_edges_FG = get<0>(res_FG);
    //     vector<pair<int, int>> edges_FG = get<2>(res_FG);
    // }
    string edge_file_name_FG = output_folder + "/added_edges_FG_" + graph_file;
    ofstream edge_file_FG(edge_file_name_FG);
    if (!edge_file_FG) {
        cerr << "Error opening file: " << edge_file_name_FG << endl;
        return 1;
    }
    for (const auto &edge : edges_FG) {
        edge_file_FG << edge.first << " " << edge.second << endl;
    }
    edge_file_FG.close();
    algo->output_data(out_file, total_t_FG, string("FG-EfficientCM"), num_edges_FG, num_nodes_FG, graph_name);

    return 0;
}