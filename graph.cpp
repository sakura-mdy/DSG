#include "graph.h"
#include <iostream>
#include <vector>
#include <memory>
#include <unordered_map>
#include <cstdio>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <psapi.h>
#else
#include <unistd.h>
#include <sys/resource.h>
#endif

using namespace std;

extern int K, b, depth;
map<int, int> coreness;
vector<int> k_core_vertices;
unordered_map<int, vector<endpoint>> g;
unordered_map<pair<int, int>, int, pair_hash> eset;
int n, m, dmax, M;

void PrintMemoryUsage() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS memInfo;
    GetProcessMemoryInfo(GetCurrentProcess(), &memInfo, sizeof(memInfo));
    // cout << "Memory usage: " << memInfo.WorkingSetSize / 1024 << " kilobytes" << endl;
#else
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    // cout << "Memory usage: " << usage.ru_maxrss << " kilobytes" << endl;
#endif
}

void ProcessEdges(const vector<pair<int, int>>& edges, size_t& total_processed_edges) {
    endpoint a, b;
    int eid = eset.size();  // Start with the current size of eset

    for (const auto& edge : edges) {
        int x = edge.first;
        int y = edge.second;

        if (eset.find({x, y}) == eset.end()) {
            eset[{x, y}] = eid;
            eset[{y, x}] = eid;
            a.u = y;
            a.eid = eid;
            g[x].push_back(a);
            b.u = x;
            b.eid = eid;
            g[y].push_back(b);
            eid++;
            m++;
        }
    }

    for (const auto& p : g) {
        dmax = max(dmax, static_cast<int>(p.second.size()));
    }

    total_processed_edges += edges.size();
    // cout << "Processed " << total_processed_edges << " edges so far." << endl;
    PrintMemoryUsage();
}

void Readin_Graph(const char *str) {
    FILE *fin = fopen(str, "r");
    if (fin == nullptr) {
        perror("Fail to open file");
        exit(EXIT_FAILURE);
    }

    int x, y;
    n = 0;
    m = 0;
    dmax = 0;

    vector<pair<int, int>> edges;
    const size_t block_size = 100000;  // Define the block size
    size_t count = 0;
    size_t total_processed_edges = 0;

    char line[256];
    while (fgets(line, sizeof(line), fin)) {
        if (line[0] == '#') {
            continue; // Skip comment lines
        }
        if (sscanf(line, "%d %d", &x, &y) != 2) {
            continue; // Skip malformed lines
        }
        if (x != y) {
            if (y < x) swap(x, y);
            edges.emplace_back(x, y);
            n = max(n, max(x, y));
            count++;
            if (count >= block_size) {
                ProcessEdges(edges, total_processed_edges);
                edges.clear();
                count = 0;
            }
        }
    }

    if (!edges.empty()) {
        ProcessEdges(edges, total_processed_edges);
    }

    n++;
    fclose(fin);

    cout << "read over" << endl;
    printf("n = %d, m = %d, dmax = %d\n", n, m, dmax);
}

void core_decompostion() {
    unique_ptr<vector<int>[]> L(new vector<int>[dmax + 1]);
    unique_ptr<int[]> degg(new int[n + 1]);
    unique_ptr<bool[]> del(new bool[n + 1]);
    int maxcore = 0;
    unordered_map<int, vector<endpoint>>::iterator p;

    for (p = g.begin(); p != g.end(); p++) {
        int x = p->first;
        degg[x] = p->second.size();
        L[degg[x]].push_back(x);
    }

    fill(del.get(), del.get() + n + 1, false);

    int md = dmax;
    for (int k = 0; k <= dmax; k++) {
        for (size_t i = 0; i < L[k].size(); i++) {
            int x = L[k][i];
            if (del[x]) {
                continue;
            }
            coreness[x] = degg[x];
            if (maxcore < coreness[x]) {
                maxcore = coreness[x];
            }
            del[x] = true;
            for (size_t j = 0; j < g[x].size(); ++j) {
                int y = g[x][j].u;
                if (del[y]) {
                    continue;
                }
                int temp = y;
                if (degg[temp] > k) {
                    degg[temp]--;
                    L[degg[temp]].push_back(temp);
                }
            }
        }
    }

    printf("maxdegree=%d, maxcore=%d \n", md, maxcore);
    for (p = g.begin(); p != g.end(); p++) {
        int i = p->first;
        if (coreness[i] >= K) {
            k_core_vertices.push_back(i);
        }
    }

    int sum = 0;
    for (int i = 1; i <= K; ++i) {
        int count = 0;
        for (const auto& node : coreness) {
            if (node.second == (K - i)) {
                count++;
            }
        }
        sum += count * i;
        if (sum > 2 * b) {
            depth = i + 2;
            break;
        }
    }
}
