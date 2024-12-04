
#ifndef _GRAPH_H_
#define _GRAPH_H_
#include <ctime>  
#include <iostream>
#include <string>
#include <memory>
#include <cstdlib> 
#include <cstdio>  
#include <cmath>   
#include <algorithm>
#include <map>
#include <unordered_map>
#include <vector>

using namespace std;

// graph 
struct endpoint
{
    int u;
    int eid;
};

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

extern unordered_map<pair<int, int>, int, pair_hash> eset;
extern map<int, int> coreness;
extern unordered_map<int, vector<endpoint>> g;
extern vector<int> k_core_vertices;
extern int n, m, dmax;

extern void Readin_Graph(const char *str);
extern void core_decompostion();

#endif  /* _GRAPH_H_ */
