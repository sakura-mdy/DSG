#include "DSG.h"
#include <chrono>
#include "graph.h"
#include <iostream>
#include <numeric>
#include <stack>
#include <cstdlib> // for rand()
#include <optional>
#include <utility>
#include <chrono>

using namespace std::chrono;
using namespace std;

extern int K, b, lambda, depth;
extern map<int, int> coreness;
extern vector<int> k_core_vertices;

pair<int, int> normalize_edge(int u, int v) {
    return (u < v) ? make_pair(u, v) : make_pair(v, u);
}

int customGreedyAlgorithm(std::vector<int> nums, int additionalValue) {
	
    nums.push_back(additionalValue);
    int initialSum = std::accumulate(nums.begin(), nums.end(), 0);
	std::sort(nums.begin(), nums.end(), std::greater<int>());
    
	while (nums.size() > 1) {
        int max1 = nums[0];
        int max2 = nums[1];
        nums.erase(nums.begin());
        nums.erase(nums.begin());

        int diff = max1 - max2;
        if (diff > 0) {
            nums.insert(nums.begin(), diff); 
        }
        if (diff < 0) {
            nums.insert(nums.begin(), -diff); 
        }
    }

    int finalResult = nums.empty() ? 0 : nums[0];
    int result = finalResult + (initialSum - finalResult) / 2;
    return result;
}

void dfs(int u, int core, vector<int>& component, unordered_map<int, bool>& visited) {
    stack<int> s;
    s.push(u);
    visited[u] = true;

    while (!s.empty()) {
        int node = s.top();
        s.pop();
        component.push_back(node);

        for (const auto& neighbor : g[node]) {
            int v = neighbor.u;
            if (!visited[v] && coreness[v] == core) {
                s.push(v);
                visited[v] = true;
            }
        }
    }
}

void EfficientCMAlgorithm::output_data(const string &out_file, const double &total_t, const string &algorithm,int &edges, int &num_followers, const string &graph_name) {
	std::fstream ff(out_file, std::ios::out|std::ios::app);
	ff << "Graph: " <<  graph_name << endl;
	ff << algorithm << ": k = " << K << " b:" << b << endl;
	ff << " the number of followers:" << num_followers << endl;
	ff << "The new edges identified by:" << algorithm << endl;
	ff << "depth is:" << depth << endl;
	ff << "The new edges identified by:" << algorithm << endl;
	ff << edges << endl;
	ff.close();
}

void EfficientCMAlgorithm::complete_conversion() {
    for (auto &cc_layer : cc) {
        int flag = cc_layer.first;
        int total_edges_added = 0; // 记录当前 flag 添加的边的总数量

        for (const auto& component : cc_layer.second) {

            unordered_set<int> nodes(component.begin(), component.end());

            // 当前连通部分的节点数和节点
            followers_number[flag].push_back(component.size());
            followers_items[flag].push_back(component);

            vector<pair<int, int>> new_edges;
            vector<int> new_biedge_nodes;
            int res = 0, m = 0;
            for (int u : component) {
			    // 计算节点 u 需要添加的边数，使其度数达到 K
			    int deg_ = max(0, K - k_degree[flag][u]);
			    // 累加所有节点的 deg_ 值
			    res += deg_;
			}
			if(res < 2*b ){
	            for (int u : component) {
	                while (k_degree[flag][u] < K) {
	                    // 查找一个在连通部分内且 k_degree 小于 K 的节点进行连接
	                    bool added = false;
	                    for (int v : component) {
	                        if (u != v && k_degree[flag][v] < K && eset.find(normalize_edge(u, v)) == eset.end() && find(new_edges.begin(), new_edges.end(), normalize_edge(u, v)) == new_edges.end()) {
	                            new_edges.emplace_back(normalize_edge(u, v));
	                            k_degree[flag][u]++;
	                            k_degree[flag][v]++;
	                            added = true;
	                            total_edges_added++; // 更新添加的边的总数量
	                            break;
	                        }
	                    }
	                    if (!added && k_degree[flag][u] < K-1) {
	                        int random_k_core_vertex;
	                        do {
	                            random_k_core_vertex = k_core_vertices[rand() % k_core_vertices.size()];
	                        } while (eset.find(normalize_edge(u, random_k_core_vertex)) != eset.end() || find(new_edges.begin(), new_edges.end(), normalize_edge(u, random_k_core_vertex)) != new_edges.end());
	
	                        new_edges.emplace_back(normalize_edge(u, random_k_core_vertex));
	                        k_degree[flag][u]++;
	                        total_edges_added++; // 更新添加的边的总数量
	                    }
	                    if (!added && k_degree[flag][u] == K-1) {
	                    	new_biedge_nodes.emplace_back(u);
	                    	m++;
	                        break;
	                    }	
	                }
	            }				
			}
			biedge_number[flag].push_back(m);
			biedge_nodes[flag].push_back(new_biedge_nodes);
            // 记录新增的边
            solution_condidates[flag].push_back(new_edges);
        }
    }
}

tuple<int, int, vector<pair<int, int>>> EfficientCMAlgorithm::complete_conversion_FC(int b, const string& graph_name) {
    unordered_map<int, vector<double>> cc_eff;

    // 计算cc_eff
    for (auto cc_it = cc.begin(); cc_it != cc.end(); ++cc_it) {
        int shell_index = cc_it->first;
        vector<vector<int>> components = cc_it->second;
        for (int i = 0; i < components.size(); ++i) {
            double eff = components[i].size() / cc_num[shell_index][i];
            cc_eff[shell_index].push_back(eff);
        }
    }

    vector<pair<int, int>> selected_components;
    vector<int> final_solution;
    int total_sum = 0;
    unordered_map<int, vector<bool>> visited_eff;
    final_solution.reserve(2000);

    // 初始化visited_eff
    for (auto& pair : cc_eff) {
        visited_eff[pair.first] = std::vector<bool>(pair.second.size(), false);
    }

    int added_nodes_count = 0;
    unordered_map<int, vector<vector<pair<int, int>>>> external_connections;
    int f_num = 0;

    auto total_time_degree_update = 0; // 更新度数总时间
    auto total_time_degree_update_2 = 0; 
    auto total_time_sort = 0; // 排序节点总时间
    auto total_time_edge_calculation = 0; // 计算边数量总时间
    auto total_time_edge_addition = 0; // 添加边的总时间
    auto total_time_edge_addition_f = 0; // 添加边的总时间
    int coreness_each;

    // 检查文件是否存在，并读取其内容
    // auto start_time = high_resolution_clock::now();
    ifstream input_file("sorted_" + graph_name);
    bool use_file_data = false;

    if (input_file.is_open()) {
        // cout<<"offline"<<endl;
        int file_b;
        input_file >> file_b;

        if (file_b >= b) {
            cout<<"exist"<<endl;
            use_file_data = true;
            int max_shell_index,max_component_index;
            unordered_map<int, int> degree_FC_copy = degree_FC;

            // 从文件中读取max_shell_index和max_component_index
            while (input_file >> max_shell_index >> max_component_index) {
                // cout<<"start"<<endl;
                // cout<<"this line"<<endl;
                // cout<<"max_shell_index is:"<<max_shell_index<<endl;
                // cout<<"max_component_index is"<<max_component_index<<endl;

                // 1. 增加度数的逻辑
                // auto start_time_in = high_resolution_clock::now();
                if(max_shell_index >= cc.size()){
                    continue;
                }
                unordered_set<int> final_solution_set(final_solution.begin(), final_solution.end());  // 使用unordered_set加速查找
                for (int u : cc[max_shell_index][max_component_index]) {
                    // cout<<"ready"<<endl;
                    int coreness_u = coreness[u];  // 提前存储coreness[u]的值
                    // coreness_each = coreness[u];
                    for (const auto& neighbor : g[u]) {
                        int v = neighbor.u;
                        // 使用unordered_set查找代替vector查找，提高效率
                        if (final_solution_set.find(v) != final_solution_set.end() && coreness[v] != coreness_u) {
                            degree_FC[u]++;
                            degree_FC[v]++;
                        }
                    }
                }
                // cout<<"one end"<<endl;
                // cout<<"this time's coreness is"<<coreness_each<<endl;
                // unordered_map<int, int> degree_FC_copy_2 = degree_FC_copy;
                // auto end_time_in = high_resolution_clock::now();
                // total_time_degree_update += duration_cast<milliseconds>(end_time_in - start_time_in).count();

                // 2. 检查加入解集所需的边的数量
                // auto start_time_in0 = high_resolution_clock::now();
                vector<int> nodes_to_connect;
                int num=0;
                for (int v : final_solution) {
                    if (degree_FC[v] < K) {
                        nodes_to_connect.push_back(v);
                        num+=(K-degree_FC[v]);
                    }
                }
                for (int v : cc[max_shell_index][max_component_index]) {
                    if (degree_FC[v] < K) {
                        nodes_to_connect.push_back(v);
                        num+=(K-degree_FC[v]);
                    }
                }
                // auto end_time_in0 = high_resolution_clock::now();
                // total_time_degree_update_2 += duration_cast<milliseconds>(end_time_in0 - start_time_in0).count();
                if(num>2*b){
                    // cout<<"no"<<endl;
;                    f_num++;
                    for (int u : cc[max_shell_index][max_component_index]) {
                        int coreness_u = coreness[u];  // 提前存储coreness[u]的值
                        for (const auto& neighbor : g[u]) {
                            int v = neighbor.u;
                            // 使用unordered_set查找代替vector查找，提高效率
                            if (final_solution_set.find(v) != final_solution_set.end() && coreness[v] != coreness_u) {
                                degree_FC[u]--;
                                degree_FC[v]--;
                            }
                        }
                    }                  
                    if(f_num==10){
                        // cout<<"end"<<endl;
                        break;
                    }
                    continue;
                }else{
                    f_num=0;
                }

                // // 3. 排序节点
                // auto start_time_in2 = high_resolution_clock::now();
                // std::sort(nodes_to_connect.begin(), nodes_to_connect.end(), [this](int a, int b) {
                //     return degree_FC[a] < degree_FC[b];
                // });
                // auto end_time_in2 = high_resolution_clock::now();
                // total_time_sort += duration_cast<milliseconds>(end_time_in2 - start_time_in2).count();

                // // 4. 计算需要添加的边的数量
                // auto start_time_in3 = high_resolution_clock::now();
                // vector<pair<int, int>> edges_to_add;
                // for (size_t i = 0; i < nodes_to_connect.size(); ++i) {
                //     for (size_t j = i + 1; j < nodes_to_connect.size(); ++j) {
                //         int u = nodes_to_connect[i];
                //         int v = nodes_to_connect[j];
                //         if (degree_FC_copy_2[u] >= K) break;
                //         if (degree_FC_copy_2[u] < K && degree_FC_copy_2[v] < K) {
                //             if (eset.find(normalize_edge(u, v)) == eset.end() && find(edges_to_add.begin(), edges_to_add.end(), normalize_edge(u, v)) == edges_to_add.end()) {
                //                 edges_to_add.push_back(normalize_edge(u, v));
                //                 degree_FC_copy_2[u]++;
                //                 degree_FC_copy_2[v]++;
                //             }
                //         }
                //     }
                // }
                // auto end_time_in3 = high_resolution_clock::now();
                // total_time_edge_calculation += duration_cast<milliseconds>(end_time_in3 - start_time_in3).count();

                // // 5. 补充添加边的逻辑，确保所有节点的度数达到K
                // auto start_time_in4 = high_resolution_clock::now();
                // for (int u : final_solution) {
                //     while (degree_FC_copy_2[u] < K) {
                //         int v = k_core_vertices[rand() % k_core_vertices.size()];
                //         if (u != v && eset.find(normalize_edge(u, v)) == eset.end() && find(edges_to_add.begin(), edges_to_add.end(), normalize_edge(u, v)) == edges_to_add.end()) {
                //             edges_to_add.push_back(normalize_edge(u, v));
                //             degree_FC_copy_2[u]++;
                //             degree_FC_copy_2[v]++;
                //         }
                //     }
                // }
                // auto end_time_in4 = high_resolution_clock::now();
                // total_time_edge_addition += duration_cast<milliseconds>(end_time_in4 - start_time_in4).count();

                // // 如果需要添加的边数量大于b，跳过该部分
                // if (edges_to_add.size() > b) {
                //     continue;
                // }

                // 将该部分加入解集中
                // cout<<"append"<<endl;
                // auto start_time_in_f = high_resolution_clock::now();
                total_sum += cc_num[max_shell_index][max_component_index];
                selected_components.push_back({max_shell_index, max_component_index});
                final_solution.insert(final_solution.end(), cc[max_shell_index][max_component_index].begin(), cc[max_shell_index][max_component_index].end());
                added_nodes_count += cc[max_shell_index][max_component_index].size();
                // auto end_time_in_f = high_resolution_clock::now();
                // total_time_edge_addition_f += duration_cast<milliseconds>(end_time_in_f - start_time_in_f).count();
            }
        }
        input_file.close();
    }
    // auto end_time = high_resolution_clock::now();
    // auto duration_file_process = duration_cast<milliseconds>(end_time - start_time).count();
    // cout << "duration_file_process: " << duration_file_process << " ms" << endl;

    // 输出每个部分的总时间
    // cout << "total_time_degree_update: " << total_time_degree_update << "ms" << endl;
    // cout << "total_time_degree_update_2: " << total_time_degree_update_2 << "ms" << endl;
    // cout << "total_time_sort: " << total_time_sort << "ms" << endl;
    // cout << "total_time_edge_calculation: " << total_time_edge_calculation << "ms" << endl;
    // cout << "total_time_edge_addition: " << total_time_edge_addition << "ms" << endl;
    // cout << "total_time_edge_addition_f: " << total_time_edge_addition_f << "ms" << endl;
    if (!use_file_data) {
        while (true) {
            double max_eff = -1;
            int max_shell_index = -1;
            int max_component_index = -1;
            for (auto& pair : cc_eff) {
                int shell_index = pair.first;
                vector<double> effs = pair.second;
                for (int i = 0; i < effs.size(); ++i) {
                    if (!visited_eff[shell_index][i] && effs[i] < 0 ) {
                        cout<<"出现异常值"<< effs[i] <<endl;
                        cin.get();
                    }
                }
            }
            for (auto& pair : cc_eff) {
                int shell_index = pair.first;
                vector<double> effs = pair.second;
                for (int i = 0; i < effs.size(); ++i) {
                    if (!visited_eff[shell_index][i]) {
                        if (effs[i] > max_eff) {
                            max_eff = effs[i];
                            max_shell_index = shell_index;
                            max_component_index = i;
                        } else if (effs[i] == max_eff && cc[shell_index][i].size() > cc[max_shell_index][max_component_index].size()) {
                            max_eff = effs[i];
                            max_shell_index = shell_index;
                            max_component_index = i;
                        }
                    }
                }
            }
            if (max_shell_index == -1 || max_component_index == -1) {
                break;
            }

            // 更新访问状态
            visited_eff[max_shell_index][max_component_index] = true;

            // 创建degree_FC的副本
            unordered_map<int, int> degree_FC_copy = degree_FC, degree_FC_copy_2 = degree_FC;

            for (int u : cc[max_shell_index][max_component_index]) {
                for (const auto& neighbor : g[u]) {
                    int v = neighbor.u;
                    if (find(final_solution.begin(), final_solution.end(), v) != final_solution.end() && coreness[v] != coreness[u]) {
                        degree_FC_copy_2[u]++;
                        degree_FC_copy_2[v]++;
                    }
                }
            }
            
            double degree_sum_FC = 0;
            for (int v : final_solution) {
                int degree = degree_FC_copy_2[v] > K ? 0 : K - degree_FC_copy_2[v];
                degree_sum_FC += degree/2.0;
            }   
            for (int v : cc[max_shell_index][max_component_index]) {
                int degree = degree_FC_copy_2[v] > K ? 0 : K - degree_FC_copy_2[v];
                degree_sum_FC += degree/2.0;
            }   		     

            // 计算cc_num加上后的值是否大于b
            if (degree_sum_FC > b) {
                f_num++;
                if(f_num == 10){
                    break;
                }
                continue; 
            }else{
                f_num = 0;
            }
            // 否则，执行某些操作
            if(degree_sum_FC > 0.99*b){
                vector<int> nodes_to_connect;
                for (int v : final_solution) {
                    if (degree_FC_copy_2[v] < K) {
                        nodes_to_connect.push_back(v);
                    }
                }
                for (int v : cc[max_shell_index][max_component_index]) {
                    if (degree_FC_copy_2[v] < K) {
                        nodes_to_connect.push_back(v);
                    }
                }
                
                // 对度数不到K的节点进行升序排列
                std::sort(nodes_to_connect.begin(), nodes_to_connect.end(), [&degree_FC_copy_2](int a, int b) {
                    return degree_FC_copy_2[a] < degree_FC_copy_2[b];
                });
        
                // 连接节点，确保其度数达到K
                vector<pair<int, int>> edges_to_add;
                for (size_t i = 0; i < nodes_to_connect.size(); ++i) {
                    for (size_t j = i + 1; j < nodes_to_connect.size(); ++j) {
                        int u = nodes_to_connect[i];
                        int v = nodes_to_connect[j];
                        if (degree_FC_copy_2[u] >= K) break;
                        if (degree_FC_copy_2[u] < K && degree_FC_copy_2[v] < K) {
                            // 检查节点之间是否有边
                            if (eset.find(normalize_edge(u, v)) == eset.end() && find(edges_to_add.begin(), edges_to_add.end(), normalize_edge(u, v)) == edges_to_add.end()) {
                                edges_to_add.push_back(normalize_edge(u, v));
                                degree_FC_copy_2[u]++;
                                degree_FC_copy_2[v]++;
                            }
                        }
                    }
                }
        
                // 如果节点之间无法连满，随意连接到k_core_vertices中的节点
                for (int u : nodes_to_connect) {
                    while (degree_FC_copy_2[u] < K) {
                        int v = k_core_vertices[rand() % k_core_vertices.size()]; // 从k_core_vertices中随机选择一个节点v
                        if (u != v && eset.find(normalize_edge(u, v)) == eset.end() && find(edges_to_add.begin(), edges_to_add.end(), normalize_edge(u, v)) == edges_to_add.end()) {
                            edges_to_add.push_back(normalize_edge(u, v));
                            degree_FC_copy_2[u]++;
                            degree_FC_copy_2[v]++;
                        }
                    }
                }
                // 检查添加的边是否大于b
                if (edges_to_add.size() > b) {
                    continue;
                }		
            }	
            // 记录已处理的节点
            unordered_set<int> processed_nodes;
            
            // 记录增加了度数的节点
            unordered_set<int> updated_nodes;
            vector<pair<int, int>> affected_indices;	//预先存储的所有有改动的索引
            vector<int> affected_node;
            
            // 第一步：记录所有增加了度数的节点
            for (int u : cc[max_shell_index][max_component_index]) {
                for (const auto& neighbor : g[u]) {
                    int v = neighbor.u;
                    // 新的判断条件：如果v在解集中且coreness不等于u
                    if (find(final_solution.begin(), final_solution.end(), v) != final_solution.end() && coreness[v] != coreness[u]) {
                        degree_FC[u]++;
                        degree_FC[v]++;
                        degree_FC_copy[u]++;
                        degree_FC_copy[v]++;
    //		            updated_nodes.insert(u);
                        if(coreness[v] < K){
                            updated_nodes.insert(v);						
                        }
                    }
                }
            }
            
            // 第二步：查找所有这些节点的邻居并更新相关数据
            for (int node : updated_nodes) {
                vector<pair<int, int>> process_core;
                for (const auto& neighbor : g[node]) {
                    int neighbor_node = neighbor.u;
                    if(coreness[neighbor_node] < K && coreness[neighbor_node] >= K-depth && find(final_solution.begin(), final_solution.end(), neighbor_node) == final_solution.end() && 
                    find(cc[max_shell_index][max_component_index].begin(), cc[max_shell_index][max_component_index].end(), neighbor_node) == cc[max_shell_index][max_component_index].end()){
                        // 存储node_partition的索引
                        int part_first = node_partition[neighbor_node].first;
                        int part_second = node_partition[neighbor_node].second;
                        pair<int, int> current_pair = make_pair(part_first, part_second);
                        if (find(process_core.begin(), process_core.end(), current_pair) != process_core.end()) {
                            continue;
                        }else{
                            process_core.push_back(current_pair);
                        }
                        int count = 0;	//neighbor_node节点所在区域与node节点的省边连接次数
                        for (const auto& p : external_connections[node_partition[neighbor_node].first][node_partition[neighbor_node].second]) {
                            if (p.first == node) {
                                count = p.second;
                                break;
                            }
                        }
                        int remaining_capacity = max(0, K - degree_FC_copy[node]);
                        if (count > remaining_capacity) {
                            int difference = count - remaining_capacity;
                            cc_num[part_first][part_second] += difference * 0.5;
                            pair<int, int> new_pair = {part_first, part_second};
                            // 将索引加入affected_indices
                            if (std::find(affected_indices.begin(), affected_indices.end(), new_pair) == affected_indices.end()) {
                                affected_indices.emplace_back(new_pair);
                            } 
        
                            // 修改external_connections
                            for (auto& p : external_connections[part_first][part_second]) {
                                if (p.first == node) {
                                    p.second = remaining_capacity;
                                    break;
                                }
                            }
                        }					
                    }
                }
            }

            // 保存每个独立区域的初始度数
            unordered_map<int, int> initial_degree_FC_copy = degree_FC_copy;

            // 将你提供的第二段代码合并进来
            for (int node : cc[max_shell_index][max_component_index]) {
                for (const auto& neighbor : g[node]) {
                    int neighbor_node = neighbor.u;
                    if (coreness[neighbor_node] < K && coreness[neighbor_node] >= K - depth && find(final_solution.begin(), final_solution.end(), neighbor_node) == final_solution.end() &&
                        coreness[neighbor_node] != coreness[node]) {
                        int part_first = node_partition[neighbor_node].first;
                        int part_second = node_partition[neighbor_node].second;
                        optional<pair<int, int>> pair_1 = nullopt;
                        optional<pair<int, int>> pair_2 = nullopt;

                        auto& vec = external_connections[part_first];
                        if (vec.size() <= static_cast<size_t>(part_second)) {
                            vec.resize(part_second + 1);
                        }

                        for (auto& p : external_connections[part_first][part_second]) {
                            if (p.first == node) {
                                pair_1 = p;
                            }
                            if (p.first == neighbor_node) {
                                pair_2 = p;
                            }
                            if (pair_1.has_value() && pair_2.has_value()) {
                                break;
                            }
                        }

                        if (!pair_1.has_value()) {
                            external_connections[part_first][part_second].emplace_back(node, 0);
                            pair_1 = external_connections[part_first][part_second].back();
                        }

                        if (!pair_2.has_value()) {
                            external_connections[part_first][part_second].emplace_back(neighbor_node, 0);
                            pair_2 = external_connections[part_first][part_second].back();
                        }

                        auto update_connection = [&](optional<pair<int, int>>& connection, int node_degree) {
                            if (degree_FC[node_degree] + connection.value().second < K) {
                                connection.value().second++;
                                cc_num[part_first][part_second] -= 0.5;
                            } else {
                                connection.value().second = max(0, K - degree_FC[node_degree]);
                            }
                        };

                        // 更新连接信息
                        update_connection(pair_1, node);
                        update_connection(pair_2, neighbor_node);

                        // 更新 external_connections 中的实际值
                        for (auto& p : external_connections[part_first][part_second]) {
                            if (p.first == pair_1.value().first) {
                                p.second = pair_1.value().second;
                            }
                            if (p.first == pair_2.value().first) {
                                p.second = pair_2.value().second;
                            }
                        }

                        pair<int, int> new_pair = {part_first, part_second};
                        if (find(affected_indices.begin(), affected_indices.end(), new_pair) == affected_indices.end()) {
                            affected_indices.emplace_back(new_pair);
                        }
                    }
                }
            }
            
            for (const auto& index : affected_indices) {
                int part_first = index.first;
                int part_second = index.second;
                if (!visited_eff[part_first][part_second]) {
                    if(cc_num[part_first][part_second] > 0){
                        cc_eff[part_first][part_second] = static_cast<double>(cc[part_first][part_second].size()) / cc_num[part_first][part_second];
                    }else{
                        cc_eff[part_first][part_second] = numeric_limits<double>::max();
                    }
                }
            }
            
            total_sum += cc_num[max_shell_index][max_component_index];
            selected_components.push_back({max_shell_index, max_component_index});
            final_solution.insert(final_solution.end(), cc[max_shell_index][max_component_index].begin(), cc[max_shell_index][max_component_index].end());
            added_nodes_count += cc[max_shell_index][max_component_index].size();
	    }
        ofstream output_file("sorted_" + graph_name);
        if (output_file.is_open()) {
            // 写入 b 的值
            output_file << b << endl;

            // 写入 selected_components 的数据
            for (const auto& component : selected_components) {
                output_file << component.first << " " << component.second << endl;
            }

            output_file.close();
        } else {
            cerr << "无法打开文件以写入数据。" << endl;
        }
    }

    // 最后处理final_solution中的节点连接，并输出添加的边
    vector<pair<int, int>> final_edges_to_add;
    sort(final_solution.begin(), final_solution.end(), [this](int a, int b) {
        return degree_FC[a] < degree_FC[b];
    });

    for (size_t i = 0; i < final_solution.size(); ++i) {
        for (size_t j = i + 1; j < final_solution.size(); ++j) {
            int u = final_solution[i];
            int v = final_solution[j];
            if (degree_FC[u] < K && degree_FC[v] < K) {
                if (eset.find(normalize_edge(u, v)) == eset.end() && find(final_edges_to_add.begin(), final_edges_to_add.end(), normalize_edge(u, v)) == final_edges_to_add.end()) {
                    final_edges_to_add.push_back(normalize_edge(u, v));
                    degree_FC[u]++;
                    degree_FC[v]++;
                }
            }
            if (degree_FC[u] >= K) {
				break;
			}
        }
    }
    
    cout<<"valid_edges:"<<final_edges_to_add.size()<<endl;
    
    // 随意连接到k_core_vertices中的节点
    for (int u : final_solution) {
        while (degree_FC[u] < K) {
            int v = k_core_vertices[rand() % k_core_vertices.size()];
            if (u != v && eset.find(normalize_edge(u, v)) == eset.end() && find(final_edges_to_add.begin(), final_edges_to_add.end(), normalize_edge(u, v)) == final_edges_to_add.end()) {
                final_edges_to_add.push_back(normalize_edge(u, v));
                degree_FC[u]++;
                degree_FC[v]++;
            }
        }
    }
    cout<<"edges:"<<final_edges_to_add.size()<<endl;
    cout << "added nodes:" << final_solution.size() << endl;
    return make_tuple(final_edges_to_add.size(), final_solution.size(), final_edges_to_add);
}

tuple<int, int, vector<pair<int, int>>> EfficientCMAlgorithm::greedy_solution_selection(int b) {
    // 创建一个参数，记录每个连通部分的效率（节点数量/边数）
    unordered_map<int, vector<double>> efficiencies;

    for (auto &layer : followers_number) {
        int flag = layer.first;
        for (int i = 0; i < layer.second.size(); ++i) {
            double edge_count = solution_condidates[flag][i].size() + static_cast<double>(biedge_number[flag][i]) / 2.0;
            if(edge_count == 0){
				double efficiency = -1;
				efficiencies[flag].push_back(efficiency);
			}else{
	            double efficiency = static_cast<double>(layer.second[i]) / edge_count;
	            efficiencies[flag].push_back(efficiency);				
			}
        }
    }

    unordered_set<int> used_nodes; // 记录已使用的节点
    vector<vector<pair<int, int>>> final_solution; // 最终解
    vector<int> final_edge_counts; // 记录最终解中每个连通部分的边数
    vector<int> final_bi_edge_counts; // 记录最终解中每个连通部分的半有效边数
    vector<vector<int>> final_bi_edge_items; // 记录最终解中每个连通部分的半有效边对应的节点
    vector<unordered_set<int>> final_nodes; // 记录最终解中每个连通部分的节点
    int total_edges = 0; // 添加的边的总数
    int count = 0;
    vector<pair<int, int>> all_edges; // 记录所有添加的边
    bool flag = false;
    int num = 0;

    while (true) {
        double max_efficiency = -1;
        int best_flag = -1, best_index = -1;

        // 找到效率最高的未访问的连通部分
        for (auto &layer : efficiencies) {
            int flag = layer.first;
            for (int i = 0; i < layer.second.size(); ++i) {
                if (layer.second[i] > max_efficiency) {
                    max_efficiency = layer.second[i];
                    best_flag = flag;
                    best_index = i;
                }
            }
        }

        // 如果没有未访问的部分，结束循环
        if (best_flag == -1) break;

        // 检查约束条件
        vector<int> &component = followers_items[best_flag][best_index];
        vector<int> bi_nodes = biedge_nodes[best_flag][best_index];
        vector<pair<int, int>> &edges = solution_condidates[best_flag][best_index];
        int bi_edge = biedge_number[best_flag][best_index];
        unordered_set<int> new_nodes(component.begin(), component.end());

        bool can_add = true;
        for (int node : new_nodes) {
            if (used_nodes.count(node)) {
                can_add = false;
                break;
            }
        }

        if (can_add && total_edges + edges.size() + customGreedyAlgorithm(final_bi_edge_counts, bi_edge) <= b) {
        	flag = true;
        	num = 0;
            // 直接加入解集
            total_edges += edges.size();
            final_solution.push_back(edges);
            final_edge_counts.push_back(edges.size());
            final_bi_edge_counts.push_back(bi_edge);
            final_bi_edge_items.push_back(bi_nodes);
            final_nodes.push_back(new_nodes);
            used_nodes.insert(new_nodes.begin(), new_nodes.end());
            all_edges.insert(all_edges.end(), edges.begin(), edges.end()); // 添加边到 all_edges
        } else {
            // 检查包含关系并处理
            vector<int> contained_indices;
            int total_existing_edges = 0;

            for (int i = 0; i < final_solution.size(); ++i) {
                unordered_set<int> &existing_nodes = final_nodes[i];
                bool is_contained = true;

                for (int node : existing_nodes) {
                    if (new_nodes.count(node) == 0) {
                        is_contained = false;
                        break;
                    }
                }

                if (is_contained) {
                    contained_indices.push_back(i);
                    total_existing_edges += final_edge_counts[i];
                }
            }

            // 对 contained_indices 进行降序排序
            sort(contained_indices.rbegin(), contained_indices.rend());

            if (!contained_indices.empty() && total_edges - total_existing_edges + edges.size() + customGreedyAlgorithm(final_bi_edge_counts, bi_edge) <= b) {
            	flag = true;
            	num = 0;
                for (int index : contained_indices) {
                    total_edges -= final_edge_counts[index];
                    for (int node : final_nodes[index]) {
                        if (used_nodes.count(node)) {
                            used_nodes.erase(node);
                        }
                    }

                    // 删除包含的边
                    for (const auto &edge : final_solution[index]) {
                        auto it = std::find(all_edges.begin(), all_edges.end(), edge);
                        if (it != all_edges.end()) {
                            all_edges.erase(it);
                        }
                    }

                    final_solution.erase(final_solution.begin() + index);
                    final_edge_counts.erase(final_edge_counts.begin() + index);
                    final_bi_edge_items.erase(final_bi_edge_items.begin()+index);
                    final_bi_edge_counts.erase(final_bi_edge_counts.begin() + index);
                    final_nodes.erase(final_nodes.begin() + index);
                }

                total_edges += edges.size();
                final_solution.push_back(edges);
                final_edge_counts.push_back(edges.size());
                final_bi_edge_counts.push_back(bi_edge);
                final_bi_edge_items.push_back(bi_nodes);
                final_nodes.push_back(new_nodes);
                used_nodes.insert(new_nodes.begin(), new_nodes.end());
                all_edges.insert(all_edges.end(), edges.begin(), edges.end()); // 添加边到 all_edges
            }else{
				flag = false;
				num++;
			}
        }

        // 标记该连通部分为已访问
        efficiencies[best_flag][best_index] = -1;
        if(num == 10){
			break;
		}
    }


        // 对final_bi_edge_counts和final_bi_edge_items进行降序排序
        vector<pair<int, vector<int>>> bi_edge_pairs;
        for (int i = 0; i < final_bi_edge_counts.size(); ++i) {
            bi_edge_pairs.push_back({final_bi_edge_counts[i], final_bi_edge_items[i]});
        }

        sort(bi_edge_pairs.begin(), bi_edge_pairs.end(), [](const pair<int, vector<int>>& a, const pair<int, vector<int>>& b) {
            return a.first > b.first;
        });

        // 尝试连接节点
        int num_ry = 0;
        unordered_set<int> processed_nodes;
        for (size_t i = 0; i < bi_edge_pairs.size(); ++i) {
            auto &current_component_nodes = bi_edge_pairs[i].second;
            for (int node : current_component_nodes) {
                if (processed_nodes.count(node)) continue;
                bool connected = false;

                // 尝试与不同component中的节点连接
                for (size_t j = i + 1; j < bi_edge_pairs.size(); ++j) {
                    auto &other_component_nodes = bi_edge_pairs[j].second;
                    for (int other_node : other_component_nodes) {
                        if (!processed_nodes.count(other_node)) {
                            all_edges.push_back({node, other_node});
                            processed_nodes.insert(other_node);
                            connected = true;
                            break;
                        }
                    }
                    if (connected) break;
                }

                // 如果无法连接，则与k_core_vertices中的任意节点连接
                if (!connected && !k_core_vertices.empty()) {
                	num_ry++;
                    all_edges.push_back({node, k_core_vertices.back()});
                }
                processed_nodes.insert(node);
            }
        }

    int total_nodes = used_nodes.size();
    return make_tuple(all_edges.size(), total_nodes, all_edges); // 返回最终解添加的边的数量和节点的数量和所有的边
}

void EfficientCMAlgorithm::partition_shell() {
    cc.clear(); // 清空之前的结果
    k_degree.clear(); // 清空k_degree

    // 对于每个核数大于等于K-i且小于K的节点集合，从1开始
    for (int i = 1; i < depth + 1; ++i) {
        unordered_map<int, bool> visited; // 记录节点是否被访问
        queue<int> q; // 用于广度优先搜索

        // 遍历每个节点，根据其核数进行处理
        for (const auto &node : coreness) {
            int u = node.first;
            int core_value = node.second;

            // 检查节点是否在当前层的范围内
            if (core_value >= K - i && core_value < K && !visited[u]) {
                // 初始化一个新的连通部分
                vector<int> component;
                q.push(u);
                visited[u] = true;

                // 广度优先搜索，找到所有连通的节点
                while (!q.empty()) {
                    int v = q.front();
                    q.pop();
                    component.push_back(v);

                    // 遍历v的所有邻居
                    for (const auto &neighbor : g[v]) {
                        int w = neighbor.u;
                        if (!visited[w] && coreness[w] >= K - i && coreness[w] < K) {
                            q.push(w);
                            visited[w] = true;
                        }
                    }
                }

                // 计算该连通部分中每个节点在由component和所有coreness大于等于K的节点构成的诱导子图中的度数
                unordered_set<int> relevant_nodes(component.begin(), component.end());
                for (const auto &node : coreness) {
                    if (node.second >= K) {
                        relevant_nodes.insert(node.first);
                    }
                }

                for (int v : component) {
                    int degree = 0;
                    for (const auto &neighbor : g[v]) {
                        int w = neighbor.u;
                        if (relevant_nodes.count(w)) {
                            degree++;
                        }
                    }
                    k_degree[i][v] = degree;
                }

                // 对找到的连通部分根据 k_degree 进行排序
                sort(component.begin(), component.end(), [i, this](int a, int b) {
                    return k_degree[i][a] < k_degree[i][b];
                });

                // 将找到的连通部分添加到结果中
                cc[i].push_back(component);
            }
        }
    }
}

void EfficientCMAlgorithm::partition_shell_FC(const string &graph_name) {
    cc.clear();
    cc_num.clear();
    degree_FC.clear();
    
    // 初始化所有节点的 node_partition 为 {-1, -1}
    for (const auto& node : g) {
        node_partition[node.first] = {-1, -1};
    }

    // 始终按照原逻辑进行初始化，不使用文件数据
    // 初始化 degree_FC
    for (const auto& node : g) {
        int u = node.first;
        int coreness_u = coreness[u];
        int degree_u = 0;
        for (const auto& neighbor : node.second) {
            int v = neighbor.u;
            if (coreness[v] == coreness_u || coreness[v] >= K) {
                degree_u++;
            }
        }
        degree_FC[u] = degree_u;
    }

    // 计算 cc 和 cc_num
    unordered_map<int, bool> visited;
    for (const auto& node : g) {
        int u = node.first;
        int core = coreness[u];

        if (core >= K || core < K - depth || visited[u]) {
            continue; // 节点的 coreness 必须小于 K 并且大于等于 K - depth
        }

        int shell_index = K - core - 1; // 节点所在的层

        if (cc.find(shell_index) == cc.end()) {
            cc[shell_index] = vector<vector<int>>();
            cc_num[shell_index] = vector<double>();
        }

        vector<int> new_component;
        dfs(u, core, new_component, visited);
        cc[shell_index].push_back(new_component);

        int degree_sum = 0;
        for (int v : new_component) {
            int degree = degree_FC[v] > K ? 0 : K - degree_FC[v];
            degree_sum += degree;
        }
        cc_num[shell_index].push_back(degree_sum / 2.0);
    }

    // 初始化 node_partition
    for (const auto &cc_entry : cc) {
        int cc_key = cc_entry.first;
        const vector<vector<int>> &cc_value = cc_entry.second;

        for (int i = 0; i < cc_value.size(); ++i) {
            const vector<int> &node_list = cc_value[i];
            for (const auto &node : node_list) {
                node_partition[node] = make_pair(cc_key, i);
            }
        }
    }
}

tuple<int, int, vector<pair<int, int>>> EfficientCMAlgorithm::EfficientCM(int b) {
	partition_shell();
    auto start_FG = std::chrono::high_resolution_clock::now();
	complete_conversion();
    auto end_FG = std::chrono::high_resolution_clock::now();
    double total_t_FG = double(std::chrono::duration_cast<std::chrono::milliseconds>(end_FG - start_FG).count()) / 1000;
    cout<<"complete growth time:"<<total_t_FG<<endl;
	return greedy_solution_selection(b);
}


tuple<int, int, vector<pair<int, int>>> EfficientCMAlgorithm::EfficientCM_FG(int b, const string &graph_name) {
	// cout<<"b num:"<< b <<endl;
    partition_shell_FC(graph_name);
    // cout<<"done"<<endl;
    auto start_FG = std::chrono::high_resolution_clock::now();
    tuple<int, int, vector<pair<int, int>>> res_FG = complete_conversion_FC(b, graph_name);
    auto end_FG = std::chrono::high_resolution_clock::now();
    double total_t_FG = double(std::chrono::duration_cast<std::chrono::milliseconds>(end_FG - start_FG).count()) / 1000;
    cout<<"partial growth time:"<<total_t_FG<<endl;
    // cout<<"ok"<<endl;
	return res_FG;
}