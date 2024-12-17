#ifndef SIMULATOR_CC_H
#define SIMULATOR_CC_H

#include <iostream>
#include <vector>
#include <assert.h>
#include <utility>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <chrono>
#include <iterator>
#include <tuple>
#include <algorithm>

#include "../core/graph.h"
#include "../core/algo.h"

class simulator_cc
{

private:
    struct Cycle
    {
        std::string cycle = "";
        double budget = 0.;
        std::unordered_set<std::string> label = {};
    };

    graph G;
    int n_drones;
    double budget;
    long seed;

    double cost_budget_sequence(graph G, std::vector<int> sequence);
    double profit_path_sequence(graph G, std::vector<int> sequence, std::unordered_set<int> &light_covered_vertices, std::vector<std::unordered_set<int>> neigh, double budget);
    std::unordered_set<int> set_difference(std::unordered_set<int> main, std::vector<int> minus);
    double cost_cycle_OP(std::vector<int> _temp, double given_budget);

    std::vector<std::vector<int>> greedy_alphabeta_innerloop(double given_budget, std::vector<std::unordered_set<int>> neigh, long given_seed, double alpha);
    std::map<std::string, Cycle> output_Cycle_temp;

public:
    simulator_cc(graph _G, int drones, double budget, long seed);
    ~simulator_cc();

    bool check_feasibility();

    void increase_seed();

    std::vector<std::vector<int>> greedy_top(double budget, long seed);
    std::vector<std::vector<int>> prim_based(double given_budget, long seed);
    std::vector<std::vector<int>> top_heur(double given_budget, long given_seed);
    std::vector<std::vector<int>> greedy_alphabeta(double given_budget, std::vector<std::unordered_set<int>> neigh, long given_seed);
    double cost_path_sequence(graph G, std::vector<int> sequence);
    std::vector<std::vector<int>> greedy_ce(double given_budget, std::vector<std::unordered_set<int>> neigh, long given_seed);

    int get_output_hashmap_size();
    void print_output_hashmap_to_file(std::string path_to_file);

    std::vector<std::string> cycle_generation(int MAX_COUNTER, int which_alg, double budget_to_test, double radius_distance, double alpha, std::string which_heur, std::string output_data_file);
};

#endif // SIMULATOR_CC_H