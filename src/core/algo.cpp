#include "algo.h"

algo::algo() {}
algo::~algo() {}

/////////////////////////////////////////////////////////////////////
//////////////////////// UTILITIES //////////////////////////////////
/////////////////////////////////////////////////////////////////////

bool compare(const std::pair<int, int> &a, const std::pair<int, int> &b)
{
    return a.second > b.second;
}

int algo::minKey_index(std::map<int, double> key, std::map<int, bool> mstSet, int INT_MAX, int number_of_nodes)
{
    double min = INT_MAX;
    int min_index = -1;
    for (auto const &pair : key)
    {
        if (mstSet[pair.first] == false and key[pair.first] <= min)
        {
            min = key[pair.first];
            min_index = pair.first;
        }
    }
    return min_index;
}

int algo::find_in_subtree(std::map<int, std::pair<double, std::unordered_set<int>>> sub_trees, int who)
{
    for (auto const &pair : sub_trees)
    {
        if (sub_trees[pair.first].second.find(who) != sub_trees[pair.first].second.end())
            return pair.first;
    }
    return -1;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////// MAIN ////////////////////////////////////
/////////////////////////////////////////////////////////////////////

std::map<int, int> algo::primMST(graph G, std::vector<int> forced_nodes, double budget, bool weighted)
{
    std::map<int, int> parent;
    std::map<int, double> key;
    std::map<int, bool> mstSet;
    std::map<int, std::pair<double, std::unordered_set<int>>> sub_trees;
    std::vector<int> id_vertices = G.get_vertices();

    double _max = budget * 2;

    for (auto const &i : id_vertices)
    {
        key[i] = _max;
        mstSet[i] = false;
        parent[i] = -1;
    }

    for (auto const &i : forced_nodes)
    {
        parent[i] = 0;
        key[i] = G.dist(0, i) + (weighted ? G.get_weight_node(i) : 0);
        std::unordered_set<int> myset = {i};
        sub_trees[i] = std::make_pair(budget - key[i], myset);
    }
    parent[0] = -1;
    mstSet[0] = forced_nodes.size() > 1 ? true : false;

    ////////////////////////////////////////////

    for (int j = 0; j < G.get_n_nodes(); j++)
    {
        int u = minKey_index(key, mstSet, _max, G.get_n_nodes());
        if (u == -1)
        {
            break;
        }

        mstSet[u] = true;

        for (auto const &pair : key)
        {
            if (pair.first == u or find(forced_nodes.begin(), forced_nodes.end(), pair.first) != forced_nodes.end())
            {
                continue;
            }
            int index = find_in_subtree(sub_trees, u);

            if (mstSet[pair.first] == false and G.dist(0, pair.first) + (weighted ? G.get_weight_node(pair.first) : 0) < key[pair.first])
            {
                double dist_2 = G.dist(u, pair.first) + (weighted ? G.get_weight_node(pair.first) : 0);

                if (parent[pair.first] != -1)
                {
                    double dist_1 = G.dist(parent[pair.first], pair.first) + (weighted ? G.get_weight_node(pair.first) : 0);
                    double dist_3 = G.dist(pair.first, parent[pair.first]) + (weighted ? G.get_weight_node(pair.first) : 0);

                    int index_p = find_in_subtree(sub_trees, pair.first);
                    if (index == index_p)
                    {

                        if (sub_trees[index].first + dist_1 - dist_2 >= 0)
                        {
                            sub_trees[index_p].first += dist_3 - dist_2;
                            parent[pair.first] = u, key[pair.first] = dist_2;
                        }
                    }
                    else
                    {
                        if (sub_trees[index].first - dist_2 >= 0)
                        {
                            sub_trees[index_p].first += dist_3;
                            sub_trees[index_p].second.erase(pair.first);
                            sub_trees[index].first -= dist_2;
                            sub_trees[index].second.insert(pair.first);
                            parent[pair.first] = u, key[pair.first] = dist_2;
                        }
                    }
                }
                else
                {
                    if (sub_trees[index].first - dist_2 >= 0)
                    {
                        sub_trees[index].first -= dist_2;
                        sub_trees[index].second.insert(pair.first);
                        parent[pair.first] = u, key[pair.first] = dist_2;
                    }
                }
            }
        }
    }

    return parent;
}

void algo::DFSUtil(graph G, int v, std::map<int, int> tree, std::map<int, bool> visited, std::vector<int> &sol, double &cost_cycle, double budget, bool weighted)
{
    if (sol.size() > 0)
    {
        double dist_1 = G.dist(sol[sol.size() - 1], v) + (weighted ? G.get_weight_node(v) : 0);

        if (cost_cycle + dist_1 + G.dist(v, 0) > budget)
            return;
        cost_cycle += dist_1;
    }
    sol.push_back(v);
    visited[v] = true;
    for (auto const &pair : tree)
    {
        if (tree[pair.first] == v and not visited[pair.first])
        {
            DFSUtil(G, pair.first, tree, visited, sol, cost_cycle, budget, weighted);
        }
    }
}

std::vector<int> algo::find_TSP(graph G, double budget, int start, std::map<int, int> tree, bool weighted)
{
    std::vector<int> sol;
    std::map<int, bool> visited;
    for (auto const &pair : tree)
    {
        visited[pair.first] = false;
    }

    double cost_cycle = 0;
    if (start != 0)
    {
        cost_cycle += G.dist(0, start) + G.get_weight_node(start);
        visited[0] = true;
    }
    DFSUtil(G, start, tree, visited, sol, cost_cycle, budget, weighted);

    return sol;
}