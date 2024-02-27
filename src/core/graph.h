#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <random>
#include <assert.h>
#include <fstream>
#include <unordered_set>
#include <map>
#include <stack>
#include <cmath>
#include <algorithm>

struct node
{
    int id;
    double x, y;
    double weak_coverage_cost = 0;
    double strong_coverage_cost = 0;
    int profit;

    node(int _id, double _x, double _y, double _weak_coverage_cost, double _strong_coverage_cost, int _profit) : id(_id), x(_x), y(_y), weak_coverage_cost(_weak_coverage_cost), strong_coverage_cost(_strong_coverage_cost), profit(_profit) {}
    node() = default;
};

class graph
{
    friend std::ostream &operator<<(std::ostream &os, const graph &G);

private:
    int area_x;
    int area_y;
    int n_nodes;
    int n_nodes_grid;
    double radius_distance;

    std::map<int, node> vertices;

    bool check_double_node(double _x, double _y);
    bool isInsidePolygon(const std::vector<node> &convexHull, double x, double y);
    std::vector<node> findConvexHull(std::map<int, node> &vertices);
    bool touch_some_node(double x, double y);

public:
    graph();
    graph(int _area_x, int _area_y);
    graph(int _area_x, int _area_y, int number_of_depot);
    ~graph();

    void set_radius_distance(double _radius_distance);

    int get_n_nodes();
    int get_n_nodes_grid();
    double get_strong_coverage_node(int id);
    double get_weight_node(int id);
    double get_weak_coverage_node(int id);
    double get_profit_node(int id);

    double get_coord_x(int id);
    double get_coord_y(int id);
    std::vector<int> get_vertices();
    std::unordered_set<int> get_real_vertices_set();
    std::unordered_set<int> get_vertices_set();

    double dist(int u, int v);
    double distw(int u, int v);
    double dist_eucl(int u, int v);
    double dist_eucl(int id_node, double node_x, double node_y);

    int read_graph_from_file(std::string file);

    void create_random_graph(int number_of_nodes, double max_weight);
    void create_random_graph(int number_of_nodes, double max_weight, long seed);
    void add_grid(double lenght_grid);

    int add_node(int id, double _x, double _y, double weak_coverage_cost, double strong_coverage_cost, int profit);
    int add_node(double _x, double _y, double weak_coverage_cost, double strong_coverage_cost, int profit);

    void erase_graph();

    void print_graph_to_file(std::string path);
};

#endif // GRAPH_H