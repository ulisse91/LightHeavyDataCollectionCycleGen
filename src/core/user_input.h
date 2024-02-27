#ifndef USERINPUT_H
#define USERINPUT_H

#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <map>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

struct input
{
    /// Default values
    int n_nodes = 10;
    int n_drones = 1;
    int grid_length = 150;
    double budget = 1;
    std::string graph_file = "";
    std::string drones_file = "";

    int seed = 0;
    int experiment = 0;

    input(int _n_nodes, int _n_drones, int _grid_length, int _n_batteries, double _budget, double _prob_sigma_prime, std::string _graph_file, std::string _distrib) : n_nodes(_n_nodes), n_drones(_n_drones), grid_length(_grid_length), budget(_budget), graph_file(_graph_file) {}
    input() = default;
};

class userinput
{

private:
public:
    userinput();
    ~userinput();

    static input read_user_input(int argc, char **argv);
    static void print_input(input n_input);
};

#endif // USERINPUT_H