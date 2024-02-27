#include <assert.h>
#include <iostream>

#include "core/user_input.h"
#include "core/graph.h"
#include "simulator/simulator_cc.h"

void sortrows(std::vector<std::vector<bool>> &matrix, int col)
{
    std::sort(matrix.begin(),
              matrix.end(),
              [col](const std::vector<bool> &lhs, const std::vector<bool> &rhs)
              {
                  return lhs[col] < rhs[col];
              });
}

void remove_duplicates_from_matrix(std::vector<std::string> &temp)
{
    std::sort(temp.begin(), temp.end());
    temp.erase(std::unique(temp.begin(), temp.end()), temp.end());
}

inline const char *const BoolToString(bool b)
{
    return b ? "1" : "0";
}

void cycle_generation_starter(graph G, int n_drones, double budget, long seed, int MAX_COUNTER, double radius_distance, std::string output_data_file, std::string path_cycles)
{
    std::vector<std::vector<std::string>> solutions;
    std::vector<std::string> algs = {
        "Greedy-TOP",      // A
        "psuedoKim",       // B
        "Greedy-ab",       // C
        "TOP-h",           // D
        "Greedy-CE",       // E
        "Greedy-ab-0",     // F
        "Greedy-ab-0.25",  // G
        "Greedy-ab-0.5",   // H
        "Greedy-ab-0.75",  // I
        "Greedy-ab-1",     // L
        "Greedy-TOP-B+10", // M
        "Greedy-TOP-B+20", // N
        "Greedy-TOP-B+30", // O
        "Greedy-TOP-B+40", // P
        "Greedy-TOP-B+50", // Q
        "psuedoKim-B-10",  // R
        "psuedoKim-B-20",  // S
        "psuedoKim-B-30"   // T
        "psuedoKim-B-40"   // U
        "psuedoKim-B-50"   // V
    };

    std::vector<std::string> which_heur_label = {
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "L",
        "M",
        "N",
        "O",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "U",
        "V"};

    std::cout << std::endl;

    // 
    simulator_cc sim = simulator_cc(G, n_drones, budget, seed);
    assert(sim.check_feasibility());

    MAX_COUNTER = 100;
    for (int i = 0; i < 5; i++)
    {
        std::vector<std::string> output_cycles = sim.cycle_generation(MAX_COUNTER, i, budget, radius_distance, 0, which_heur_label[i], output_data_file);

        remove_duplicates_from_matrix(output_cycles);
        solutions.push_back(output_cycles);

        std::ofstream outfile;
        outfile.open(output_data_file, std::ios_base::app);

        outfile << " | unici: " << output_cycles.size() << " | s: " << std::to_string(seed) << std::endl;

        outfile.close();
        std::cout << " | unici: " << output_cycles.size() << std::endl;
    }

    double alpha = 0;
    for (int i = 5; i < 10; i++)
    {

        std::vector<std::string> output_cycles = sim.cycle_generation(MAX_COUNTER, 5, budget, radius_distance, alpha, which_heur_label[i], output_data_file);

        remove_duplicates_from_matrix(output_cycles);
        solutions.push_back(output_cycles);

        std::ofstream outfile;
        outfile.open(output_data_file, std::ios_base::app);

        outfile << " | unici: " << output_cycles.size() << " | s: " << std::to_string(seed) << std::endl;

        outfile.close();
        alpha += 0.25;
        std::cout << " | unici: " << output_cycles.size() << std::endl;
    }

    int increment = 0;
    for (int i = 10; i < 15; i++)
    {
        double new_budget = budget + budget * (0.1 * (++increment));

        std::vector<std::string> output_cycles = sim.cycle_generation(MAX_COUNTER, 0, new_budget, radius_distance, 0, which_heur_label[i], output_data_file);

        remove_duplicates_from_matrix(output_cycles);
        solutions.push_back(output_cycles);

        std::ofstream outfile;
        outfile.open(output_data_file, std::ios_base::app);

        outfile << " | unici: " << output_cycles.size() << " | " << std::to_string(seed) << std::endl;

        outfile.close();
        std::cout << " | unici: " << output_cycles.size() << std::endl;
    }

    increment = 0;
    for (int i = 15; i < 20; i++)
    {

        double new_budget = budget - budget * (0.1 * (++increment));

        std::vector<std::string> output_cycles = sim.cycle_generation(MAX_COUNTER, 1, new_budget, radius_distance, 0, which_heur_label[i], output_data_file);

        remove_duplicates_from_matrix(output_cycles);
        solutions.push_back(output_cycles);
        std::ofstream outfile;
        outfile.open(output_data_file, std::ios_base::app);

        outfile << " | unici: " << output_cycles.size() << " | s: " << std::to_string(seed) << std::endl;

        outfile.close();
        std::cout << " | unici: " << output_cycles.size() << std::endl;
    }

    std::vector<std::string> all_together;
    std::cout << std::endl
              << std::endl;
    for (int i = 0; i < (int)solutions.size(); i++)
    {
        for (int j = 0; j < (int)solutions.size(); j++)
        {
            std::vector<std::string> temp;
            temp.reserve(temp.size() + solutions[i].size() + solutions[j].size());
            temp.insert(temp.end(), solutions[i].begin(), solutions[i].end());
            temp.insert(temp.end(), solutions[j].begin(), solutions[j].end());

            remove_duplicates_from_matrix(temp);

            double i_plus_j = solutions[i].size() + solutions[j].size();
            double i_U_j = temp.size();
            double i_int_j = (i_plus_j - i_U_j);

            std::cout.precision(2);
            std::cout << (i_int_j / i_U_j) * 100 << "\t";
        }
        std::cout << std::endl;

        all_together.reserve(all_together.size() + solutions[i].size());
        all_together.insert(all_together.end(), solutions[i].begin(), solutions[i].end());
    }

    std::cout << std::endl
              << all_together.size() << " ";
    remove_duplicates_from_matrix(all_together);
    std::cout << all_together.size() << std::endl;

    std::cout << sim.get_output_hashmap_size() << std::endl;

    assert((int)sim.get_output_hashmap_size() == (int)all_together.size());

    sim.print_output_hashmap_to_file(path_cycles);
}

void simulation_compute_cycles(input n_input)
{
    int radius_int = 150;
    double radius_distance = radius_int / 1000.0;
    std::cout << "radius distance: " << radius_distance << std::endl;
    int MAX_COUNTER = 100;
    std::string path_graph = "data/graph/graph_n" + std::to_string(n_input.n_nodes) + "_s" + std::to_string(n_input.seed) + ".csv";
    std::string path_cycles = "data/output/cycles_n" + std::to_string(n_input.n_nodes) + "_b" + std::to_string((int)n_input.budget) + "_q" + std::to_string(n_input.n_drones) + "_mc" + std::to_string(MAX_COUNTER) + "_r" + std::to_string(radius_int) + "_s" + std::to_string(n_input.seed) + ".csv";
    std::string output_data_file = "data/output/output_data" + std::to_string(n_input.n_nodes) + "_b" + std::to_string((int)n_input.budget) + "_q" + std::to_string(n_input.n_drones) + "_mc" + std::to_string(MAX_COUNTER) + "_r" + std::to_string(radius_int) + ".csv";

    graph G = graph(1, 1);
    G.erase_graph();
    G.create_random_graph(n_input.n_nodes, 0, n_input.seed);
    G.print_graph_to_file(path_graph);

    cycle_generation_starter(G, n_input.n_drones, n_input.budget, n_input.seed, MAX_COUNTER, radius_distance, output_data_file, path_cycles);
}

void simulation_compute_cycles_with_grids(input n_input)
{
    int radius_int = 150;
    double radius_distance = radius_int / 1000.0;
    std::cout << "radius distance: " << radius_distance << std::endl;
    int MAX_COUNTER = 100;

    double length_grid = n_input.grid_length / 1000.0;

    graph G = graph(1, 1);
    G.erase_graph();
    G.read_graph_from_file(n_input.graph_file);
    G.set_radius_distance(radius_distance);
    G.add_grid(length_grid);
    std::cout << G << std::endl;

    std::string path_graph = "data/graph/grid_graph_n" + std::to_string(G.get_n_nodes() - G.get_n_nodes_grid() - 1) + "_gl" + std::to_string(n_input.grid_length) + "_s" + std::to_string(n_input.seed) + ".csv";
    std::string path_cycles = "data/output/grid_cycles_n" + std::to_string(G.get_n_nodes() - G.get_n_nodes_grid() - 1) + "_b" + std::to_string((int)n_input.budget) + "_q" + std::to_string(n_input.n_drones) + "_mc" + std::to_string(MAX_COUNTER) + "_r" + std::to_string(radius_int) + "_gl" + std::to_string(n_input.grid_length) + "_s" + std::to_string(n_input.seed) + ".csv";
    std::string output_data_file = "data/output/grid_output_data" + std::to_string(G.get_n_nodes() - G.get_n_nodes_grid() - 1) + "_b" + std::to_string((int)n_input.budget) + "_q" + std::to_string(n_input.n_drones) + "_mc" + std::to_string(MAX_COUNTER) + "_r" + std::to_string(radius_int) + "_gl" + std::to_string(n_input.grid_length) + ".csv";

    G.print_graph_to_file(path_graph);

    cycle_generation_starter(G, n_input.n_drones, n_input.budget, n_input.seed, MAX_COUNTER, radius_distance, output_data_file, path_cycles);
}

int main(int argc, char **argv)
{
    input n_input = userinput::read_user_input(argc, argv);

    switch (n_input.experiment)
    {

    case 1: // --generate-trips
        simulation_compute_cycles(n_input);
        break;
    case 2: // --generate-trips grid
        simulation_compute_cycles_with_grids(n_input);
        break;
    default:
        std::cerr << "\n[ERROR main] experiment to run not present\n";
        assert(0);
        break;
    }
    return 0;
}