#include "simulator_cc.h"

simulator_cc::simulator_cc(graph _G, int _drones, double _budget, long _seed)
{
    this->G = _G;
    this->n_drones = _drones;
    this->budget = _budget;
    this->seed = _seed;
}

simulator_cc::~simulator_cc() {}

void simulator_cc::increase_seed() {
    this->seed++;
}

bool simulator_cc::check_feasibility()
{
    for (auto &v : this->G.get_vertices())
        if (2 * this->G.distw(0, v) > this->budget)
            return false;
    return true;
}

double simulator_cc::cost_path_sequence(graph G, std::vector<int> sequence)
{
    double sum_of_elems = 0;
    if (sequence.size() > 1)
    {
        for (size_t i = 0; i < sequence.size() - 1; i++)
        {
            sum_of_elems += G.dist(sequence[i], sequence[i + 1]);
        }
        sum_of_elems += G.dist(sequence[sequence.size() - 1], 0);
    }
    return sum_of_elems;
}

double simulator_cc::profit_path_sequence(graph G, std::vector<int> sequence, std::unordered_set<int> &light_covered_vertices, std::vector<std::unordered_set<int>> neigh, double budget)
{
    double profit_gained = 0;
    double remaining_budget = budget;
    sequence.push_back(0);
    double current_path_cost = 0;
    int temp_size = light_covered_vertices.size();

    if (sequence.size() > 1)
    {
        for (size_t i = 0; i < sequence.size() - 1; i++)
        {
            std::vector<int> how_many_light_covered = {};
            for (int element : neigh[sequence[i + 1]])
            {
                if (light_covered_vertices.count(element) > 0)
                {
                    how_many_light_covered.push_back(element);
                }
            }
            if (light_covered_vertices.count(sequence[i + 1]) > 0)
            {
                how_many_light_covered.push_back(sequence[i + 1]);
            }

            remaining_budget -= (G.dist(sequence[i], sequence[i + 1]) + G.get_weak_coverage_node(sequence[i + 1]) * how_many_light_covered.size());
            current_path_cost += G.dist(sequence[i], sequence[i + 1]);
            light_covered_vertices = set_difference(light_covered_vertices, how_many_light_covered);
        }

        std::unordered_set<int> nodes_set(sequence.begin(), sequence.end());
        while (remaining_budget > 0 and nodes_set.size() > 0)
        {
            double max_ratio = 0;
            int which_node = -1;
            for (const int &i : nodes_set)
            {
                if (max_ratio < G.get_profit_node(i) / G.get_strong_coverage_node(i) and remaining_budget - G.get_strong_coverage_node(i) >= 0)
                {
                    max_ratio = G.get_profit_node(i) / G.get_strong_coverage_node(i);
                    which_node = i;
                }
            }

            if (which_node != -1)
            {
                remaining_budget -= G.get_strong_coverage_node(which_node);
                profit_gained += G.get_profit_node(which_node);
                nodes_set.erase(which_node);
            }
            else
            {
                break;
            }
        }
    }
    return profit_gained;
}

////////////////////////////////////
////////////////////////////////////
////////////////////////////////////

std::vector<std::vector<int>> simulator_cc::greedy_top(double given_budget, long given_seed)
{
    std::unordered_set<int> graph_vertices = this->G.get_vertices_set();
    std::vector<std::vector<int>> curr_sol(this->n_drones, {0});
    std::vector<double> residual_budget(this->n_drones, given_budget);

    bool round_finished = false;

    std::uniform_int_distribution<int> unif_4(0, 2);
    std::mt19937 re(given_seed);

    while (!round_finished)
    {
        round_finished = true;
        for (int drone = 0; drone < this->n_drones; drone++)
        {
            int next_step = -1;
            int last_step = curr_sol[drone][curr_sol[drone].size() - 1];
            double min_dist = 0;

            std::vector<std::pair<int, double>> temp_neigh = {{-1, 0}, {-1, 0}, {-1, 0}};

            for (auto const &node : graph_vertices)
            {
                double frac = 1;
                if ((this->G.get_profit_node(node) / (frac * this->G.dist(last_step, node) + this->G.get_weight_node(node)) > min_dist) and this->G.dist(last_step, node) + this->G.get_weight_node(node) + this->G.dist(0, node) < residual_budget[drone])
                {
                    temp_neigh[2].first = node;
                    temp_neigh[2].second = this->G.get_profit_node(node) / (frac * this->G.dist(last_step, node) + this->G.get_weight_node(node));
                    min_dist = temp_neigh[2].second;

                    std::sort(temp_neigh.begin(), temp_neigh.end(), [](auto &left, auto &right)
                              { return left.second > right.second; });
                }
            }

            int coin = unif_4(re);
            assert(0 <= coin and coin < (int)temp_neigh.size());
            next_step = temp_neigh[coin].first;

            if (next_step != -1)
            {
                round_finished = false;
                curr_sol[drone].push_back(next_step);
                if (graph_vertices.find(next_step) != graph_vertices.end())
                    graph_vertices.erase(graph_vertices.find(next_step));
                residual_budget[drone] -= (this->G.dist(last_step, next_step) + this->G.get_weight_node(next_step));
            }
        }
    }
    return curr_sol;
}

/////////////////////////////////////////////
/////////////////////////////////////////////
///////////// KIM

double simulator_cc::cost_budget_sequence(graph G, std::vector<int> sequence)
{
    double sum_of_elems = 0;
    if (sequence.size() > 1)
    {
        for (size_t i = 0; i < sequence.size() - 1; i++)
        {
            sum_of_elems += G.dist(sequence[i], sequence[i + 1]) + G.get_weight_node(i + 1);
        }
    }
    return sum_of_elems;
}

std::vector<std::vector<int>> simulator_cc::prim_based(double given_budget, long given_seed)
{
    std::unordered_set<int> graph_vertices = this->G.get_vertices_set();
    std::vector<std::vector<int>> curr_sol(this->n_drones, {0});

    std::map<int, int> tree;
    std::vector<int> centers;

    centers = G.get_vertices();
    std::shuffle(centers.begin(), centers.end(), std::default_random_engine(given_seed));

    centers.erase(centers.begin() + this->n_drones, centers.end());

    if (this->n_drones > 1)
        centers.push_back(0);

    tree = algo::primMST(G, centers, given_budget, false);

    if (this->n_drones > 1)
        centers.erase(std::remove(centers.begin(), centers.end(), 0), centers.end());

    for (int i = 0; i < (int)centers.size(); i++)
    {
        std::vector<int> tsp_i = algo::find_TSP(G, given_budget, centers[i], tree, false);
        curr_sol[i].insert(curr_sol[i].end(), tsp_i.begin(), tsp_i.end());
    }

    return curr_sol;
}

/////////////////////////////////////////////
/////////////////////////////////////////////
/////////////////////////////////////////////

std::vector<std::vector<int>> simulator_cc::greedy_alphabeta(double given_budget, std::vector<std::unordered_set<int>> neigh, long given_seed)
{
    std::uniform_real_distribution<double> unif_3(0, 1);
    // Mersenne Twister: Good quality random number generator
    std::mt19937 re(given_seed);

    double alpha = unif_3(re);

    return greedy_alphabeta_innerloop(given_budget, neigh, given_seed, alpha);
}

std::vector<std::vector<int>> simulator_cc::greedy_alphabeta_innerloop(double given_budget, std::vector<std::unordered_set<int>> neigh, long given_seed, double alpha)
{
    std::unordered_set<int> graph_vertices = this->G.get_vertices_set();
    std::vector<std::vector<int>> curr_sol(this->n_drones, {0});
    std::vector<double> residual_budget(this->n_drones, given_budget);
    std::unordered_set<int> light_covered_vertices = this->G.get_real_vertices_set();

    std::uniform_int_distribution<int> unif_4(0, 2);
    // Mersenne Twister: Good quality random number generator
    std::mt19937 re(given_seed);

    double beta = 1 - alpha;
    bool round_finished = false;

    while (!round_finished)
    {
        round_finished = true;
        for (int drone = 0; drone < this->n_drones; drone++)
        {
            int next_step = -1;
            int last_step = curr_sol[drone][curr_sol[drone].size() - 1];
            double min_dist = 0;

            std::vector<std::pair<int, double>> temp_neigh = {{-1, 0}, {-1, 0}, {-1, 0}};

            for (auto const &node : graph_vertices)
            {
                std::unordered_set<int> how_many_light_covered = {};
                for (int element : neigh[node])
                {
                    if (light_covered_vertices.count(element) > 0)
                    {
                        how_many_light_covered.insert(element);
                    }
                }
                if (light_covered_vertices.count(node) > 0)
                    how_many_light_covered.insert(node);

                double ratio_num = alpha * this->G.get_profit_node(node) + beta * how_many_light_covered.size();

                double ratio_den = this->G.dist(last_step, node) + this->G.get_weight_node(node) + this->G.get_weak_coverage_node(node) * how_many_light_covered.size();

                double budget_to_use = this->G.dist(last_step, node) + this->G.get_weight_node(node) + this->G.dist(0, node) + this->G.get_weak_coverage_node(node) * how_many_light_covered.size();

                if ((ratio_num / ratio_den > min_dist) and budget_to_use < residual_budget[drone])
                {

                    temp_neigh[2].first = node;
                    temp_neigh[2].second = ratio_num / ratio_den;
                    min_dist = temp_neigh[2].second;

                    std::sort(temp_neigh.begin(), temp_neigh.end(), [](auto &left, auto &right)
                              { return left.second > right.second; });
                }
            }

            int coin = unif_4(re);
            assert(0 <= coin and coin < (int)temp_neigh.size());

            next_step = temp_neigh[coin].first;

            if (next_step != -1)
            {
                std::vector<int> how_many_light_covered = {};
                for (int element : neigh[next_step])
                {
                    if (light_covered_vertices.count(element) > 0)
                    {
                        how_many_light_covered.push_back(element);
                    }
                }
                if (light_covered_vertices.count(next_step) > 0)
                {
                    how_many_light_covered.push_back(next_step);
                }

                light_covered_vertices = set_difference(light_covered_vertices, how_many_light_covered);

                round_finished = false;
                curr_sol[drone].push_back(next_step);
                if (graph_vertices.find(next_step) != graph_vertices.end())
                    graph_vertices.erase(graph_vertices.find(next_step));
                residual_budget[drone] -= (this->G.dist(last_step, next_step) + this->G.get_weight_node(next_step) + this->G.get_weak_coverage_node(next_step) * how_many_light_covered.size());
            }
        }
    }
    return curr_sol;
}

/////////////////////////////////////////////
/////////////////////////////////////////////
/////////////////////////////////////////////

std::vector<std::vector<int>> simulator_cc::greedy_ce(double given_budget, std::vector<std::unordered_set<int>> neigh, long given_seed)
{
    std::unordered_set<int> graph_vertices = this->G.get_vertices_set();
    std::vector<std::vector<int>> curr_sol(this->n_drones, {0});
    std::vector<double> residual_budget(this->n_drones, given_budget);
    std::unordered_set<int> light_covered_vertices = this->G.get_real_vertices_set();
    std::unordered_set<int> real_vertices_set = this->G.get_real_vertices_set();

    bool round_finished = false;
    std::uniform_int_distribution<int> unif_4(0, 2);
    // Mersenne Twister: Good quality random number generator
    std::mt19937 re(given_seed);

    while (!round_finished)
    {
        round_finished = true;
        for (int drone = 0; drone < this->n_drones; drone++)
        {
            int next_step = -1;
            int last_step = curr_sol[drone][curr_sol[drone].size() - 1];
            double min_dist = 0;

            std::vector<std::pair<int, double>> temp_neigh = {{-1, 0}, {-1, 0}, {-1, 0}};

            for (auto const &node : graph_vertices)
            {
                std::unordered_set<int> how_many_light_covered = {};
                for (int element : neigh[node])
                {
                    if (light_covered_vertices.count(element) > 0)
                    {
                        how_many_light_covered.insert(element);
                    }
                }
                if (light_covered_vertices.count(node) > 0)
                    how_many_light_covered.insert(node);

                double ratio_num = how_many_light_covered.size();

                double ratio_den = this->G.dist(last_step, node) + this->G.get_weak_coverage_node(node) * how_many_light_covered.size();

                double budget_to_use = this->G.dist(last_step, node) + this->G.dist(0, node) + this->G.get_weak_coverage_node(node) * how_many_light_covered.size();

                if ((ratio_num / ratio_den > min_dist) and budget_to_use < residual_budget[drone])
                {

                    temp_neigh[2].first = node;
                    temp_neigh[2].second = ratio_num / ratio_den;
                    min_dist = temp_neigh[2].second;

                    std::sort(temp_neigh.begin(), temp_neigh.end(), [](auto &left, auto &right)
                              { return left.second > right.second; });
                }
            }

            int coin = unif_4(re);
            assert(0 <= coin and coin < (int)temp_neigh.size());

            next_step = temp_neigh[coin].first;

            if (next_step != -1)
            {
                std::vector<int> how_many_light_covered = {};
                for (int element : neigh[next_step])
                {
                    if (light_covered_vertices.count(element) > 0)
                    {
                        how_many_light_covered.push_back(element);
                    }
                }
                how_many_light_covered.push_back(next_step);

                light_covered_vertices = set_difference(light_covered_vertices, how_many_light_covered);
                graph_vertices = set_difference(graph_vertices, how_many_light_covered);
                real_vertices_set = set_difference(real_vertices_set, how_many_light_covered);

                round_finished = false;
                curr_sol[drone].push_back(next_step);
                if (graph_vertices.find(next_step) != graph_vertices.end())
                    graph_vertices.erase(graph_vertices.find(next_step));
                residual_budget[drone] -= (this->G.dist(last_step, next_step) + this->G.get_weak_coverage_node(next_step) * how_many_light_covered.size());

                if (real_vertices_set.size() == 0)
                {
                    round_finished = true;
                }
            }
        }
    }
    return curr_sol;
}

/////////////////////////////////////////////////////
///////////////////// TOP 2 /////////////////////////
/////////////////////////////////////////////////////

std::vector<std::vector<int>> simulator_cc::top_heur(double given_budget, long given_seed)
{
    std::unordered_set<int> graph_vertices = this->G.get_vertices_set();
    std::vector<std::vector<int>> curr_sol(this->n_drones, {0});
    std::map<int, int> tree;
    graph temp_graph = graph(1, 1);

    for (auto const &i : graph_vertices)
    {
        temp_graph.add_node(i, this->G.get_coord_x(i), this->G.get_coord_y(i), this->G.get_weak_coverage_node(i), this->G.get_strong_coverage_node(i), this->G.get_profit_node(i));
    }

    std::vector<int> centers = G.get_vertices();
    std::shuffle(centers.begin(), centers.end(), std::default_random_engine(given_seed));

    centers.erase(centers.begin() + this->n_drones, centers.end());

    tree = algo::primMST(temp_graph, centers, given_budget, true);
    if (this->n_drones > 1)
        centers.erase(std::remove(centers.begin(), centers.end(), 0), centers.end());

    for (int i = 0; i < (int)centers.size(); i++)
    {
        std::vector<int> tsp_i = algo::find_TSP(temp_graph, given_budget, centers[i], tree, true);
        curr_sol[i].insert(curr_sol[i].end(), tsp_i.begin(), tsp_i.end());
        curr_sol[i].push_back(0);
        graph_vertices = set_difference(graph_vertices, curr_sol[i]);
    }

    for (int i = 0; i < (int)centers.size(); i++)
    {
        for (int j = 1; j < (int)curr_sol[i].size() - 1; j++)
        {
            double prev_cost = cost_cycle_OP(curr_sol[i], given_budget);
            int new_node = curr_sol[i][j];
            for (auto const &node : graph_vertices)
            {
                curr_sol[i][j] = node;

                if (cost_budget_sequence(this->G, curr_sol[i]) < given_budget and cost_cycle_OP(curr_sol[i], given_budget) > prev_cost)
                {
                    new_node = node;
                }
            }
            curr_sol[i][j] = new_node;
            graph_vertices = set_difference(graph_vertices, {new_node});
        }
    }
    return curr_sol;
}

std::unordered_set<int> simulator_cc::set_difference(std::unordered_set<int> main, std::vector<int> minus)
{
    for (auto const &i : minus)
    {
        if (main.find(i) != main.end())
        {
            main.erase(main.find(i));
        }
    }
    return main;
}

double simulator_cc::cost_cycle_OP(std::vector<int> _temp, double given_budget)
{
    double sum_of_elems = 0;
    for (auto const &i : _temp)
    {
        sum_of_elems += this->G.get_profit_node(i);
    }
    return sum_of_elems;
}

/////////////////////////////////////////////////////
///////////////////// SIMULATOR /////////////////////
/////////////////////////////////////////////////////

inline const char *const BoolToString(bool b)
{
    return b ? "1" : "0";
}

std::vector<std::string> simulator_cc::cycle_generation(int MAX_COUNTER, int which_alg, double budget_to_test, double radius_distance, double alpha, std::string which_heur, std::string output_data_file)
{
    int counter = 0, counter_feasible = 0;
    long seed = this->seed;
    double cycle_length = 0;
    double cycle_cost_path = 0;
    double percentage_coverage = 0;
    double percentage_coverage_max = 0;
    double min_cycle = (this->G.get_n_nodes() + 1) * 10, avg_cycle = 0, max_cycle = 0;

    double min_cycle_feasible = (this->G.get_n_nodes() + 1) * 10, avg_cycle_feasible = 0, max_cycle_feasible = 0;

    std::vector<std::string> algs = {"Greedy-TOP", "psuedoKim", "Greedy-ab", "TOP-h", "Greedy-CE", "Greedy-ab-"};
    std::vector<std::string> output_vector;
    std::vector<std::unordered_set<int>> neigh(G.get_n_nodes() + 1, std::unordered_set<int>());

    for (int i = 1; i < G.get_n_nodes(); i++)
    {
        for (int j = 1; j < G.get_n_nodes() - G.get_n_nodes_grid(); j++)
        {
            if (i != j and G.dist_eucl(i, j) <= radius_distance)
            {
                neigh[i].insert(j);
            }
        }
    }

    while (counter < MAX_COUNTER)
    {
        std::vector<std::vector<int>> sol;
        switch (which_alg)
        {
        case 0:
            sol = greedy_top(budget_to_test, seed);
            break;

        case 1:
            sol = prim_based(budget_to_test, seed);
            break;

        case 2:
            sol = greedy_alphabeta(budget_to_test, neigh, seed);
            break;

        case 3:
            sol = top_heur(budget_to_test, seed);
            break;

        case 4:
            sol = greedy_ce(budget_to_test, neigh, seed);
            break;

        case 5:
            sol = greedy_alphabeta_innerloop(budget_to_test, neigh, seed, alpha);
            break;

        default:
            exit(EXIT_FAILURE);
            break;
        }

        std::unordered_set<int> how_many_light_covered = {};

        for (size_t i = 0; i < sol.size(); i++)
        {
            double current_path_cost = cost_path_sequence(G, sol[i]);
            cycle_length += sol[i].size();
            cycle_cost_path += current_path_cost;

            std::vector<bool> temp_vector(G.get_n_nodes(), 0);
            for (size_t j = 0; j < sol[i].size(); j++)
            {
                temp_vector[sol[i][j]] = true;
            }

            std::string cycle_string = "";
            for (int j = 0; j < (int)temp_vector.size(); j++)
            {

                cycle_string.append(BoolToString(temp_vector[j]));
            }

            output_vector.push_back(cycle_string);

            Cycle temp_c;
            temp_c.budget = current_path_cost;
            temp_c.cycle = cycle_string;
            temp_c.label = {which_heur};

            if (this->output_Cycle_temp.find(cycle_string) == this->output_Cycle_temp.end())
            {
                // not found
                this->output_Cycle_temp.insert({cycle_string, temp_c});
            }
            else
            {
                // found
                std::map<std::string, Cycle>::iterator i = this->output_Cycle_temp.find(cycle_string);
                assert(i != this->output_Cycle_temp.end());
                i->second.label.insert({which_heur});
            }

            for (size_t t_node = 0; t_node < sol[i].size(); t_node++)
            {
                if (sol[i][t_node] < G.get_n_nodes() - G.get_n_nodes_grid())
                    how_many_light_covered.insert(sol[i][t_node]);

                for (int element : neigh[sol[i][t_node]])
                    if (element < G.get_n_nodes() - G.get_n_nodes_grid())
                        how_many_light_covered.insert(element);
            }
        }

        double current_sol_profit = 0;
        for (size_t z = 0; z < 4; z++)
        {
            double temp_current_sol_profit = 0;
            std::unordered_set<int> light_covered_vertices = this->G.get_real_vertices_set();
            for (size_t i = 0; i < sol.size(); i++)
            {
                temp_current_sol_profit += profit_path_sequence(G, sol[(i + z) % 4], light_covered_vertices, neigh, this->budget);
            }
            current_sol_profit = std::max(current_sol_profit, temp_current_sol_profit);
        }

        if (current_sol_profit < min_cycle)
        {
            min_cycle = current_sol_profit;
        }
        if (current_sol_profit > max_cycle)
        {
            max_cycle = current_sol_profit;
        }
        avg_cycle += current_sol_profit;

        percentage_coverage += how_many_light_covered.size();
        if (how_many_light_covered.size() > percentage_coverage_max)
            percentage_coverage_max = how_many_light_covered.size();
        if ((int)how_many_light_covered.size() == this->G.get_n_nodes() - this->G.get_n_nodes_grid())
        {
            if (current_sol_profit < min_cycle_feasible)
            {
                min_cycle_feasible = current_sol_profit;
            }
            if (current_sol_profit > max_cycle_feasible)
            {
                max_cycle_feasible = current_sol_profit;
            }
            avg_cycle_feasible += current_sol_profit;
            counter_feasible++;
        }

        counter++;
        seed++;
    }

    std::cout << algs[which_alg];
    if (which_alg == 5)
        std::cout << alpha;
    if (budget_to_test != 2250000)
    {
        std::cout << " B: " << budget_to_test;
    }

    std::cout << " #iterazioni: " << counter << " | #cicli: " << output_vector.size() << " | avg #nodes in path: " << cycle_length / output_vector.size() << " | avg used budget (path only): " << cycle_cost_path / output_vector.size();

    std::cout << " | p_min: " << min_cycle << " p_avg: " << avg_cycle / counter << " p_max: " << max_cycle;
    percentage_coverage = std::ceil(((percentage_coverage / counter) / (this->G.get_n_nodes() - this->G.get_n_nodes_grid())) * 100.0) / 100.0;
    percentage_coverage_max = std::ceil(((percentage_coverage_max) / (this->G.get_n_nodes() - this->G.get_n_nodes_grid())) * 100.0) / 100.0;
    std::cout << " | n_feasible: " << counter_feasible << " (avg: " << percentage_coverage * 100 << "% - max: " << percentage_coverage_max * 100 << "%) | p_min_f: " << min_cycle_feasible << " p_avg_f: " << avg_cycle_feasible / counter_feasible << " p_max_f: " << max_cycle_feasible;

    std::ofstream outfile;
    outfile.open(output_data_file, std::ios_base::app);

    outfile << algs[which_alg];
    if (which_alg == 5)
        outfile << alpha;
    outfile << " B: " << budget_to_test;

    outfile << " #iterazioni: " << counter << " | #cicli: " << output_vector.size() << " | avg #nodes in path: " << cycle_length / output_vector.size() << " | avg used budget (path only): " << cycle_cost_path / output_vector.size();

    outfile << " | p_min: " << min_cycle << " p_avg: " << avg_cycle / counter << " p_max: " << max_cycle;
    outfile << " | n_feasible: " << counter_feasible << " (avg: " << percentage_coverage * 100 << "% - max: " << percentage_coverage_max * 100 << "%) | p_min_f: " << min_cycle_feasible << " p_avg_f: " << avg_cycle_feasible / counter_feasible << " p_max_f: " << max_cycle_feasible;

    outfile.close();

    return output_vector;
}

void simulator_cc::print_output_hashmap_to_file(std::string path_to_file)
{
    int id = 0;
    std::ofstream outfile;
    outfile.open(path_to_file, std::ios_base::app);

    outfile << "I " << this->G.get_n_nodes() << std::endl;
    outfile << "U " << this->n_drones << std::endl;
    outfile << "K " << this->output_Cycle_temp.size() << std::endl;
    outfile << "B " << this->budget << std::endl;

    outfile << std::fixed;
    outfile.precision(2);
    for (std::map<std::string, Cycle>::const_iterator it = this->output_Cycle_temp.begin();
         it != this->output_Cycle_temp.end(); ++it)
    {
        outfile << id++ << " ";
        for (int i = 0; i < (int)it->first.size(); i++)
        {
            outfile << it->first[i] << " ";
        }
        double value = it->second.budget;
        value = std::ceil(value * 100.0) / 100.0;
        outfile << value << " ";
        std::copy(it->second.label.begin(),
                  it->second.label.end(),
                  std::ostream_iterator<std::string>(outfile, ""));
        outfile << std::endl;
    }
    outfile.close();
}

int simulator_cc::get_output_hashmap_size()
{
    return this->output_Cycle_temp.size();
}
