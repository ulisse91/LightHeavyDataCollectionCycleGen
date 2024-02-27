#include "user_input.h"

input userinput::read_user_input(int argc, char **argv)
{
    std::map<std::string, int> experiment_map =
        {{"generate-trips", 1}};

    time_t my_time = time(NULL);
    std::cout << std::endl
              << ctime(&my_time) << std::endl;

    input n_input = input();

    po::options_description desc{"Options"};
    desc.add_options()("help,h", "Help screen")("nodes,n", po::value<int>(), "[int] Number of nodes (default = 10)")("drones,q", po::value<int>(), "[int] Number of drones (default = 1)")("budget,b", po::value<double>(), "[double] Budget (default = 1)")("seed,s", po::value<int>(), "[int] seed for random graph generator")("graphfile,f", po::value<std::string>(), "[path] Graph file")("simulation", po::value<std::string>(), "[string] Which simulation/experiment to launch (default = \"generate-trips\"). Possible values: generate-trips (1)")("e", po::value<int>(), "[int] Which simulation/experiment to launch (default = 0). Possible values: see \"simulation\"")("grid-length", po::value<int>(), "[int] grid-length");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << '\n';
        exit(0);
    }
    else
    {
        if (vm.count("simulation"))
        {
            if (experiment_map.find(vm["simulation"].as<std::string>()) == experiment_map.end())
            {
                std::cout << "[ERROR] Simulation: \"" << vm["simulation"].as<std::string>() << "\" does not exist" << std::endl;
                exit(EXIT_FAILURE);
            }
            else
            {
                n_input.experiment = experiment_map[vm["simulation"].as<std::string>()];
            }
            std::cout << "Experiment: " << vm["simulation"].as<std::string>() << std::endl;
        }
        else
        {
            if (vm.count("e"))
            {
                n_input.experiment = vm["e"].as<int>();
                std::cout << "Experiment: " << vm["e"].as<int>() << std::endl;
            }
            else
            {
                std::cout << "Experiment was not set. ";
                std::cout << "Set Default value: " << n_input.experiment << std::endl;
            }
        }

        if (n_input.experiment == 1)
        {

            if (vm.count("graphfile"))
            {
                n_input.graph_file = vm["graphfile"].as<std::string>();
                std::cout << "Graph-file: " << vm["graphfile"].as<std::string>() << std::endl;
            }
            else
            {
                std::cout << "Graph-File was not set." << std::endl;
                if (vm.count("nodes"))
                {
                    n_input.n_nodes = vm["nodes"].as<int>();
                    std::cout << "Number of nodes: " << vm["nodes"].as<int>() << std::endl;
                }
                else
                {
                    std::cout << "Number of nodes was not set. ";
                    std::cout << "Set Default value: " << n_input.n_nodes << std::endl;
                }
            }
            if (vm.count("drones"))
            {
                n_input.n_drones = vm["drones"].as<int>();
                std::cout << "Number of drones: " << vm["drones"].as<int>() << std::endl;
            }
            else
            {
                std::cout << "Number of drones was not set. ";
                std::cout << "Set Default value: " << n_input.n_drones << std::endl;
            }

            if (vm.count("budget"))
            {
                n_input.budget = vm["budget"].as<double>();
                std::cout << "Budget: " << vm["budget"].as<double>() << std::endl;
            }
            else
            {
                std::cout << "Budget was not set. ";
                std::cout << "Set Default value: " << n_input.budget << std::endl;
            }
            if (vm.count("seed"))
            {
                n_input.seed = vm["seed"].as<int>();
                std::cout << "Seed: " << vm["seed"].as<int>() << std::endl;
            }
            else
            {
                std::cout << "Seed: 0 (default value)" << std::endl;
            }
        }
        else if (n_input.experiment == 2)
        {
            if (vm.count("graphfile"))
            {
                n_input.graph_file = vm["graphfile"].as<std::string>();
                std::cout << "Graph-file: " << vm["graphfile"].as<std::string>() << std::endl;
            }
            else
            {
                std::cout << "[ERROR] Graph-File was not set." << std::endl;
                std::cerr << "[ERROR] Graph-File is required for --exp 2" << std::endl;
                exit(EXIT_FAILURE);
            }
            if (vm.count("drones"))
            {
                n_input.n_drones = vm["drones"].as<int>();
                std::cout << "Number of drones: " << vm["drones"].as<int>() << std::endl;
            }
            else
            {
                std::cout << "Number of drones was not set. ";
                std::cout << "Set Default value: " << n_input.n_drones << std::endl;
            }
            if (vm.count("grid-length"))
            {
                n_input.grid_length = vm["grid-length"].as<int>();
                std::cout << "Number of grid-length: " << vm["grid-length"].as<int>() << std::endl;
            }
            else
            {
                std::cout << "grid-length was not set. ";
                std::cout << "Set Default value: " << n_input.grid_length << std::endl;
            }
            if (vm.count("budget"))
            {
                n_input.budget = vm["budget"].as<double>();
                std::cout << "Budget: " << vm["budget"].as<double>() << std::endl;
            }
            else
            {
                std::cout << "Budget was not set. ";
                std::cout << "Set Default value: " << n_input.budget << std::endl;
            }
            if (vm.count("seed"))
            {
                n_input.seed = vm["seed"].as<int>();
                std::cout << "Seed: " << vm["seed"].as<int>() << std::endl;
            }
            else
            {
                std::cout << "Seed: 0 (default value)" << std::endl;
            }
        }
        else
        {
            std::cerr << "[user_input.cpp] Something wrong. Please check!" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    std::cout << std::endl;

    return n_input;
}

void userinput::print_input(input n_input)
{
    std::cout << "n_nodes: " << n_input.n_nodes << "\n"
              << "n_drones:" << n_input.n_drones << "\n"
              << "budget: " << n_input.budget << "\n"
              << "seed: " << n_input.seed << "\n"
              << "graph_file: " << n_input.graph_file << "\n"
              << "experiment: " << n_input.experiment << "\n";
}