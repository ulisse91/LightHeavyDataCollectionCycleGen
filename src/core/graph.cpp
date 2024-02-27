#include "graph.h"

graph::graph(int _area_x, int _area_y)
{
    this->area_x = _area_x;
    this->area_y = _area_y;

    // DEPOT:
    this->vertices[0] = node(0, 0, 0, 0, 0, 0); // v_0 = depot
    this->n_nodes = 1;
    this->n_nodes_grid = 0;
    this->radius_distance = 0;
}

graph::graph() {}
graph::~graph() {}

void graph::create_random_graph(int number_of_nodes, double max_weight, long seed)
{
    erase_graph();
    this->vertices[0] = node(0, 0, 0, 0, 0, 0); // v_0 = depot
    this->n_nodes = 1;
    this->n_nodes_grid = 0;

    if (seed == -1)
        seed = std::random_device{}();

    std::uniform_real_distribution<double> unif_1(0, this->area_x);
    std::uniform_real_distribution<double> unif_2(0, this->area_y);
    std::uniform_int_distribution<int> unif_3(100, 300);
    std::uniform_int_distribution<int> unif_4(1, 10);

    // Mersenne Twister: Good quality random number generator
    std::mt19937 re(seed);

    for (size_t i = 0; i < (unsigned)number_of_nodes; i++)
    {
        double _x = unif_1(re);
        double _y = unif_2(re);
        double _strong_cov = 700 * unif_3(re);
        double weak_cov = 700;
        int _p = unif_4(re);

        this->add_node(_x, _y, weak_cov, _strong_cov, _p);
    }
    assert((int)this->vertices.size() == this->n_nodes);
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
double crossProduct(const node &p1, const node &p2, const node &p3)
{
    return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
}

bool compare(const node &p1, const node &p2)
{
    double cp = crossProduct(node(0, 0, 0, 0, 0, 0), p1, p2);
    if (cp == 0)
    {
        // If the points are collinear, sort by distance from the starting point
        double dist1 = (p1.x - 0) * (p1.x - 0) + (p1.y - 0) * (p1.y - 0);
        double dist2 = (p2.x - 0) * (p2.x - 0) + (p2.y - 0) * (p2.y - 0);
        return dist1 < dist2;
    }
    return cp > 0;
}

std::vector<node> graph::findConvexHull(std::map<int, node> &vertices)
{
    //   // Find the point with the smallest y-coordinate (and smallest x-coordinate if there's a tie)
    //   auto startIt = std::min_element(vertices.begin(), vertices.end(), [](const auto& p1, const auto& p2) {
    //     return (p1.second.y < p2.second.y) || (p1.second.y == p2.second.y && p1.second.x < p2.second.x);
    //   });
    node start = node(0, 0, 0, 0, 0, 0);

    // Sort the points in counterclockwise order around the starting point
    std::vector<node> sortedPoints;
    for (const auto &entry : vertices)
    {
        if (entry.second.id == start.id)
        {
            continue;
        }
        sortedPoints.push_back(entry.second);
    }
    std::sort(sortedPoints.begin(), sortedPoints.end(), compare);

    // Perform the Graham scan
    std::vector<node> convexHull;
    convexHull.push_back(start);
    convexHull.push_back(sortedPoints[0]);

    for (size_t i = 1; i < sortedPoints.size(); ++i)
    {
        while (convexHull.size() >= 2 &&
               crossProduct(convexHull[convexHull.size() - 2], convexHull[convexHull.size() - 1], sortedPoints[i]) <= 0)
        {
            convexHull.pop_back();
        }
        convexHull.push_back(sortedPoints[i]);
    }

    return convexHull;
}

bool graph::isInsidePolygon(const std::vector<node> &convexHull, double x, double y)
{
    int n = convexHull.size();
    int i = 0, j = n - 1;
    bool isInside = false;

    // Iterate through each edge of the convex hull
    while (i < n)
    {
        if (((convexHull[i].y > y) != (convexHull[j].y > y)) &&
            (x < (convexHull[j].x - convexHull[i].x) * (y - convexHull[i].y) / (convexHull[j].y - convexHull[i].y) + convexHull[i].x))
        {
            // Toggle the inside flag if the point intersects the edge
            isInside = !isInside;
        }
        j = i++;
    }

    return isInside;
}

bool graph::touch_some_node(double x, double y)
{
    for (int j = 1; j < this->n_nodes - this->n_nodes_grid; j++)
    {
        if (dist_eucl(j, x, y) <= this->radius_distance)
        {
            return true;
        }
    }
    return false;
}

void graph::add_grid(double lenght_grid)
{
    int n_nodes_original = this->n_nodes;
    // find convex hull
    std::vector<node> convexHull = findConvexHull(vertices);

    for (size_t i = 0; i < convexHull.size(); i++)
    {
        std::cout << "(" << convexHull[i].x << ", " << convexHull[i].y << "), ";
    }
    std::cout << std::endl;

    // Grid Creation
    double x = 0, y = 0;
    this->n_nodes_grid = 0;
    while (true)
    {
        double _x = x;
        double _y = y;
        double _strong_cov = 0;
        double weak_cov = 700;
        int _p = 0;

        if (isInsidePolygon(convexHull, _x, _y) and touch_some_node(_x, _y))
        {
            this->add_node(_x, _y, weak_cov, _strong_cov, _p);
            this->n_nodes_grid++;
        }
        x += lenght_grid;
        if (x > 1)
        {
            x = 0;
            y += lenght_grid;
            if (y > 1)
                break;
        }
    }

    std::cout << std::endl;
    for (size_t i = n_nodes_original; i < this->vertices.size(); i++)
    {
        std::cout << "(" << vertices[i].x << ", " << vertices[i].y << "), ";
    }
    std::cout << std::endl;
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

int graph::add_node(double _x, double _y, double weak_coverage_cost, double strong_coverage_cost, int profit)
{
    if (not check_double_node(_x, _y))
    {
        return -1;
    }

    if (_x > this->area_x or _y > this->area_y)
    {
        return -2;
    }

    this->vertices[this->n_nodes] = node(this->n_nodes, _x, _y, weak_coverage_cost, strong_coverage_cost, profit);
    this->n_nodes++;

    return 1;
}

void graph::print_graph_to_file(std::string path)
{
    std::ofstream myfile;
    myfile.open(path, std::ios::trunc);
    myfile << "id, x, y, weak_coverage_cost, strong_coverage_cost, profit" << std::endl;
    for (auto const &v : this->vertices)
        myfile << v.second.id << ", " << v.second.x << ", " << v.second.y << ", " << v.second.weak_coverage_cost << ", " << v.second.strong_coverage_cost << ", " << v.second.profit
               << std::endl;
    myfile.close();
}

int graph::add_node(int id, double _x, double _y, double weak_coverage_cost, double strong_coverage_cost, int profit)
{
    if (not check_double_node(_x, _y))
    {
        return -1;
    }
    if (_x > this->area_x or _y > this->area_y)
    {
        return -2;
    }

    this->vertices[id] = node(id, _x, _y, weak_coverage_cost, strong_coverage_cost, profit);
    this->n_nodes++;

    return 1;
}

bool graph::check_double_node(double _x, double _y)
{
    for (auto const &pair : this->vertices)
    {
        if (_x == pair.second.x and _y == pair.second.y)
        {
            return false;
        }
    }
    return true;
}

/////////////////////////////////////////////////
//////////////////// I/O ////////////////////////
/////////////////////////////////////////////////

int graph::read_graph_from_file(std::string file)
{
    std::fstream fin;
    fin.open(file, std::ios::in);

    if (fin.is_open())
    {
        int count_lines = 1;
        std::string _id, _x, _y, _weak_cov, _strong_cov, _p;

        while (fin.good())
        {
            getline(fin, _id, ',');
            getline(fin, _x, ',');
            getline(fin, _y, ',');
            getline(fin, _weak_cov, ',');
            getline(fin, _strong_cov, ',');
            getline(fin, _p);
            std::cout << _id << " " << _x << " " << _y << " " << _weak_cov << " " << _strong_cov << " " << _p << "\n";

            if (++count_lines <= 2 /* how many lines in the file to skip (starting from 1) */ or _id == "" /* skip blank lines (last line)*/)
                continue;
            this->add_node(stod(_x), stod(_y), stod(_weak_cov), stod(_strong_cov), stod(_p));
        }
        fin.close();
        return 1;
    }
    return -1; // Unable to open file
}

void graph::erase_graph()
{
    this->n_nodes = 0;
    this->n_nodes_grid = 0;
    this->vertices = std::map<int, node>();
}

/////////////////////////////////////////////////
/////////////// GETTER & SETTER /////////////////
/////////////////////////////////////////////////

int graph::get_n_nodes()
{
    return this->n_nodes;
}

int graph::get_n_nodes_grid()
{
    return this->n_nodes_grid;
}

std::vector<int> graph::get_vertices()
{
    std::vector<int> nodes_id;
    for (auto const &pair : this->vertices)
    {
        assert(pair.first == pair.second.id);
        nodes_id.push_back(pair.first);
    }
    return nodes_id;
}

std::unordered_set<int> graph::get_real_vertices_set()
{
    std::unordered_set<int> nodes_id;
    for (auto const &pair : this->vertices)
    {
        if (pair.first < this->n_nodes - this->n_nodes_grid)
        {
            assert(pair.first == pair.second.id);
            nodes_id.insert(pair.first);
        }
    }
    return nodes_id;
}

double graph::get_strong_coverage_node(int id)
{
    return this->vertices[id].strong_coverage_cost;
}

double graph::get_weight_node(int id)
{
    return get_strong_coverage_node(id);
}

double graph::get_weak_coverage_node(int id)
{
    return this->vertices[id].weak_coverage_cost;
}

void graph::set_radius_distance(double _radius_distance)
{
    this->radius_distance = _radius_distance;
}

double graph::get_profit_node(int id)
{
    return this->vertices[id].profit;
}

double graph::get_coord_x(int id)
{
    return this->vertices[id].x;
}

double graph::get_coord_y(int id)
{
    return this->vertices[id].y;
}

std::unordered_set<int> graph::get_vertices_set()
{
    std::unordered_set<int> nodes_id;
    for (auto const &pair : this->vertices)
    {
        assert(pair.first == pair.second.id);
        nodes_id.insert(pair.first);
    }
    return nodes_id;
}

/////////////////////////////////////////////////
//////////////// DISTANCES //////////////////////
/////////////////////////////////////////////////

double graph::distw(int u, int v)
{
    return dist(u, v) + this->vertices[u].strong_coverage_cost / 2 + this->vertices[v].strong_coverage_cost / 2;
}

double graph::dist(int u, int v)
{
    // 200 = 1000 meters * 0.2, i.e., 1km (square side) times 200 J/m
    return 200 * 1000 * (sqrt(pow(this->vertices[u].x - this->vertices[v].x, 2) + pow(this->vertices[u].y - this->vertices[v].y, 2)));
}

double graph::dist_eucl(int u, int v)
{
    return (sqrt(pow(this->vertices[u].x - this->vertices[v].x, 2) + pow(this->vertices[u].y - this->vertices[v].y, 2)));
}

double graph::dist_eucl(int id_node, double node_x, double node_y)
{
    return (sqrt(pow(this->vertices[id_node].x - node_x, 2) + pow(this->vertices[id_node].y - node_y, 2)));
}

std::ostream &operator<<(std::ostream &os, const graph &G)
{
    std::cout << "***** GRAPH *****" << std::endl;
    for (auto const &v : G.vertices)
        std::cout << v.second.id << ": (" << v.second.x << ", " << v.second.y << ") wc:" << v.second.weak_coverage_cost << " sc: " << v.second.strong_coverage_cost << " p:"
                  << v.second.profit << std::endl;
    std::cout << "*****************" << std::endl;

    return os;
}