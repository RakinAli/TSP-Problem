#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono> // Include the chrono library
#include <limits>
#include <stack>
#include <unordered_map>
#include <unordered_set>

// Define a point structure
struct Point
{
  double x, y;
};

// Read input from stdin
std::vector<Point> readInput()
{
  int total_nodes;
  std::cin >> total_nodes;
  std::vector<Point> nodes;

  for (int i = 0; i < total_nodes; i++)
  {
    Point p;
    std::cin >> p.x >> p.y;
    nodes.push_back(p);
  }

  return nodes;
}

// Calculate the Euclidean distance between two points
double distance(const Point &p1, const Point &p2)
{
  double dx = p1.x - p2.x;
  double dy = p1.y - p2.y;
  return std::sqrt(dx * dx + dy * dy);
}

// Calculate the total distance of a tour
double tourDistance(const std::vector<Point> &nodes, const std::vector<int> &tour)
{
  double total_distance = 0.0;
  for (int i = 0; i < tour.size(); i++)
  {
    total_distance += distance(nodes[tour[i]], nodes[tour[(i + 1) % tour.size()]]);
  }
  return total_distance;
}

// Build a distance matrix for all nodes
std::vector<std::vector<double>> buildDistanceMatrix(const std::vector<Point> &nodes)
{
  std::vector<std::vector<double>> distance_matrix(nodes.size(), std::vector<double>(nodes.size(), 0.0));

  for (int i = 0; i < nodes.size(); i++)
  {
    for (int j = 0; j < nodes.size(); j++)
    {
      distance_matrix[i][j] = distance(nodes[i], nodes[j]);
    }
  }

  return distance_matrix;
}

// Naive TSP algorithmc--> gives 3 points in kattis
std::vector<int> greedyTour(const std::vector<Point> &nodes, int n)
{
  std::vector<int> tour(n, -1);
  std::vector<bool> used(n, false);
  tour[0] = 0;
  used[0] = true;

  for (int i = 1; i < n; i++)
  {
    int best = -1;
    for (int j = 0; j < n; j++)
    {
      if (!used[j] && (best == -1 || distance(nodes[tour[i - 1]], nodes[j]) < distance(nodes[tour[i - 1]], nodes[best])))
      {
        best = j;
      }
    }
    tour[i] = best;
    used[best] = true;
  }

  return tour;
}

// Function to find the nearest unvisited node from the current node and pick one of the "neighbors" closest nodes randomly
int RandomNN(int current_node, const std::vector<bool> &visited, const std::vector<std::vector<double>> &distance_matrix, int neighbors)
{
  int node_count = distance_matrix.size();
  std::vector<int> candidates;

  // Find all unvisited nodes
  for (int i = 0; i < node_count; ++i)
  {
    if (!visited[i] && i != current_node)
    {
      candidates.push_back(i);
    }
  }

  // If there are no unvisited nodes, return to the starting node
  if (candidates.size() <= neighbors)
  {
    // If the number of candidates is less than or equal to magic_number, return one of them randomly
    return candidates[rand() % candidates.size()];
  }

  else
  {
    // Select "neighbors" closest candidates
    std::vector<std::pair<int, double>> closest_candidates;
    for (int candidate : candidates)
    {
      double distance = distance_matrix[current_node][candidate];
      closest_candidates.push_back(std::make_pair(candidate, distance));
    }

    // Sort the candidates by distance in ascending order
    std::sort(closest_candidates.begin(), closest_candidates.end(), [](const auto &a, const auto &b)
              { return a.second < b.second; });

    // Choose one of the magic_number closest candidates randomly
    int random_index = rand() % neighbors;
    return closest_candidates[random_index].first;
  }
}

// Builds a random greedy tour. Randomly picks the "magic_number" of the11 closest nodes and then builds a greedy tour from there.
std::vector<int> randomGreedyTour(const std::vector<Point> &nodes, std::vector<std::vector<double>> distance_matrix, int magic_number)
{
  int node_count = nodes.size();
  std::vector<bool> visited(node_count, false);
  std::vector<int> tour;

  // Start at a random node
  int current_node = rand() % node_count;
  tour.push_back(current_node);
  visited[current_node] = true;

  // Starting the random greedy tour
  while (tour.size() < node_count)
  {
    // Find the nearest unvisited node
    int next_node = RandomNN(current_node, visited, distance_matrix, magic_number);
    tour.push_back(next_node);
    visited[next_node] = true;
    current_node = next_node;
  }
  return tour;
}

// Print the tour
void printTour(const std::vector<int> &tour)
{
  for (int i = 0; i < tour.size(); i++)
  {
    std::cout << tour[i] << std::endl;
  }
}

// Two-opt algorithm. Returns a better tour if one is found. --> 18 points in kattis on it's own
std::vector<int> twoOpt(const std::vector<Point> &nodes, std::vector<int> tour)
{
  bool finding_tour = true;
  // Look for a better tour
  while (finding_tour)
  {
    finding_tour = false;
    // Pick two edges
    for (int i = 0; i < tour.size() - 1; i++)
    {
      // Pick the second edge
      for (int j = i + 1; j < tour.size(); j++)
      {
        //
        double d1 = distance(nodes[tour[i]], nodes[tour[i + 1]]);                     // distance between 1 and 2
        double d2 = distance(nodes[tour[j]], nodes[tour[(j + 1) % tour.size()]]);     // distance between 3 and 4
        double d3 = distance(nodes[tour[i]], nodes[tour[j]]);                         // distance between 1 and 3
        double d4 = distance(nodes[tour[i + 1]], nodes[tour[(j + 1) % tour.size()]]); // distance between 2 and 4

        // If the distance between 1-2 and 3-4 is greater than the distance between 1-3 and 2-4 --> reverse the tour
        if (d1 + d2 > d3 + d4)
        {
          std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
          finding_tour = true;
        }
      }
    }
  }
  std::vector<int> newTour = tour;
  return newTour;
}

// Repeated Randomized Nearest Neighbour with 2-opt. Returns the best tour found within a time limit. 25 points in kattis
std::vector<int> RRNN(const std::vector<Point> &nodes, std::vector<std::vector<double>> distance_matrix, int magic_number, double time_limit)
{
  int node_count = nodes.size();
  std::vector<int> best_tour;
  double best_distance = 0.0;

  // Get the current time
  auto start_time = std::chrono::high_resolution_clock::now();

  // Worst case --> 2 OPT on the greedy tour
  std::vector<int> greedy_tour = greedyTour(nodes, node_count);
  std::vector<int> new_greedy_tour = twoOpt(nodes, greedy_tour);

  best_distance = tourDistance(nodes, new_greedy_tour);
  best_tour = new_greedy_tour;

  // Repeat the algorithm until the time limit is reached
  while (true)
  {
    // Check if the time limit has been reached
    auto current_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(current_time - start_time);
    if (elapsed_time.count() >= time_limit)
    {
      break; // Time limit reached, exit the loop
    }

    // Build a random greedy tour
    std::vector<int> tour = randomGreedyTour(nodes, distance_matrix, magic_number);

    // Improve the tour with 2-opt
    std::vector<int> new_tour = twoOpt(nodes, tour);

    // Calculate the distance of the new tour
    double new_distance = tourDistance(nodes, new_tour);

    // If the new tour is better than the best tour, update the best tour
    if (best_tour.empty() || new_distance < best_distance)
    {
      best_tour = new_tour;
      best_distance = new_distance;
    }
  }
  return best_tour;
}

std::vector<std::pair<int, int>> findMinimumSpanningTree(const std::vector<std::vector<double>> &distance_matrix)
{
  int n = distance_matrix.size();
  std::vector<bool> in_mst(n, false);                                    // Tracks whether a vertex is in the MST
  std::vector<double> min_weight(n, std::numeric_limits<double>::max()); // Minimum weight to connect to MST
  std::vector<int> parent(n, -1);                                        // Tracks the parent node in MST
  std::vector<std::pair<int, int>> mst_edges;                            // Stores the result MST as edges

  // Starting with the first node
  min_weight[0] = 0;
  parent[0] = -1; // The root of MST has no parent

  for (int count = 0; count < n - 1; count++)
  {
    // Find the vertex with minimum weight edge that's not already in the MST
    double min = std::numeric_limits<double>::max();
    int min_index = -1;

    for (int v = 0; v < n; v++)
    {
      if (!in_mst[v] && min_weight[v] < min)
      {
        min = min_weight[v];
        min_index = v;
      }
    }

    // Include this vertex in the MST
    in_mst[min_index] = true;

    // Update the min_weight and parent index of the adjacent vertices
    for (int v = 0; v < n; v++)
    {
      if (distance_matrix[min_index][v] && !in_mst[v] && distance_matrix[min_index][v] < min_weight[v])
      {
        parent[v] = min_index;
        min_weight[v] = distance_matrix[min_index][v];
      }
    }
  }

  // Construct the MST edges from the parent array
  for (int i = 1; i < n; ++i)
  {
    if (parent[i] != -1)
    {
      mst_edges.push_back({parent[i], i});
    }
  }
  return mst_edges;
}

// Finds all the odds degree vertices in the MST
std::vector<int> findOddDegreeVertices(const std::vector<std::pair<int, int>> &mst_edges, int num_nodes)
{
  // Degree array to count the degree of each vertex
  std::vector<int> degree(num_nodes, 0);

  // Count the degree of each vertex
  for (const auto &edge : mst_edges)
  {
    degree[edge.first]++;
    degree[edge.second]++;
  }

  // Find vertices with odd degree
  std::vector<int> odd_degree_vertices;
  for (int i = 0; i < num_nodes; ++i)
  {
    if (degree[i] % 2 != 0)
    {
      odd_degree_vertices.push_back(i);
    }
  }

  return odd_degree_vertices;
}

// Finds the minimum weight matching on the subgraph induced by the odd degree vertices
std::vector<std::pair<int, int>> findMinimumWeightMatching(const std::vector<int> &odd_degree_vertices, const std::vector<std::vector<double>> &distance_matrix)
{
  std::vector<bool> matched(odd_degree_vertices.size(), false); // Keep track of matched vertices
  std::vector<std::pair<int, int>> matching;                    // Store pairs of matched vertices

  // Perform greedy matching
  for (int i = 0; i < odd_degree_vertices.size(); ++i)
  {
    if (!matched[i])
    {
      double min_distance = std::numeric_limits<double>::max();
      int min_index = -1;

      // Find the closest unmatched vertex
      for (int j = 0; j < odd_degree_vertices.size(); ++j)
      {
        if (!matched[j] && i != j)
        {
          double weight = distance_matrix[odd_degree_vertices[i]][odd_degree_vertices[j]];
          if (weight < min_distance)
          {
            min_distance = weight;
            min_index = j;
          }
        }
      }

      // If a match is found, mark both vertices as matched and add them to the matching
      if (min_index != -1)
      {
        matching.emplace_back(odd_degree_vertices[i], odd_degree_vertices[min_index]);
        matched[i] = true;
        matched[min_index] = true;
      }
    }
  }
  return matching;
}


// Add the MST edges to the Eulerian graph
std::vector<std::pair<int, int>> combineMSTAndMatching(const std::vector<std::pair<int, int>> &mst_edges, const std::vector<std::pair<int, int>> &matching)
{
  std::vector<std::pair<int, int>> eulerian_graph;

  // Add the MST edges to the Eulerian graph
  for (const auto &edge : mst_edges)
  {
    eulerian_graph.push_back(edge);
  }

  // Add the matching edges to the Eulerian graph
  for (const auto &edge : matching)
  {
    eulerian_graph.push_back(edge);
  }

  return eulerian_graph;
}

// Find an Eulerian tour in the Eulerian graph
std::vector<int> findEulerianTour(const std::vector<std::pair<int, int>> &edges, int num_nodes)
{
  std::vector<int> eulerian_tour;
  std::unordered_map<int, std::vector<int>> adj_list;
  std::stack<int> current_path;
  std::vector<int> circuit;

  // Create the adjacency list
  for (const auto &edge : edges)
  {
    adj_list[edge.first].push_back(edge.second);
    adj_list[edge.second].push_back(edge.first);
  }

  // Ensure all vertices have even degree
  for (auto &pair : adj_list)
  {
    if (pair.second.size() % 2 != 0)
    {
      std::cerr << "Graph is not Eulerian: all vertices must have even degree." << std::endl;
      return {};
    }

    // Sort the adjacency list to get the same result consistently
    std::sort(pair.second.begin(), pair.second.end());
  }

  int current_vertex = edges.begin()->first; // Start from the first vertex
  current_path.push(current_vertex);

  while (!current_path.empty())
  {
    if (!adj_list[current_vertex].empty())
    {
      // If the current vertex has neighbors, push it onto the stack and move to a neighbor
      current_path.push(current_vertex);
      int next_vertex = adj_list[current_vertex].back();
      adj_list[current_vertex].pop_back();

      // Remove the edge in the opposite direction
      adj_list[next_vertex].erase(std::find(adj_list[next_vertex].begin(), adj_list[next_vertex].end(), current_vertex));
      current_vertex = next_vertex;
    }
    else
    {
      // If the current vertex has no neighbors, add it to the circuit and pop back to the previous vertex
      circuit.push_back(current_vertex);
      current_vertex = current_path.top();
      current_path.pop();
    }
  }

  // Reverse the circuit to get the Eulerian tour
  std::reverse(circuit.begin(), circuit.end());
  return circuit;
}

// Shortcut the Eulerian tour to get a Hamiltonian circuit
std::vector<int> shortcutEulerianTour(const std::vector<int> &eulerian_tour)
{
  std::vector<int> hamiltonian_circuit;
  std::unordered_set<int> visited;

  for (int vertex : eulerian_tour)
  {
    // If we have not visited this vertex before, add it to the Hamiltonian circuit
    if (visited.insert(vertex).second)
    {
      hamiltonian_circuit.push_back(vertex);
    }
  }

  // Add the starting vertex to the end to form a circuit
  hamiltonian_circuit.push_back(hamiltonian_circuit.front());

  return hamiltonian_circuit;
}

std::vector<int> christopides(const std::vector<Point> &nodes, std::vector<std::vector<double>> distance_matrix)
{
  // 1. Find the minimum spanning tree
  std::vector<std::pair<int, int>> mst_edges = findMinimumSpanningTree(distance_matrix);

  // 2. Find all vertices with odd degree in the MST. This is done naively and can be improved
  std::vector<int> odd_degree_vertices = findOddDegreeVertices(mst_edges, nodes.size());

  // 3. Find minimum weight perfect matching on the subgraph induced by the odd degree vertices
  std::vector<std::pair<int, int>> matching = findMinimumWeightMatching(odd_degree_vertices, distance_matrix);

  // 4. Combine the edges of the MST and the matching to form a Eulerian graph
  std::vector<std::pair<int, int>> eulerian_graph = combineMSTAndMatching(mst_edges, matching);

  // 5. Find an Eulerian tour in the Eulerian graph
  std::vector<int> eulerian_tour = findEulerianTour(eulerian_graph, nodes.size());
  // Print the eulerian tour

  // 6. Convert Eulerian tour to Hamiltonian circuit (shortcutting)
  std::vector<int> hamiltonian_circuit = shortcutEulerianTour(eulerian_tour);

  return hamiltonian_circuit;
}

int main(void)
{
  // Read input from stdin and build a distance matrix
  std::vector<Point> nodes = readInput();
  std::vector<std::vector<double>> distance_matrix = buildDistanceMatrix(nodes);

  // Do Christofides
  std::vector<int> tour2 = christopides(nodes, distance_matrix);

  // 2-opt on the Christofides tour
  std::vector<int> new_tour2 = twoOpt(nodes, tour2);

  // RRNN
  std::vector<int> tour3 = RRNN(nodes, distance_matrix, 3, 1.75);

  // Output the tour
  std::cout << "Results for Christofides " << std::endl;
  std::cout << tourDistance(nodes, tour2) << std::endl;

  std::cout << "Results for RRNN " << std::endl;
  std::cout << tourDistance(nodes, tour3) << std::endl;

  std::cout << "Results for 2-opt with Christofides " << std::endl;
  std::cout << tourDistance(nodes, new_tour2) << std::endl;
}
