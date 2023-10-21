#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

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

// Builds a random greedy tour. Randomly picks the "magic_number" of the closest nodes and then builds a greedy tour from there.
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

// Repeated Randomized Nearest Neighbour with 2-opt. Returns the best tour found.
std::vector<int> RRNN(const std::vector<Point> &nodes, std::vector<std::vector<double>> distance_matrix, int magic_number, int iterations)
{
  int node_count = nodes.size();
  std::vector<int> best_tour;
  double best_distance = 0.0;

  // Repeat the algorithm for "iterations" times
  for (int i = 0; i < iterations; i++)
  {
    // Build a random greedy tour
    std::vector<int> tour = randomGreedyTour(nodes, distance_matrix, magic_number);

    // Improve the tour with 2-opt
    std::vector<int> new_tour = twoOpt(nodes, tour);

    // Calculate the distance of the new tour
    double new_distance = tourDistance(nodes, new_tour);

    // If the new tour is better than the best tour, update the best tour
    if (i == 0 || new_distance < best_distance)
    {
      best_tour = new_tour;
      best_distance = new_distance;
    }
  }
  return best_tour;
}

int main(void)
{
  // Read input from stdin and build a distance matrix
  std::vector<Point> nodes = readInput();
  std::vector<std::vector<double>> distance_matrix = buildDistanceMatrix(nodes);

  // Sanity check
  // std::cout << "Total nodes: " << nodes.size() << std::endl;

  // Do a greedy tour
  std::vector<int> tour = randomGreedyTour(nodes, distance_matrix, 3);

  // printTour(tour);

  // Improve the tour with 2-opt
  std::vector<int> tour2 = twoOpt(nodes, tour);

  // RRNN
  std::vector<int> tour3 = RRNN(nodes, distance_matrix, 3, 50);

  // Output the tour
  printTour(tour3);
}
