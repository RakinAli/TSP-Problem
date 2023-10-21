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

// Naive TSP algorithm
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

// Print the tour
void printTour(const std::vector<int> &tour)
{
  for (int i = 0; i < tour.size(); i++)
  {
    std::cout << tour[i] << std::endl;
  }
}

std::vector<int> twoOpt(const std::vector<Point> &nodes, std::vector<int> tour)
{
  // Can you build the 2-opt algorithm here. It should take the tour and return a new tour.

  bool finding_tour = true;

  while (finding_tour)
  {
    finding_tour = false;
    for (int i = 0; i < tour.size() - 1; i++)
    {
      for (int j = i + 1; j < tour.size(); j++)
      {
        double d1 = distance(nodes[tour[i]], nodes[tour[i + 1]]);
        double d2 = distance(nodes[tour[j]], nodes[tour[(j + 1) % tour.size()]]);
        double d3 = distance(nodes[tour[i]], nodes[tour[j]]);
        double d4 = distance(nodes[tour[i + 1]], nodes[tour[(j + 1) % tour.size()]]);

        if (d1 + d2 > d3 + d4)
        {
          std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
          finding_tour = true;
        }
      }
    }
  }
  return tour;
}

int main(void)
{
  // std::cout << "Hello, World!" << std::endl;

  std::vector<Point> nodes = readInput();
  std::cout << "Total nodes: " << nodes.size() << std::endl;

  // Do a greedy tour
  std::vector<int> tour = greedyTour(nodes, nodes.size());

  printTour(tour);

  // std::cout << "2-opt tour: " << std::endl;

  // Print the tour
  std::vector<int> tour2 = twoOpt(nodes, tour);

  printTour(tour2);
}
