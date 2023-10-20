#include <iostream>
#include <vector>
#include <cmath>

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
double dist(const Point &p1, const Point &p2)
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
      if (!used[j] && (best == -1 || dist(nodes[tour[i - 1]], nodes[j]) < dist(nodes[tour[i - 1]], nodes[best])))
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


int main(void)
{
  std::cout << "Hello, World!" << std::endl;

  std::vector<Point> nodes = readInput();
  std::cout << "Total nodes: " << nodes.size() << std::endl;

  // Do a greedy tour
  std::vector<int> tour = greedyTour(nodes, nodes.size());

  // Print the tour
  printTour(tour);
}
