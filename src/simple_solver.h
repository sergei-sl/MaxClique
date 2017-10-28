#pragma once

#include "graph.h"
#include <vector>
class SimpleSolver
{
public:
    size_t solve(const Graph graph, Graph& out_solution);
private:
    void RunGreedyHeuristic();
    void findClique(Vertices& clique, const std::vector<Vertex>& P) noexcept;
    std::vector<Vertex> IntersectWithNeighbors(Vertex v, std::vector<Vertex>::const_iterator begin, std::vector<Vertex>::const_iterator end) noexcept;
    std::vector<Vertex> Lemma1(std::vector<Vertex>::const_iterator begin, std::vector<Vertex>::const_iterator end, size_t) noexcept;
    std::vector<Vertex> Lemma2(std::vector<Vertex>::const_iterator begin, std::vector<Vertex>::const_iterator end, Vertices&) noexcept;
protected:
    Vertices m_curr_clique;
   // Vertices m_p;
    std::vector<Vertex> m_p;
    Graph m_current_graph;
};
