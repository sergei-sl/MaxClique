#pragma once

#include "graph.h"
#include <vector>
class SimpleSolver
{
public:
    size_t solve(const Graph graph, Graph& out_solution);
private:
    void RunGreedyHeuristic();
    std::vector<Vertex> IntersectWithNeighbors(Vertex v, std::vector<Vertex>::const_iterator begin, std::vector<Vertex>::const_iterator end) noexcept;

    void findClique(Vertex* const clique, size_t clique_size, Vertex* P, size_t size) noexcept;
    Vertex* IntersectWithNeighbors(Vertex v, Vertex* in_P, size_t& size) noexcept;
    Vertex* Lemma1(Vertex v, Vertex* in_P, size_t& size, size_t clique_size) noexcept;
    Vertex* Lemma2(Vertex v, Vertex* in_P, size_t& size, Vertex* const clique, size_t& clique_size) noexcept;
protected:
    Vertices m_curr_clique;
   // Vertices m_p;
    std::vector<Vertex> m_p;
    Graph m_current_graph;
};
