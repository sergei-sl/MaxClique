#pragma once

#include "graph.h"
#include <deque>

class SimpleSolver
{
public:
    size_t solve(const Graph graph, Graph& out_solution);
private:
    void RunHeuristics();
    void findClique(Vertices clique, std::deque<Vertex> P);
    std::deque<Vertex> IntersectWithNeighbors(Vertex v, std::deque<Vertex>& P);
    void Lemma1(std::deque<Vertex>&, size_t);
    void Lemma2(std::deque<Vertex>&, Vertices&);
protected:
    Vertices m_curr_clique;
   // Vertices m_p;
    std::deque<Vertex> m_p;
    Graph m_current_graph;
    std::vector<uint8_t> m_visited;
};
