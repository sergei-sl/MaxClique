#include "simple_solver.h"

size_t SimpleSolver::solve(const Graph in_graph, Graph& out_solutions)
{
    auto buf = in_graph.GetVerticesOrderedByDegree();
    m_p.clear();
    m_visited.resize(buf.size(), 0);
    for (Vertex v : buf)
    {
        m_p.push_back(v);
    }

    m_current_graph = in_graph;
    RunHeuristics();
    std::vector<Vertex> init;
    init.clear();
    findClique(init, m_p);
    out_solutions.SetVertexSet(m_curr_clique);
    return m_curr_clique.size();
}

void SimpleSolver::RunHeuristics()
{
    std::vector<uint8_t> forbidden(m_current_graph.GetVertexCount(), 0);
    size_t max_times = m_current_graph.GetVertexCount() / 10;
    size_t sum = 0;
    for (size_t i = 0; i < max_times; ++i)
    {
        std::deque<Vertex> D;
        Vertices clique(0);
        for (Vertex v : m_p)
        {
            if (!forbidden[v])
                D.push_back(v);
        }
        if (m_curr_clique.size()> D.size() || forbidden.size()-sum < m_curr_clique.size())
            break;
        while (D.size())
        {
            Vertex k = D.front();
            D.pop_front();
            clique.push_back(k);
            D = IntersectWithNeighbors(k, D);
        }
        if (clique.size() > m_curr_clique.size())
        {
            m_curr_clique = clique;
        }
        for (Vertex v : clique)
        {
            sum++;
            forbidden[v] = 1;
        }
    }

}

void SimpleSolver::findClique(Vertices clique, std::deque<Vertex> P)
{
    if (clique.size() > m_curr_clique.size())
    {
        m_curr_clique = std::vector<Vertex>(clique);
    }
    if (clique.size() + P.size() > m_curr_clique.size())
    {
        while (P.size())
        {
            auto loc_clique = clique;
            Vertex p = P.front();
            P.pop_front();
            loc_clique.push_back(p);
            if (P.size() > 20)
            {
                Lemma1(P, loc_clique.size());
                if (P.size() > 10)
                {
                    Lemma2(P, loc_clique);
                }
            }
            findClique(loc_clique, IntersectWithNeighbors(p, P));
        }
    }

}

std::deque<Vertex> SimpleSolver::IntersectWithNeighbors(Vertex v, std::deque<Vertex>& in_P)
{
    std::deque<Vertex> out_P;
    auto adjacent = m_current_graph.GetAdjacencyMatrixRow(v);
    size_t s = 0;
    for (Vertex i : in_P)
    {
        if (i != v && adjacent[i])
        {
            out_P.push_back(i);
        }
    }
    return out_P;
}

void SimpleSolver::Lemma1(std::deque<Vertex>& in_P, size_t C_s)
{
    for (int j = 1; j <= 3; ++j)
    {
        std::deque<Vertex> out_P;
        for (Vertex i : in_P)
        {
            auto adjacent = m_current_graph.GetAdjacencyMatrixRow(i);
            size_t s = 0;
            for (Vertex v : in_P)
                s += adjacent[v];
            if (s + C_s >= m_curr_clique.size())
            {
                out_P.push_back(i);
            }
        }
        in_P.swap(out_P);
    }
}

void SimpleSolver::Lemma2(std::deque<Vertex>& in_P, Vertices& clique)
{
    for (int j = 1; j <= 3; ++j)
    {
        size_t size = in_P.size() - 1;
        std::deque<Vertex> out_P;
        for (Vertex i : in_P)
        {
            auto adjacent = m_current_graph.GetAdjacencyMatrixRow(i);
            size_t s = 0;
            for (Vertex v : in_P)
                s += adjacent[v];
            if (s == size)
            {
                clique.push_back(i);
            }
            else
            {
                out_P.push_back(i);
            }
        }
        in_P.swap(out_P);
    }

}
