#include "simple_solver.h"

void checkClique(Vertices & c, const Graph& g)
{
    for (auto v : c)
    {
        for (auto v2 : c)
        {
            if (v2 != v && !g.HasEdge(v, v2))
                std::cout << "UJAS!"<<std::endl;
        }

    }
}
size_t SimpleSolver::solve(const Graph in_graph, Graph& out_solutions)
{
    auto buf = in_graph.GetVerticesOrderedByDegree();
    m_p.clear();
    for (Vertex v : buf)
    {
        m_p.push_back(v);
    }

    m_current_graph = in_graph;
    RunGreedyHeuristic();
    std::vector<Vertex> init;
    init.clear();
    init.reserve(m_curr_clique.size() * 2);
    findClique(init, m_p);
    checkClique(m_curr_clique, in_graph);
    out_solutions.SetVertexSet(m_curr_clique);
    return m_curr_clique.size();
}

void SimpleSolver::RunGreedyHeuristic()
{
    std::vector<uint8_t> forbidden(m_current_graph.GetVertexCount(), 0);
    size_t max_times = m_current_graph.GetVertexCount() / 10;
    size_t sum = 0;
    for (size_t i = 0; i < max_times; ++i)
    {
        std::vector<Vertex> D;
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
            D.erase(D.begin());
            clique.push_back(k);
            D = IntersectWithNeighbors(k, D.begin(), D.end());
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

void SimpleSolver::findClique(Vertices& clique, const std::vector<Vertex>& P) noexcept
{
    static size_t count = 0;
    size_t g = clique.size();
    if (g > m_curr_clique.size())
    {
        m_curr_clique = clique;
    }
    if ((g + P.size()) <= m_curr_clique.size())
        return;
    auto x = P.begin();
    while (x != P.end())
    {
        auto loc_clique = clique;
        loc_clique.push_back(*x);

        if (0/*std::distance(x, P.end())> 15 && P.size() < m_p.size() - 15*/)//Disable Domain Filtering 
        {
            Vertices buf = (IntersectWithNeighbors(*x, x+1, P.end()));
            if (buf.empty())
                continue;
            buf = (Lemma1(buf.begin(), buf.end(), loc_clique.size()));
            if (buf.empty())
                return;
            buf = (Lemma2(buf.begin(), buf.end(), loc_clique));
            if (buf.empty())
                return;
            findClique(loc_clique, buf);

            count++;
        }
        else
        {
            findClique(loc_clique, IntersectWithNeighbors(*x, x+1, P.end()));

            count++;
        }
        x++;
    }
}

std::vector<Vertex> SimpleSolver::IntersectWithNeighbors(Vertex v, std::vector<Vertex>::const_iterator begin, std::vector<Vertex>::const_iterator end) noexcept
{
    std::vector<Vertex> out_P;
    out_P.reserve(std::distance(begin, end));
    const auto& adjacent = m_current_graph.GetAdjacencyMatrixRow(v);
    for (auto& i = begin; i<end;++i)
    {
        if (adjacent[*i])
        {
            out_P.push_back(*i);
        }
    }
    return out_P;
}

std::vector<Vertex> SimpleSolver::Lemma1(std::vector<Vertex>::const_iterator begin, std::vector<Vertex>::const_iterator end, size_t C_s) noexcept
{
    std::vector<Vertex> out_P;
    out_P.reserve(std::distance(begin, end));
    for (auto i = begin; i<end; ++i)
    {
        const auto& adjacent = m_current_graph.GetAdjacencyMatrixRow(*i);
        size_t s = 0;
        for (auto v = begin; v<end; ++v)
            s += adjacent[*v];
        if (s + C_s >= m_curr_clique.size())
        {
            out_P.push_back(*i);
        }
    }
    return out_P;
}

std::vector<Vertex> SimpleSolver::Lemma2(std::vector<Vertex>::const_iterator begin, std::vector<Vertex>::const_iterator end, Vertices& clique) noexcept
{
        size_t size = std::distance(begin, end) - 1;
        std::vector <Vertex> out_P;
        for (auto i = begin; i<end; ++i)
        {
            auto adjacent = m_current_graph.GetAdjacencyMatrixRow(*i);
            size_t s = 0;
            for (auto v = begin; v<end; ++v)
                s += adjacent[*v];
            if (s == size)
            {
                clique.push_back(*i);
            }
            else
            {
                out_P.push_back(*i);
            }
        }
        return out_P;
}

