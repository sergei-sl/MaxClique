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
    init.reserve(m_p.size());
    findClique(init.data(), init.size(), buf.data(), buf.size());
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


void SimpleSolver::findClique(Vertex* const clique, size_t clique_size, Vertex* P, size_t size) noexcept
{
    static size_t count = 0;
    if (clique_size > m_curr_clique.size())
    {
        m_curr_clique.assign(clique, clique + clique_size);
    }
    if ((clique_size + size) <= m_curr_clique.size())
        return;
    auto x = P;
    size_t it = 0;
    Vertex* buf = nullptr;
    while (it < size)
    {
        Vertex* const loc_clique = clique;
        loc_clique[clique_size] = *x;
        size_t new_size = size - it - 1;
        size_t new_clique_size = clique_size + 1;

        Vertex* buf = IntersectWithNeighbors(*x, x + 1, new_size);
        if (size - it> 15 && size < m_p.size() - 15)//Enable Domain Filtering 
        {
            buf = Lemma1(*x, buf, new_size, new_clique_size);
            if (0 == new_size)
            {
                if (buf)
                    delete buf;
                buf = nullptr;
                x++;
                it++;
                continue;
            }
            buf = Lemma2(*x, buf, new_size, clique, new_clique_size);
            if (0 == new_size)
            {
                if (buf)
                    delete buf;
                buf = nullptr;
                x++;
                it++;
                continue;
            }
            findClique(loc_clique, new_clique_size, buf, new_size);
            if (buf)
                delete buf;
            buf = nullptr;
            count++;
        }
        else
        {           
            findClique(loc_clique, new_clique_size, buf, new_size);
            if (buf)
                delete buf;
            buf = nullptr;
            count++;
        }
        x++;
        it++;
    }
}

Vertex* SimpleSolver::IntersectWithNeighbors(Vertex v, Vertex* in_P, size_t& size) noexcept
{
    Vertex* out_P = new Vertex[size];
    size_t loc_size = 0;
    const auto& adjacent = m_current_graph.GetAdjacencyMatrixRow(v);
    for (auto i = 0; i<size; ++i)
    {
        if (adjacent[in_P[i]])
        {
            out_P[loc_size] = in_P[i];
            loc_size++;
        }
    }
    size = loc_size;
    return out_P;
}

Vertex* SimpleSolver::Lemma1(Vertex v, Vertex* in_P, size_t& size, size_t C_s) noexcept
{
    Vertex* out_P = new Vertex[size];
    size_t loc_size = 0;
    for (size_t i = 0; i<size; ++i)
    {
        const auto& adjacent = m_current_graph.GetAdjacencyMatrixRow(in_P[i]);
        size_t s = 0;
        for (size_t v = 0; v<size; ++v)
            s += adjacent[in_P[v]];
        if (s + C_s >= m_curr_clique.size())
        {
            out_P[loc_size] = in_P[i];
            loc_size++;
        }
    }
    size = loc_size;
    delete in_P;
    return out_P;
}

Vertex* SimpleSolver::Lemma2(Vertex v, Vertex* in_P, size_t& size, Vertex* const clique, size_t& clique_size) noexcept
{
    Vertex* out_P = new Vertex[size];
    size_t loc_size = 0;
    for (size_t i = 0; i<size; ++i)
    {
        auto adjacent = m_current_graph.GetAdjacencyMatrixRow(in_P[i]);
        size_t s = 0;
        for (size_t v = 0; v<size; ++v)
            s += adjacent[in_P[v]];
        if (s == size)
        {
            clique[clique_size] = in_P[i];
            clique_size++;
        }
        else
        {
            out_P[loc_size] = in_P[i];
            loc_size++;
        }
    }
    size = loc_size;
    delete in_P;
    return out_P;
}