#include "simple_solver.h"
#include <algorithm>
#include <cassert>

bool validateColoring(const Graph& graph, const GraphColoring& coloring)
{
    for (Color c = 0; c < coloring.getColorCount(); ++c) {
        for (Vertex v1 : coloring.getVerticesByColor(c)) {
            for (Vertex v2 : coloring.getVerticesByColor(c)) {
                if (graph.HasEdge(v1, v2)) {
                    std::cout << "CHECK FAILED: EDGE IN COLOR" << c << ":" << v1 << v2;
                    assert(false);
                    return false;
                }
            }
        };
    }
    for (Vertex v = 0; v < graph.GetVertexCount(); ++v) {
        if (coloring.getVertexColor(v) >= coloring.getColorCount()) {
            std::cout << "CHECK FAILED: UNCOLORED VERTEX" << v;
            assert(false);
            return false;
        };
    }
    return true;
}

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

static bool almostEqual(double x, double y)
{
    double maxXYOne = std::max({ 1.0, std::fabs(x) , std::fabs(y) });

    return std::fabs(x - y) <= std::numeric_limits<double>::epsilon()*maxXYOne;
}

static bool isIntegerVar(double var)
{
    return almostEqual(var, 1.) || almostEqual(var, 0.);
}

size_t SimpleSolver::solve(const Graph in_graph, Graph& out_solutions)
{
    auto buf = in_graph.GetVerticesOrderedByDegree();
    m_p.clear();
    for (Vertex v : buf)
    {
        m_p.push_back(v);
    }
    min_color_num = m_p.size();
    SIZE = m_p.size();
    m_current_graph = in_graph;
    RunGreedyHeuristic();
    std::vector<Vertex> init;
    init.clear();
    init.reserve(m_p.size());
    //findClique(init.data(), init.size(), buf.data(), buf.size());
    solveLP();
    checkClique(m_curr_clique, in_graph);
    out_solutions.SetVertexSet(m_curr_clique);
    return m_curr_clique.size();
}

void SimpleSolver::addConstraint(const std::vector<double>& constr_coef)
{
    size_t count = 0;
    for (Vertex v = 0; v < constr_coef.size(); ++v)
    {
        if (almostEqual(constr_coef[v], 1.))
        {
            count++;
            rmatval.push_back(1);
            rmatind.push_back(v);
        }
    }
    rmatbeg.push_back(m_prev_offset);
    m_prev_offset += count;
}

void SimpleSolver::solveLP() noexcept
{
    std::vector<Vertex> init;
    init.clear();
    init.reserve(m_p.size());
    generateConstraints();
    if (m_curr_clique.size() == min_color_num)
        return;
    m_obj.resize(m_p.size(), 1.);
    /* Declare and allocate space for the variables and arrays where we
    will store the optimization results including the status, objective
    value, variable values, dual values, row slacks and variable
    reduced costs. */

    int      solstat;
    double   objval;
    std::vector<double> solution;
    std::vector<double> pi;
    std::vector<double> slack;
    std::vector<double> dj;
    
    int           status = 0;
    int           i, j;
    int           cur_numrows, cur_numcols;
    /* Initialize the CPLEX environment */
    env = CPXopenCPLEX(&status);

    //For debug
    //CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
    //CPXsetintparam(env, CPXPARAM_Read_DataCheck, CPX_DATACHECK_WARN);
    //
    lp = CPXcreateprob(env, &status, "lpex1");
    CPXchgobjsen(env, lp, CPX_MAX);

    CPXnewcols(env, lp, SIZE, m_obj.data(), m_lower_bounds.data(), m_upper_bounds.data(), NULL, NULL);

    senses.resize(NUMROWS);
    std::generate(senses.begin(), senses.end(), []() { return 'L'; });

    size_t NUMNZ = rmatval.size();
    CPXaddrows(env, lp, 0, NUMROWS, NUMNZ, m_main_c.data(), senses.data(), (const int*)rmatbeg.data(),
        (const int*)rmatind.data(), (const double*)rmatval.data(), NULL, NULL);

    status = CPXlpopt(env, lp);
    cur_numrows = CPXgetnumrows(env, lp);
    cur_numcols = CPXgetnumcols(env, lp);

    solution.resize(cur_numcols);
    slack.resize(cur_numrows);
    dj.resize(cur_numcols);
    pi.resize(cur_numrows);
    status = CPXgetx(env, lp, solution.data(), 0, cur_numcols - 1);
    status = CPXgetobjval(env, lp, &objval);
    if (objval <= m_curr_clique.size())
    {
        return;
    }

    Vertex for_branch;
    if (!getBranchVariable(solution, for_branch))
    {
        // all variables are integer, update clique
        m_curr_clique.clear();
        for (Vertex v : m_p)
        {
            if (almostEqual(solution[v], 1.))
            {
                m_curr_clique.push_back(v);
            }
        }
        return;
    }
    //branch = 1
    branchLP(solution, for_branch, 1.);

    //branch = 0
    branchLP(solution, for_branch, 0.);
}

void SimpleSolver::generateConstraints() noexcept
{
    m_lower_bounds.resize(m_p.size(), 0.);
    m_upper_bounds.resize(m_p.size(), 1.);
    generateEdgeConstraints();
    generateColoringConstraints();
    generateIndependentSetsConstraints();

    for each (auto& var in m_main_constraints)
    {
        NUMROWS++;
        addConstraint(var);
    }
    m_main_c.resize(m_main_constraints.size(), 1.);
}

void SimpleSolver::generateEdgeConstraints() noexcept
{
    for (Vertex i = 0; i < m_p.size(); ++i)
    {

        const auto& adjacent = m_current_graph.GetAdjacencyMatrixRow(i);
        for (Vertex j = i + 1; j < m_p.size(); ++j)
        {
            if (!adjacent[j])
            {
                std::vector<double> constr_coef;
                constr_coef.resize(m_p.size(), 0.);
                constr_coef[i] = 1;
                constr_coef[j] = 1;
                m_main_constraints.push_back(std::move(constr_coef));
            }
        }
    }
}

void SimpleSolver::generateColoringConstraints() noexcept
{
    ColoringSolver solver;
    std::vector<GraphColoring> colorings;
    solver.solve(m_current_graph, colorings);
    for (auto& col : colorings)
    {
        validateColoring(m_current_graph, col);
        size_t col_num = col.getColorCount();
        if (col_num < min_color_num)
            min_color_num = col_num;
        for (size_t color = 0; color < col_num; ++color)
        {
            const auto & vert_by_col = col.getVerticesByColor(color);
            std::vector<double> constr_coef;
            constr_coef.resize(m_p.size(), 0.);
            for (auto& vert : vert_by_col)
            {
                constr_coef[vert] = 1;
            }
            m_main_constraints.push_back(std::move(constr_coef));
        }
    }
}

void SimpleSolver::generateIndependentSetsConstraints() noexcept
{
}

bool SimpleSolver::getBranchVariable(const std::vector<double>& solution, Vertex & branch_var) noexcept
{
    double max_val = 0.0;
    bool found = false;
    branch_var = 0;
    for (Vertex& v : m_p)
    {
        if (!isIntegerVar(solution[v]) && solution[v] > max_val - std::numeric_limits<double>::epsilon())
        {
            branch_var = v;
            max_val = solution[v];
            found = true;
        }
    }
    return found;
}

void SimpleSolver::branchLP(std::vector<double>& solution, Vertex bvar, double bound) noexcept
{
    int           status = 0;
    double   objval;
    const char type_B[] = { 'B' };
    const char type_L[] = { 'L' };
    const char type_U[] = { 'U' };
    double normal_bounds[2] = { 0., 1. };
    status = CPXchgbds(env, lp, 1, (int*)&bvar, type_B, &bound);

    status = CPXlpopt(env, lp);
    CPXgetx(env, lp, solution.data(), 0, SIZE - 1);
    CPXgetobjval(env, lp, &objval);
    if (objval <= m_curr_clique.size())
    {
        return;
    }
    Vertex for_branch;
    if (!getBranchVariable(solution, for_branch))
    {
        // all variables are integer, update clique
        m_curr_clique.clear();
        for (Vertex v : m_p)
        {
            if (almostEqual(solution[v], 1.))
                m_curr_clique.push_back(v);
        }
        return;
    }
    //branch = 1
    branchLP(solution, for_branch, 1.);

    //branch = 0
    branchLP(solution, for_branch, 0.);


    status = CPXchgbds(env, lp, 1, (int*)&for_branch, type_U, normal_bounds+1);
    status = CPXchgbds(env, lp, 1, (int*)&for_branch, type_L, normal_bounds);

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