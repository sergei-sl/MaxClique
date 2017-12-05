#pragma once

#include "graph.h"
#include "coloring_problem.h"
#include <vector>
#include <ilcplex/cplex.h>
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

    void addConstraint(const std::vector<double>& constr_coef);
    void solveLP() noexcept;
    void generateConstraints() noexcept;
    void generateEdgeConstraints() noexcept;
    void generateColoringConstraints() noexcept;
    void generateIndependentSetsConstraints() noexcept;
    bool getBranchVariable(const std::vector<double>&, Vertex&) noexcept;
    void branchLP(std::vector<double>& solution, Vertex to_branch, double bound) noexcept;
protected:
    Vertices m_curr_clique;
   // Vertices m_p;
    std::vector<Vertex> m_p;
    Graph m_current_graph;
    //objective function coeffitients
    std::vector<double> m_obj;

    std::vector<double> m_lower_bounds;

    std::vector<double> m_upper_bounds;

    std::vector<std::vector<double>> m_main_constraints;

    size_t min_color_num;
    std::vector<int> rmatbeg;
    std::vector<int> rmatind;
    std::vector<double> rmatval;
    std::vector<double> rhs;
    std::vector<char> senses;
    size_t SIZE, NUMROWS = 0;
    size_t m_prev_offset = 0;
    //Cplex environment
    CPXENVptr     env = NULL;
    //Cplex LP problem
    CPXLPptr      lp = NULL;

    //array of ones
    std::vector<double> m_main_c;
};
