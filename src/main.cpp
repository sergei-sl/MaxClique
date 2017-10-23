#include <cstdio>
#include <iostream>
#include <string>
#include <map>
#include <thread>
#include <chrono>

#include "graph.h"
#include "simple_solver.h"
#include "cliquer_wrapper.h"
#include "utils.h"


void StartCheckThread(Graph graph)
{
    std::thread([=]() {
        auto result = CliquerWrapper::unweightedFindMax(graph);
        cout << "\nRight Answer:" << endl;
        cout << "Max Clique Size:"<< result.size()<<endl;
        cout << "Max Clique:";
        for (auto vertex : result)
            cout << vertex+1 << "\t";

    }).detach();

}

void StartTimeoutWatchdog(uint32_t timeout_secs)
{
    std::thread([timeout_secs]() {
        std::this_thread::sleep_for(std::chrono::seconds(timeout_secs));
        std::cout << "0" << std::endl;
        if (!timeout_secs)
            cout << "WARNING: Timeout limit was zero or was not provided!" << endl;
        exit(1);
    }).detach();
}

int main(int argc, char* argv[])
{
    uint32_t timeout = 0;
    try {
        timeout = std::stoul(argv[2]);
    }
    catch (std::invalid_argument e) {
        cout << "WARNING: Could not get timeout limit!" << endl;
        return 1;
    }
    StartTimeoutWatchdog(timeout);
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
	Graph in_graph = readDimacs(argv[1]);
    Graph out_graph(in_graph.GetVertexCount());
#ifdef DEBUG
	//StartCheckThread(in_graph);
#endif
    SimpleSolver solver;

    solver.solve(in_graph, out_graph);
    Vertices clique = out_graph.GetVertexSet();
    cout << clique.size() << endl;
    for (auto vertex : clique)
    {
        cout << vertex + 1 << " ";
    }

	return 0;
}
