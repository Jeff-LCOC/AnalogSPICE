// Author: ShaoqunLi

#ifndef BIPARTITE_GRAPH_H_
#define BIPARTITE_GRAPH_H_

#include <vector>

class BipartiteGraph {
public:
    BipartiteGraph(int n)
        : size(n), edges(n), pre(n, -1), visited(n, 0), result(n, -1) {}
    ~BipartiteGraph() {}
    void addEdge(int u, int v);
    void clear();
    bool hungarianForPerfectMatch();
    void transPreToResult();
    std::vector<int> getResult() { return result; }
private:
    bool match(int u);
    bool hasEdge(int u, int v) const;

    const int size;
    std::vector<std::vector<int>> edges;
    std::vector<int> pre; // u->v pairs, pre means u of v
    std::vector<int> visited; // if v is visited
    std::vector<int> result; // v of each u in u->v
};

#endif // BIPARTITE_GRAPH_H_
