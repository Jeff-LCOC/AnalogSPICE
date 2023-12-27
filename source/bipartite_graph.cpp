// Author: ShaoqunLi

#include "bipartite_graph.h"

void BipartiteGraph::addEdge(int u, int v) {
    edges[u].push_back(v);
}

void BipartiteGraph::clear() {
    for (auto& i : edges) {
        i.clear();
    }
    for (auto& i : pre) {
        i = -1;
    }
    for (auto& i : visited) {
        i = 0;
    }
    for (auto& i : result) {
        i = -1;
    }
}

bool BipartiteGraph::hungarianForPerfectMatch() {
    int cnt = 0;
    for (int i = 0; i < size; i++) {
        // reset visited vector
        for (auto& j : visited) {
            j = 0;
        }
        if (match(i)) {
            cnt++;
        }
    }
    return cnt == size;
}

// get result(u->v) from pre(v->u)
void BipartiteGraph::transPreToResult() {
    for (int i = 0; i < size; i++) {
        result[pre[i]] = i;
    }
}

// if can find an augmenting path from node u
bool BipartiteGraph::match(int u) {
    for (int i = 0; i < size; i++) {
        if (hasEdge(u, i) && !visited[i]) {
            visited[i] = 1;
            if (pre[i] == -1 || match(pre[i])) {
                pre[i] = u;
                return true;
            }
        }
    }
    return false;
}

bool BipartiteGraph::hasEdge(int u, int v) const {
    for (auto i : edges[u]) {
        if (i == v) {
            return true;
        }
    }
    return false;
}
