/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2; -*- */
#include "graph.hpp"
Graph::Graph(int n)
{
    this->n = n;
    coord.reserve(n);
}
Graph::~Graph()
{
}
double Graph::kruskal(vector<int>& solution)
{
    vector<ARISTA_PONDERADA> G; //Grafo
    vector<ARISTA_PONDERADA> S; // Soluci��n
    ARISTA_PONDERADA ap;//arista ponderada
    int u, v;//vertices
    int N = solution.size();
    ARISTA arc;
    vector<int> fake(n, 0);
    double prizes = 0;
    vector<int> rank(N);
    vector<int> parent(N);
    boost::disjoint_sets<int*, int*> ds(&rank[0], &parent[0]);
    //representaci��n del grafo con un vector de aristas
    for (auto i = 0; i < N; i++)
    {
        u = solution[i];
        prizes += prize[u];
        fake[u] = i;
        ds.make_set(i);
        for (auto j = i + 1; j < N; j++)
        {
            v = solution[j];
            //if (vertexCluster[u] != vertexCluster[v]) //calculating without pre-processing method
            if (edge_adjacent[u][v] == true)
            {
                arc.first = u;
                arc.second = v;
                ap.first = cost[u][v];
                ap.second = arc;
                G.emplace_back(ap);
            }
        }
    }
    sort(G.begin(), G.end());
    double cost = 0;
    for (unsigned i = 0; i < G.size(); i++)
    {
        arc = G[i].second;
        u = arc.first;
        v = arc.second;
        if (ds.find_set(fake[u]) != ds.find_set(fake[v]))
        {
            S.emplace_back(pair<int, ARISTA>(G[i].first, arc));
            ds.union_set(fake[u], fake[v]);
            cost += G[i].first;
        }
    }
    return cost - prizes;
}
