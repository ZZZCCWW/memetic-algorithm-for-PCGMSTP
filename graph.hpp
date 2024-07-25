#ifndef _GRAPH_H_
#define _GRAPH_H_
#include <cstdio>
#include <vector>
#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <random>
#include <bits/stdc++.h>
#include <chrono>
#include <boost/pending/disjoint_sets.hpp>
#define ARISTA pair<int, int>
#define ARISTA_PONDERADA pair<int, ARISTA>
#define Infinite 9999999.9
using namespace std;
class Graph {
public:
	int n; // vertices number
	int m; // clusters number
	vector<pair<int, int> > coord; // coordenates of edges (x, y)
	vector<vector<int> > cluster; // cluster id with its vertices
	vector<vector<bool> > edge_adjacent; // cluster id with its vertices
	vector<int> vertexCluster; // vertex i belongs to id cluster
	vector<int> prize; // prize of vertices
	vector<vector<int> > cost; // cost matrix
	vector<vector<int> > index_random;
	Graph(int n);
	virtual ~Graph();
	double kruskal(vector <int>& solution);
};
#endif /* _GRAPH_H_ */