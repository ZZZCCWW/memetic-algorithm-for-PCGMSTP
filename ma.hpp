/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2; -*- */
#ifndef _MA_H_
#define _MA_H_
#include "graph.hpp"
class MA {
public:
	int iteration_Max,sizePop,sizePop_Max,sizeChampionship;
	double currentCost=99999999,bestCost = 99999999;
	double alpha,alpha0,beta,ASD,fitness,timelimit;
	Graph* g;
	default_random_engine generator;
	uniform_real_distribution<double> randValue;
	multimap<int,vector <int>> population;
	multimap<int, vector <int>> population1; //population after crossover and mutation
	multimap<int, vector <int>> population2; //population1 after localsearch
	multimap<int, vector <int>> championship;
	multimap<int, vector <int>>::iterator it;
	vector <int> parent1;
	vector <int> parent2;
	vector <int> solutionBest;
	vector <int> solution;
	vector <int> solution2;
	//LS2
	vector <int> clusterDegree;
	vector <int> parentCluster;
	vector<vector<int>> bestEdge;
	vector<vector<int>> childrenCluster;
	vector<vector<int>> childrenVertex;
	vector <int> X;
	vector <int> R;
	vector <int> D;
	struct S_edge
	{
		int start;
		int end;
		int dis;
		bool operator<(const S_edge& n)const {
			return dis < n.dis;
		}
	};
	S_edge* s_edge;
	// Functions
	MA(int seed, Graph* h,int pre_sizePop,int pre_iteration,int pre_sizeChampionship,double pre_alpha,double pre_beta, double pre_timelimit);
	virtual ~MA();
	int Crossover();
	int Mutation();
	int Selection();
	int Localsearch();
	int Alpha_adjust();
	double LS1(double SolCost, vector <int>& solution);
	double LS2(double SolCost, vector <int>& solution);
	int LS2_Kruskal(vector <int>& s_node); //Turning a solution into an MST
	int LS2_buildGlobalTree(int root); // Turning an MST into a global tree
	int LS2_getsolution(int head, vector <int>& solution); //Turning the improved MST into an improved solution
	int LS2_Kruskalfind(int city);
	void LS2_Kruskalmerge(int p1, int p2);
	void print_sol(vector <int> &solution);
	int rnd(unsigned low, unsigned high);
	int InitialSolution(int N);
	int maforpcgmstp(std::chrono::time_point<std::chrono::high_resolution_clock> start,int bkv);
};
#endif /* _MA_H_ */