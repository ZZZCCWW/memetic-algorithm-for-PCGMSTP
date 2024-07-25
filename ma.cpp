#include "ma.hpp"
MA::MA(int seed, Graph* h,int pre_sizePop,int pre_iteration,int pre_sizeChampionship,double pre_alpha,double pre_beta, double pre_timelimit)
{
	g = h;
	generator.seed(seed);
	randValue = uniform_real_distribution<double>(0.0, 1.0);
	sizePop = pre_sizePop;
	iteration_Max = pre_iteration;
	timelimit = pre_timelimit;
	sizeChampionship = pre_sizeChampionship;
	alpha0=pre_alpha;
	beta=pre_beta;
	solution = vector <int>(g->cluster.size(), 0);
	solution2 = vector <int>(g->cluster.size(), 0);
}
MA::~MA()
{
}
int MA::rnd(unsigned low, unsigned high)
{
	return low + generator() % ((high + 1) - low);
}
void MA::print_sol(vector <int> &solution)
{
  for(unsigned i = 0; i < solution.size(); i++)
    cout << "solution["<<i<<"]="<<solution[i]<<";"<<endl;
  cout << endl;
}
int MA::InitialSolution(int N)
{
	for (unsigned i = 0; i < N; i++)
	{
		if (g->cluster[i].size() > sizePop_Max)
		{
			sizePop_Max = g->cluster[i].size();
		}
	}
	while(population.size()<sizePop)
	{
		for (auto i = 0; i < N; i++)
		{
			int pos = rnd(0, sizePop_Max - 1);
			solution[i] = g->cluster[i][pos% g->cluster[i].size()];
		}
		currentCost = g->kruskal(solution);
		Localsearch();
		population.insert(make_pair(currentCost,solution));
		if (bestCost > currentCost)
		{
			bestCost = currentCost;
			solutionBest = solution;
		}
	}
	return bestCost;
}
double MA::LS1(double SolCost, vector <int>& solution)
{
	float actual = SolCost, nuevo;
	for (unsigned i = 0; i < g->cluster.size(); i++)
	{
		shuffle(g->index_random[i].begin(), g->index_random[i].end(), generator);
		for (unsigned t = 0; t < g->cluster[i].size(); t++)
		{
			int j = g->index_random[i][t];
			int templocal = solution[i];
			if (solution[i] != g->cluster[i][j])
			{
				solution[i] = g->cluster[i][j];
				nuevo = g->kruskal(solution);
				if (actual > nuevo)
				{
					actual = nuevo;
					break;
				}
				else
				{
					solution[i] = templocal;
				}
			}
		}
	}
	return actual;
}
double MA::LS2(double SolCost, vector <int>& solution)
{
	int leaf_num = 0, origin = 0, minW = 0, maxW = 99999;
	vector <int> leaf = vector <int>(g->n, 0);
	vector <int> W = vector <int>(g->n, 0);
	bestEdge = vector< vector<int> >(g->n, vector<int>(g->cluster.size(), -1));
	childrenCluster = vector< vector<int> >(g->cluster.size(), vector<int>(g->cluster.size() - 1));
	childrenVertex = vector< vector<int> >(g->cluster.size(), vector<int>(g->n - 1));
	parentCluster = vector <int>(g->cluster.size(), 0);
	clusterDegree = vector <int>(g->cluster.size(), 0);
	LS2_Kruskal(solution);
	for (int i = 0; i < g->cluster.size(); i++)
	{
		for (int j = 0; j < clusterDegree[i]; j++)
			parentCluster[childrenCluster[i][j]] = i;
	}
	for (int i = 0; i < g->cluster.size(); i++)
	{
		if (clusterDegree[i] == 0)
		{
			leaf[leaf_num] = i;
			leaf_num++;
		}
		for (int j = 0; j < g->cluster[i].size(); j++)
			W[g->cluster[i][j]] -= g->prize[g->cluster[i][j]];
	}
	while (true)
	{
		for (int l = 0; l < leaf_num; l++)
		{
			int bestCost, bestVertex, vertex1, vertex2, cluster1 = leaf[l], cluster2 = parentCluster[leaf[l]];
			for (int i = 0; i < g->cluster[cluster2].size(); i++)
			{
				vertex2 = g->cluster[cluster2][i];
				bestVertex = solution[cluster1];
				bestCost = g->cost[bestVertex][vertex2] + W[bestVertex];
				for (int j = 0; j < g->cluster[cluster1].size(); j++)
				{
					vertex1 = g->cluster[cluster1][j];
					if (g->edge_adjacent[vertex1][vertex2] == 1) continue;
					if (g->cost[vertex1][vertex2] + W[vertex1] < bestCost)
					{
						bestCost = g->cost[vertex1][vertex2] + W[vertex1];
						bestVertex = vertex1;
					}
				}
				W[vertex2] += bestCost;
				bestEdge[vertex2][cluster1] = bestVertex;
			}
			clusterDegree[cluster1]--;
			clusterDegree[cluster2]--;
		}
		if (clusterDegree[0] == 0) break;
		leaf_num = 0;
		for (int i = 0; i < g->cluster.size(); i++)
		{
			if (clusterDegree[i] == 0)
			{
				leaf[leaf_num] = i;
				leaf_num++;
			}
		}
	}
	for (int i = 0; i < g->cluster[0].size(); i++)
	{
		minW = W[g->cluster[0][i]];
		if (minW < maxW) {
			maxW = minW;
			origin = g->cluster[0][i];
		}
	}
	LS2_getsolution(origin,solution);
	SolCost = g->kruskal(solution);
	return SolCost;
}
int MA::LS2_Kruskalfind(int city) {
	if (X[city] != city) X[city] = LS2_Kruskalfind(X[city]);
	return X[city];
}
void MA::LS2_Kruskalmerge(int p1, int p2) {
	int r1 = LS2_Kruskalfind(p1);
	int r2 = LS2_Kruskalfind(p2);
	if (r1 != r2) {
		if (R[r1] > R[r2]) X[r2] = r1;
		else {
			X[r1] = r2;
			if (R[r1] == R[r2]) R[r2]++;
		}
	}
}
int MA::LS2_Kruskal(vector <int>& solution)
{
	int e_num = 0, edgenum = 0;
	X = vector <int>(g->n, 0);
	R = vector <int>(g->n, 0);
	D = vector <int>(g->n, 0);
	s_edge = new S_edge[g->n * g->n];
	for (int i = 0; i < g->n; i++) {
		X[i] = i;
		R[i] = 0;
		D[i] = 0;
	}
	for (int i = 0; i < g->cluster.size(); i++) {
		for (int j = i + 1; j < g->cluster.size(); j++) {
			if (g->edge_adjacent[solution[i]][solution[j]]) {
				s_edge[e_num].start = solution[i];
				s_edge[e_num].end = solution[j];
				s_edge[e_num].dis = g->cost[solution[i]][solution[j]];
				e_num++;
			}
		}
	}
	sort(s_edge, s_edge + e_num);
	for (int i = 0; i < e_num; i++) {
		if (LS2_Kruskalfind(s_edge[i].start) != LS2_Kruskalfind(s_edge[i].end)) {
			LS2_Kruskalmerge(s_edge[i].start, s_edge[i].end);
			childrenVertex[g->vertexCluster[s_edge[i].start]][clusterDegree[g->vertexCluster[s_edge[i].start]]] = g->vertexCluster[s_edge[i].end];
			clusterDegree[g->vertexCluster[s_edge[i].start]]++;
			childrenVertex[g->vertexCluster[s_edge[i].end]][clusterDegree[g->vertexCluster[s_edge[i].end]]] = g->vertexCluster[s_edge[i].start];
			clusterDegree[g->vertexCluster[s_edge[i].end]]++;
			edgenum++;
			if (edgenum == g->cluster.size() - 1) break;
		}
	}
	delete[]s_edge;
	LS2_buildGlobalTree(0);
	return 0;
}
int MA::LS2_buildGlobalTree(int root)
{
	int c = 0;
	parentCluster[root] = 1;
	for (int j = 0; j < clusterDegree[root]; j++)
	{
		if (parentCluster[childrenVertex[root][j]] == 0) {
			childrenCluster[root][c] = childrenVertex[root][j];
			c++;
			LS2_buildGlobalTree(childrenVertex[root][j]);
		}
	}
	clusterDegree[root] = c;
	return 0;
}
int MA::LS2_getsolution(int head, vector <int>& solution)
{
	solution[g->vertexCluster[head]] = head;
	for (int i = 0; i < g->cluster.size(); i++) {
		if (bestEdge[head][i] != -1) {
			LS2_getsolution(bestEdge[head][i], solution);
		}
	}
	return 0;
}
int MA::Crossover(){
	int racer, num=0, pos1 = rnd(0, g->cluster.size() - 1), pos2 = rnd(0, g->cluster.size() - 1);
	while (pos2 == pos1)
		pos2 = rnd(0, g->cluster.size() - 1);
	if (pos2 < pos1) {
		int pos3 = pos2;
		pos2 = pos1;
		pos1 = pos3;
	}
	while (championship.size() < sizeChampionship) {
		racer = rnd(0, population.size() - 1);
		num = 0;
		for (it = population.begin(); it != population.end(); it++)
		{
			if (num == racer)
			{
				championship.insert(make_pair(it->first, it->second));
				break;
			}
			num++;
		}
	}
	parent1 = championship.begin()->second;
	championship.clear();
	while (championship.size() < sizeChampionship) {
		racer = rnd(0, population.size() - 1);
		num = 0;
		for (it = population.begin(); it != population.end(); it++)
		{
			if (num == racer)
			{
				championship.insert(make_pair(it->first, it->second));
				break;
			}
			num++;
		}

	}
	parent2 = championship.begin()->second;
	championship.clear();
	for (auto i = 0; i < g->cluster.size(); i++)
	{
		if (pos1 <= i && i <= pos2)
		{
			solution[i] = parent1[i];
			solution2[i] = parent2[i];
		}
		else
		{
			solution[i] = parent2[i];
			solution2[i] = parent1[i];
		}
	}
	currentCost = g->kruskal(solution);
	population1.insert(make_pair(currentCost,solution));
	if (bestCost > currentCost)
	{
		bestCost = currentCost;
		solutionBest = solution;
	}
	currentCost = g->kruskal(solution2);
	population1.insert(make_pair(currentCost,solution2));
	if (bestCost > currentCost)
	{
		bestCost = currentCost;
		solutionBest = solution2;
	}
	return 0;
}
int MA::Mutation()
{
	int racer,num,pos1=rnd(0, g->cluster.size() - 1);
	while (championship.size() < sizeChampionship) {
		racer = rnd(0, population.size() - 1);
		num = 0;
		for (it = population.begin(); it != population.end(); it++)
		{
			if (num == racer)
			{
				championship.insert(make_pair(it->first, it->second));
				break;
			}
			num++;
		}
	}
	solution = championship.begin()->second;
	solution[pos1] = g->cluster[pos1][rnd(0, sizePop_Max - 1) % g->cluster[pos1].size()];
	currentCost = g->kruskal(solution);
	population1.insert(make_pair(currentCost,solution));
	championship.clear();
	while (championship.size() < sizeChampionship) {
		racer = rnd(0, population.size() - 1);
		num = 0;
		for (it = population.begin(); it != population.end(); it++)
		{
			if (num == racer)
			{
				championship.insert(make_pair(it->first, it->second));
				break;
			}
			num++;
		}
	}
	currentCost = championship.begin()->first;
	solution = championship.begin()->second;
	pos1 = rnd(0, g->cluster.size() - 1);
	solution[pos1] = g->cluster[pos1][rnd(0, sizePop_Max - 1) % g->cluster[pos1].size()];
	currentCost = g->kruskal(solution);
	population1.insert(make_pair(currentCost,solution));
	championship.clear();
	return 0;
}
int MA::Localsearch()
{
	double rate = randValue(generator);
	if(rate<beta)
		currentCost = LS1(currentCost, solution);
	else
		currentCost = LS2(currentCost, solution);
	if (bestCost > currentCost)
	{
		bestCost = currentCost;
		solutionBest = solution;
	}
	return 0;
}
int MA::Selection()
{
	int elite_num = 0;
	for (it = population.begin(); it != population.end(); it++)
	{
		currentCost = it->first;
		solution = it->second;
		population2.insert(make_pair(currentCost,solution));
		population2.erase(--population2.end());
		elite_num++;
		if (elite_num == sizePop / 10)
			break;
	}
	population.swap(population2);
	return 0;
}
int MA::Alpha_adjust()
{
	double min_value,range=0.0;
	ASD=0,fitness=0;
	it = population.end();
	it--;
	min_value=it->first;
	range=min_value-population.begin()->first;
	if(range==0)
	{
		alpha=alpha0;
		return 0;
	}
	for (it = population.begin(); it != population.end(); it++)
		fitness+=100*(min_value-it->first)/range;
	fitness=fitness/population.size();
	for (it = population.begin(); it != population.end(); it++)
		ASD+=(100*(min_value-it->first)/range-fitness)*(100*(min_value-it->first)/range-fitness);
	ASD=sqrt(ASD)/population.size();
	alpha=alpha0/(1+ASD);
	return 0;
}
int MA::maforpcgmstp(std::chrono::time_point<std::chrono::high_resolution_clock> start,int bkv)
{
	InitialSolution(g->cluster.size());
	double lastimprove=0.0;
	for (auto iteration = 1; iteration <= iteration_Max; iteration++)
	{
		// if(chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start).count()>timelimit*1e9)
		//  	break;
		//cout<<"iter "<<iteration<<", best current solution is "<<bestCost<<" and the best known solution is "<<bkv<<endl;
		// if(bestCost==bkv) break;
		Alpha_adjust();
		while (population1.size() < sizePop)
		{
			if (randValue(generator) < alpha)
				Mutation();
			else
				Crossover();
		}
		for (it = population1.begin(); it != population1.end(); it++)
		{
				currentCost = it->first;
				solution = it->second;
				if (currentCost!=lastimprove)
				{
					lastimprove=currentCost;
					Localsearch();
				}
				population2.insert(make_pair(currentCost,solution));
		}
		population1.clear();
		Selection();
		population2.clear();
	}
	//print_sol(solutionBest);
	return bestCost;
}
