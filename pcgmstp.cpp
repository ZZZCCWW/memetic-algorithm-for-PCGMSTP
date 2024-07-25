#include "ma.hpp"
#define _CRT_SECURE_NO_WARNINGS
int n; // number of vertices
int m; // number of clusters
Graph* g;
vector< pair<string, int> > instances;
int dist(int i, int j)
{
    int xd = g->coord[i].first - g->coord[j].first;
    int yd = g->coord[i].second - g->coord[j].second;
    return (int)(sqrt(xd * xd + yd * yd) + 0.5);
}
void read(string instance)
{
    FILE* file;
    int temp, size, x, y;
    if ((file = fopen(instance.c_str(), "r")) == NULL)
        cout << "Error to read instance1: " << instance << endl;
    if (fscanf(file, "%d", &n) == EOF)
        cout << "Error to read number of vertices" << endl;
    g = new Graph(n);
    for (auto i = 0; i < n; i++)
    {
        if (fscanf(file, "%d", &x) == EOF)
            cout << "Error to read coord x" << endl;
        if (fscanf(file, "%d", &y) == EOF)
            cout << "Error to read coord y" << endl;
        g->coord.emplace_back(pair<int, int>(x, y));
    }
    if (fscanf(file, "%d", &m) == EOF)
        cout << "Error to read number of clusters" << endl;
    g->m = m;
    for (auto i = 0; i < m; i++)
    {
        if (fscanf(file, "%d", &size) == EOF)
            cout << "Error to read size of cluster " << i << endl;
        vector<int> vertex;
        for (auto j = 0; j < size; j++)
        {
            if (fscanf(file, "%d", &temp) == EOF)
                cout << "Error to read vertices of cluster" << endl;
            vertex.emplace_back(temp - 1);
        }
        g->cluster.emplace_back(vertex);
        sort(g->cluster[i].begin(), g->cluster[i].end());
    }
    for (auto i = 0; i < m; i++)
    {
        if (fscanf(file, "%d", &temp) == EOF)
            cout << "Error to read clusters" << endl;
    }
    if (fscanf(file, "%d", &temp) == EOF)
        cout << "Error to read prizes" << endl;
    if (temp == -999)
    {
        for (auto i = 0; i < n; i++)
        {
            if (fscanf(file, "%d", &temp) == EOF)
                cout << "Error to read prizes" << endl;
            g->prize.emplace_back(temp);
        }
        if (fscanf(file, "%d", &temp) == EOF)
            cout << "Error to read prizes" << endl;
        if (temp != -999)
            cout << "cannot read the prizes, second. temp:  " << temp << endl;
    }
    for (auto i = 0; i < n; i++)
        g->vertexCluster.emplace_back(0);
    for (auto i = 0; i < m; i++)
    {
        g->index_random.emplace_back(vector<int>());
        for (unsigned j = 0; j < g->cluster[i].size(); j++)
        {
            g->vertexCluster[g->cluster[i][j]] = i;
            g->index_random[i].emplace_back(j);
        }
    }
    for (auto i = 0; i < n; i++)
    {
        vector<int> temp2;
        for (auto j = 0; j < n; j++)
        {
            temp2.emplace_back(dist(i, j));
        }
        g->cost.emplace_back(temp2);
    }
    for (auto i = 0; i < n; i++)
    {
        vector<bool> temp3;
        for (auto j = 0; j < n; j++)
        {
            temp3.emplace_back(true);
        }
        g->edge_adjacent.emplace_back(temp3);
    }
    for (auto i = 0; i < m; i++)
    {
        for (unsigned j = 0; j < g->cluster[i].size(); j++)
        {
            for (unsigned k = 0; k < g->cluster[i].size(); k++)
            {
                g->edge_adjacent[g->cluster[i][j]][g->cluster[i][k]] = false;
            }
        }
    }
    bool judge = true;
    for (int  cluster1 = 0; cluster1 < m; cluster1++)
    {
        for (unsigned vertex1 = 0; vertex1 < g->cluster[cluster1].size(); vertex1++)
        {
            for (int cluster2 = cluster1 + 1; cluster2 < m; cluster2++)
            {
                for (unsigned vertex2 = 0; vertex2 < g->cluster[cluster2].size(); vertex2++)
                {
                    if (g->edge_adjacent[g->cluster[cluster1][vertex1]][g->cluster[cluster2][vertex2]] == false)
                        continue;
                    for (int cluser3 = 0; cluser3 < m; cluser3++)
                    {
                        if (g->edge_adjacent[g->cluster[cluster1][vertex1]][g->cluster[cluster2][vertex2]] == false)
                            break;
                        if (cluser3 == cluster1 || cluser3 == cluster2)
                            continue;
                        judge = false;
                        for (unsigned vertex3 = 0; vertex3 < g->cluster[cluser3].size(); vertex3++)
                        {                           
                            if (g->cost[g->cluster[cluster1][vertex1]][g->cluster[cluser3][vertex3]]>= dist(g->cluster[cluster1][vertex1], g->cluster[cluster2][vertex2])-g->prize[g->cluster[cluster2][vertex2]]|| dist(g->cluster[cluster2][vertex2], g->cluster[cluser3][vertex3])>= dist(g->cluster[cluster1][vertex1], g->cluster[cluster2][vertex2])-g->prize[g->cluster[cluster1][vertex1]])
                            {
                                    judge = true;
                                    break;
                            }
                        }
                        if (judge == false)
                            break;
                    }
                    if (judge == false) {
                        g->edge_adjacent[g->cluster[cluster1][vertex1]][g->cluster[cluster2][vertex2]] = false;
                        g->edge_adjacent[g->cluster[cluster2][vertex2]][g->cluster[cluster1][vertex1]] = false;
                    }
                }
            }
        }
    }
    fclose(file);
}
int Execution(string instances_set,int sizeInstances,int test_time,int iteration_Max,int sizePop,int sizeChampionship,
    double alpha0,double beta,double timelimit,bool cover)
{
    string f,bestvalue,instance, output;
    string path = "instances/";
    string calpath ="res/cal_"+instances_set+".txt";
    string timepath ="res/time_"+instances_set+".txt";
    int cost,bkv,seed=0;
	double total_time=0.0;
    ifstream ifs,ifs1;
    ofstream ofs,ofs1;
    if(cover)
    {
        remove(calpath.c_str());
        remove(timepath.c_str());
    }
    ifs.open("in/instances_"+instances_set+".txt", ios::in);
    ifs1.open("in/BKS_"+instances_set+".txt", ios::in);
    ofs.open("res/cal_"+instances_set+".txt", ios::app);
    ofs1.open("res/time_"+instances_set+".txt", ios::app);
	cout<<"instances_set is "<<instances_set<<" with "<<test_time<<" test time"<<endl;
    for (int j = 0; j < sizeInstances; j++)
    {
        getline(ifs, f);
        instance = f + ".GTP";
		getline(ifs1, bestvalue);
		bkv=stoi(bestvalue);
	    total_time=0;
		cout<<"running on instance "<<f<<", it's bkv is "<<bkv<<endl;
        for (int i = 0; i < test_time; i++)
        {
            std::random_device rd;
            std::mt19937 mt(rd());
            seed = rd() % 9999;
            double timeRun = 0.0;
            read(path + instance);
            MA* mma;
            mma = new MA(seed,g,sizePop,iteration_Max,sizeChampionship,alpha0,beta,timelimit);
            auto start = chrono::high_resolution_clock::now();
            ios_base::sync_with_stdio(false);
            cost=mma->maforpcgmstp(start,bkv);
            auto end = chrono::high_resolution_clock::now();
            timeRun = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
            timeRun *= 1e-9;
            cout<<cost<<endl;
            total_time+=timeRun;
            delete mma;
            delete g;
            if(i!=test_time-1)
                ofs << cost << " ";
            else
            {
                ofs << cost << endl;
                ofs1 << total_time/test_time << endl;
            }
        }
    }
    return 0;
}

int main()
{
    int test_time=10; //Repeat time of each instance
    int iteration_Max=50;
    int sizePop=500;
    double timelimit=100; //Seconds
    double sizeChampionship=5; //The size of Championship in Crossover and mutation
    double alpha0=0.125;
    double beta=0.7;
    bool cover=true;
    string instances_set="cc"; //g1,cc,u3,u5,u7,u9,126
    int sizeInstances;
    if(instances_set=="g1")
        sizeInstances=30;
    else if(instances_set=="cc")
    {
        sizeInstances=16;
    }
    else if(instances_set=="u3"||instances_set=="u5"||instances_set=="u7"||instances_set=="u9")
    {
        sizeInstances=20;
    }
    else
    {
        sizeInstances=126;
    }
    Execution(instances_set,sizeInstances,test_time,iteration_Max,sizePop,sizeChampionship,alpha0,beta,timelimit,cover);
    return 0;
}