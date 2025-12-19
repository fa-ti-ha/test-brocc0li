#include <bits/stdc++.h>

using namespace std;


struct edge{
  int u,v;
};


edge randomized(vector<edge>&graph)
{
 ofstream fout("output.txt",ios::app);
  random_device rd;
  mt19937 gen(rd());
 uniform_int_distribution<> dist(0,graph.size()-1);

    int random=dist(gen);
    fout<<graph[random].u<<" "<<graph[random].v<<endl;
    return {graph[random].u,graph[random].v};

}
int mincut(vector<edge>&graph)
{

   edge e=randomized(graph);
    for(int i=0;i<graph.size();i++)
   {
       if(graph[i].u==e.u && graph[i].v==e.v)
       {
           graph.erase(graph.begin()+i);
       }


   }
   for(int i=0;i<graph.size();i++)
   {
       if(graph[i].u==e.v)
        graph[i].u=e.u;
       if(graph[i].v==e.v)
        graph[i].v=e.u;

   }


}




int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt",ios::app);
    edge e;
    vector<edge>graph;
//    e.u=1;
//    e.v=2;
    fin>>e.u>>e.v;
    graph.push_back(e);
//    e.u=2;
//    e.v=4;
    fin>>e.u>>e.v;
    graph.push_back(e);
//    e.u=4;
//    e.v=3;
    fin>>e.u>>e.v;
    graph.push_back(e);
//    e.u=3;
//    e.v=1;
    fin>>e.u>>e.v;
    graph.push_back(e);
//    e.u=1;
//    e.v=4;
    fin>>e.u>>e.v;
    graph.push_back(e);
    int cut=mincut(graph);


    return 0;
}
