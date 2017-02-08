#include "./stdc++.h" // include standard c++ library, this file is included in the project folder
# define INF 0x3f3f3f3f  // make a large numeric number
#include <cstdlib> // import randn() c standard library
#include <iostream>
#include <typeinfo>
#include "graph.h"

using namespace std;

// This class represents a directed graph using
// adjacency list representation
// Allocates memory for adjacency list
Graph::Graph(int V)
{
    this->V = V;
    adj = new list<iPair> [V];
}
int Graph::adjacent(int u, int v) // this function iterate the element in the graph, and if there is connection between u and v, return 1; otherwise return -1;
{
    std::list<iPair>::const_iterator iterator;
    for(iterator = adj[u].begin(); iterator != adj[u].end(); ++iterator){
        if ((*iterator).first == v){
             return 1;
            }
    }
    return -1;  
} 
 
void Graph::addEdge(int u, int v, int w) // this function
{
    adj[u].push_back(make_pair(v, w));
    adj[v].push_back(make_pair(u, w));
}

list<pair<int, int> >* Graph::getAdjList(){
    return adj;
    }

int Graph::Ve(){
    return V;
    }
    
int Graph::E(){
    int count = 0;
    for(int i=0; i<V; i++){
        std::list<iPair>::const_iterator iterator;
        for(iterator = adj[i].begin(); iterator != adj[i].end(); ++iterator){
            count++;
         }
        }
    return count;
    }    
    
list<iPair> Graph::neighbors(int x){
    return *adj;
    }
    
void Graph::deleteG(int x, int y){
    std::list<iPair>::const_iterator iterator;
    for(iterator = adj[x].begin(); iterator != adj[x].end(); ++iterator){
        if((*iterator).first == y){
            adj[x].remove(*iterator);
            }
    }
}

int  Graph::get_node_value(int x){
        return this->DST[x];
}
 
// Prints shortest paths from source to all other vertices
void Graph::shortestPath(int src)
{
    // Creating a priority queue 
    priority_queue< iPair, vector <iPair> , greater<iPair> > pq;
    vector<int> dist(V, INF);
    this->DST = dist;
    // Insert source itself in priority queue and initialize
    // its distance as 0.
    pq.push(make_pair(0, src));
    dist[src] = 0;
    
    /* Looping till priority queue becomes empty (or all distances are not finalized) */
    while (!pq.empty())
    {
        // The first vertex in pair is the minimum distance vertex, extract it from priority queue.
        // vertex label is stored in second of pair (it
        // has to be done this way to keep the vertices
        // sorted distance (distance must be first item
        // in pair)
        int u = pq.top().second;
        pq.pop();
 
        // 'i' is used to get all adjacent vertices of a vertex
        list< pair<int, int> >::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
        {
            // Get vertex label and weight of current adjacent
            // of u.
            int v = (*i).first;
            int weight = (*i).second;
 
            //  If there is shorted path to v through u.
            if (dist[v] > dist[u] + weight)
            {
                // Updating distance of v
                dist[v] = dist[u] + weight;
                pq.push(make_pair(dist[v], v));
            }
        }
    }
    // Print shortest distances stored in dist[]
    printf("Vertex   Distance from Source\n");
    for (int i = 0; i < V; ++i)
        printf(" From the source node %d to the target node %d \t The shortest distance is: \t %d\n",src, i, int(dist[i]));
}
