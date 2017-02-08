#ifndef Graph_H
#define Graph_H
#include "./stdc++.h"
using namespace std;
typedef pair<int, int> iPair; // iPair ==>  Integer Pair
class Graph
{
private:
    int V;    // No. of vertices
    // In a weighted graph, we need to store vertex
    // and weight pair for every edge
    list<pair<int, int> > *adj;
    // Create a vector for distances and initialize all
    // distances as infinite (INF)
    vector<int> DST;
    
public:
    Graph(int V);  // Constructor
    // function to add an edge to graph
    void addEdge(int u, int v, int w);
    //add (G, x, y): adds to G the edge from x to y, if it is not there.
    // prints shortest path from s
    void shortestPath(int s);
    int adjacent(int u, int v);
    int Ve();
    //V (G): returns the number of vertices in the graph
    int E();
    //E (G): returns the number of edges in the graph
    list<iPair > neighbors(int x);
    //neighbors (G, x): lists all nodes y such that there is an edge from x to y.
    void deleteG(int x, int y);
    //delete (G, x, y): removes the edge from x to y, if it is there.
    int get_node_value(int x);
    //get_node_value (G, x): returns the value associated with the node x.
    void set_node_value(int x, int a);
    //set_node_value( G, x, a): sets the value associated with the node x to a.
    int get_edge_value(int x, int y);
    //get_edge_value( G, x, y): returns the value associated to the edge (x,y).
    void set_edge_value(int x, int y, int z);
    //set_edge_value (G, x, y, v): sets the value associated to the edge (x,y) to v.
    list<pair<int, int> >* getAdjList();
    // get the list value from the graph
};
 
#endif