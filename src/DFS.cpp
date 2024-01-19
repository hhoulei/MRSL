#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <bits/stdc++.h>
#include <list>
#include <stack>

using namespace std;
using namespace arma;
using namespace Rcpp;



// Class to represent a graph
class Graph {
  // No. of vertices'
  int V;

  int* orderdfs[1000];

  // Pointer to an array containing adjacency listsList
  list<int>* adj;

  // A function used by topologicalSort
  void topologicalSortUtil(int v, bool visited[],
                           stack<int>& Stack);

public:

  //int ordd[10];
  int ordd[10];
  // Constructor
  Graph(int V);

  // function to add an edge to graph
  void addEdge(int v, int w);

  // prints a Topological Sort of
  // the complete graph
  int*  topologicalSort();

  //T[] ToArray()
};

Graph::Graph(int V)
{
  this->V = V;
  adj = new list<int>[V];
}

void Graph::addEdge(int v, int w)
{
  // Add w to v’s list.
  adj[v].push_back(w);
}

// A recursive function used by topologicalSort
void Graph::topologicalSortUtil(int v, bool visited[],
                                stack<int>& Stack)
{
  // Mark the current node as visited.
  visited[v] = true;

  // Recur for all the vertices
  // adjacent to this vertex
  list<int>::iterator i;
  for (i = adj[v].begin(); i != adj[v].end(); ++i)
    if (!visited[*i])
      topologicalSortUtil(*i, visited, Stack);

    // Push current vertex to stack
    // which stores result
    Stack.push(v);
}

// The function to do Topological Sort.
// It uses recursive topologicalSortUtil()
int* Graph::topologicalSort()
{
  stack<int> Stack;

  // Mark all the vertices as not visited
  bool* visited;
  visited = new bool[V];
  for (int i = 0; i < V; i++)
    visited[i] = false;

  // Call the recursive helper function
  // to store Topological
  // Sort starting from all
  // vertices one by one

  int rt=0;
  for (int i = 0; i < V; i++)
    if (visited[i] == false)
      topologicalSortUtil(i, visited, Stack);


    // Print contents of stack
    while (Stack.empty() == false) {
      ordd[rt]=Stack.top();
      cout << Stack.top() << " ";
      rt++;
      Stack.pop();
    }

    return(ordd);

}

// Driver Code
// [[Rcpp::export]]
Rcpp::IntegerVector DFS(SEXP amattt,SEXP n_nodes)
{
  arma::mat amatt = as<arma::mat>(amattt);
  const int n=as<int>(n_nodes);

  Graph g(n);
  for(int t = 0; t < n; t++){
    for(int r=0; r < n; r++){
      if(amatt[t,r]==1) g.addEdge(t, r);
    }
  }
  /*
  // Create a graph given in the above diagram
  Graph g(6);
  g.addEdge(5, 2);
  g.addEdge(5, 0);
  g.addEdge(4, 0);
  g.addEdge(4, 1);
  g.addEdge(2, 3);
  g.addEdge(3, 1);
   */

  cout << "Following is a Topological Sort of the given "
  "graph \n";

  //vec ordd1(10);
  // Function Call
  //g.topologicalSort();
  int* ordd2;
  ordd2=g.topologicalSort();
  //cout << ordd2[1] << endl;
  //ordd1=g.topologicalSort();
  //cout << ordd0 << endl;

  /*
  int ordd3[10];
  for (int i = 0; i < 10; i++){
    ordd3[i]=*ordd2;
    ordd2++;
    cout << ordd3[i] << endl;
  }

 */


  //return ordd;
  return Rcpp::IntegerVector::create(ordd2[0],ordd2[1],ordd2[2],ordd2[3],
                                    ordd2[4],ordd2[5],ordd2[6],ordd2[7],
                                    ordd2[8],ordd2[9]);
  //return ordd2[0];
}

