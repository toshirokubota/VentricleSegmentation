#ifndef ___KRUSKAL_H___
#define ___KRUSKAL_H___

#include <vector>
#include <algorithm>
using namespace std;
#include <struct_edge.h>
#include <DisjointSet.h>

vector<struct_edge>
kruskal(vector<struct_edge>& e, int numNodes);

#endif /* ___KRUSKAL_H___ */