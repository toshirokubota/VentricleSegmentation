#include <StronglyConnectedComponent.h>

vector<int>
StronglyConnectedComponents(const AdjacencyMatrix& adj)
{
	vector<Vertex<int>*> vnodes(adj.size());
	for(int i=0; i<vnodes.size(); ++i)
	{
		vnodes[i] = new Vertex<int>(i);
	}
	vector<int> comp(adj.size());

	//Step 1: call DFS to compute finishing times for each vertex
	for(int i=0; i<vnodes.size(); ++i)
	{
		vnodes[i]->color = White;
	}
	int time = 0;
	for(int i=0; i<vnodes.size(); ++i)
	{
		if(vnodes[i]->color == White)
		{
			DFS(i, adj, time, vnodes);
		}
	}
	//Step 2: compute G^T by transposing the adjacency matrix
	AdjacencyMatrix adjT(adj.size());
	for(int i=0; i<adj.size(); ++i)
	{
		for(int j=0; j<adj.size(); ++j)
		{
			adjT.Set(i, j, adj.Get(j, i));
		}
	}
	//Step 3: call DFS on G^T.  Consider the vertices in order of decreasing finishing time
	for(int i=0; i<vnodes.size(); ++i)
	{
		vnodes[i]->color = White;
	}
	time = 0;
	int count = 1;
	for(int i=0; i<vnodes.size(); ++i)
	{
		//find a white vertex with the largest finishing time
		int maxF = 0;
		int index = -1;;
		for(int j=0; j<vnodes.size(); ++j)
		{
			if(vnodes[j]->color == White && vnodes[j]->f > maxF)
			{
				index = j;
				maxF = vnodes[j]->f;
			}
		}
		if(index >= 0)
		{
			vector<Vertex<int>*> trace = DFS(index, adjT, time, vnodes);
			for(int t=0; t<trace.size(); t++)
			{
				comp[trace[t]->key] = count;
			}
			count++;
		}
		else
		{
			break;
		}
	}
	for(int i=0; i<vnodes.size(); ++i)
	{
		delete vnodes[i];
	}
	return comp;
}
