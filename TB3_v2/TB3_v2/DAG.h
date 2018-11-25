#pragma once
#include <iostream>
#include <vector>
#include <stack>
#include <iterator>

using namespace std;

class DAG {
private:
	vector<bool>marked;
	vector<int> edgeTo;
	vector<bool> onStack;
	stack<int> Cycle;

public:

	DAG(vector<vector<int>> G) {
		marked.resize(G.size());
		onStack.resize(G.size());
		edgeTo.resize(G.size());

		for (int i = 0; i < G.size(); i++) {
			if (!marked[i] && Cycle.empty()) {
				DFS(G, i);
			}
		}
	}

	vector<int> GetNeighbors(vector<vector<int>>G, int v) {
		vector<int> Childs;
		for (int i = 0; i < G.size(); i++) {
			if (G[i][v] != 0) {
				Childs.push_back(i);
			}
		}

		vector<int> Nums(Childs.size());
		for (int i = 0; i < Childs.size(); i++) {
			Nums[i] = Childs[i];
		}

		return Nums;
	}

	void DFS(vector<vector<int>> G, int v) {
		onStack[v] = true;
		marked[v] = true;
		vector<int> neighbour = GetNeighbors(G, v);

		for (int w : neighbour) {
			if (!Cycle.empty()) {
				return;
			}
			else if (!marked[w]) {
				edgeTo[w] = v;
				DFS(G, w);
			}
			else if (onStack[w]) {
				for (int x = v; x != w; x = edgeTo[x]) {
					Cycle.push(x);
				}
				Cycle.push(w);
				Cycle.push(v);
			}
		}
		onStack[v] = false;
	}

	bool HasCycle() {
		if (Cycle.empty()) {
			return false;
		}
		else {
			return true;
		}
	}
};