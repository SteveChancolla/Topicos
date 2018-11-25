#pragma once
#include <iostream>
#include <vector>
#include<fstream>
#include <algorithm>

using namespace std;

class Tools
{
private:

public:
	void PrintMatrix(vector<vector<int>> G) {
		for (int i = 0; i < G.size(); i++) {
			for (int j = 0; j < G[0].size(); j++) {
				cout << G[i][j];
			}
			cout << "\n";
		}
	};

	void PrintMatrix1(vector<vector<string>> G) {
		for (int i = 0; i < G.size(); i++) {
			for (int j = 0; j < G[0].size(); j++) {
				cout << G[i][j] << " \t ";
			}
			cout << "\n";
		}
	};
	void PrintMatrix2(vector<vector<double>> G) {
		for (int i = 0; i < G.size(); i++) {
			for (int j = 0; j < G[0].size(); j++) {
				cout << G[i][j] << " \t ";
			}
			cout << "\n";
		}
	};

	static bool sortcol(const vector<double>& v1, const vector<double>& v2) {
		return v1[0] > v2[0];
	}


	vector<vector<int>> MakeMatrix(vector<int> Temp, int N) {
		vector<vector<int>> Result(N, vector<int>(N));
		int k = 0;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (j == i) {
					Result[i][j] = 0;
				}
				else {
					Result[j][i] = Temp[k];
					k++;
				}
			}
		}

		return Result;
	}

	vector<vector<double>> MakeMatrix2(vector<double> Temp, double N) {
		vector<vector<double>> Result(N, vector<double>(N));
		double k = 0;
		for (double i = 0; i < N; i++) {
			for (double j = 0; j < N; j++) {
				if (j > i) {
					if (j == i) {
						Result[i][j] = 0;
					}
					else {
						Result[i][j] = Temp[k];
						Result[j][i] = Temp[k];
						k++;
					}
				}
			}
		}

		return Result;
	}

	vector<vector<string>> ReadingData(string FileName) {
		ifstream in(FileName);
		if (!in.is_open())
		{
			cout << "Failed to open" << endl;
		}
		else {
			string line, c;

			vector<vector<string>> Data;
			for (int i = 0; getline(in, line); i++) {
				c.clear();
				stringstream ss(line);
				vector<string> temp;
				for (int j = 0; getline(ss, c, ','); j++) {
					temp.push_back(c);
				}
				Data.push_back(temp);
			}
			return  Data;
		}
	};

	vector<vector<int>> Dataint(vector<vector<string>> Data, vector<vector<string>>varvals) {
		vector<vector<int>> DataF;
		for (int i = 0; i < Data.size(); i++) {
			vector<int> Temp;
			for (int j = 0; j < Data[i].size(); j++) {
				int pos = find(varvals[j].begin(), varvals[j].end(), Data[i][j]) - varvals[j].begin();
				Temp.push_back(pos);
			}
			DataF.push_back(Temp);
		}
		return DataF;
	}

	vector<int> Join(int i, vector<int> Padres) {
		int N = Padres.size() + 1;
		vector<int> Arr(N);
		int k = 0;
		Arr[k++] = i;
		for (int x : Padres) {
			Arr[k++] = x;
		}
		return Arr;
	};

	vector<int> GetPadres(int i, vector<vector<int>> G) {
		vector<int> Padres;
		for (int j = 0; j < G.size(); j++) {
			if (G[i][j] == 1) {
				Padres.push_back(j);
			}
		}
		return Padres;
	};

	vector<vector<int>> GetEdge(vector<vector<int>> G) {
		vector<vector<int>> Edges;
		vector<int> Edge;
		for (int i = 0; i < G.size(); i++) {
			for (int j = 0; j < G.size(); j++) {
				if (G[i][j] == 1 && i != j && i < j) {
					Edge.push_back(i);
					Edge.push_back(j);
					Edges.push_back(Edge);

				}
				Edge.clear();
			}
		}

		return Edges;
	}

	vector<vector<int>> GetVars(int N, vector<vector<int>> G) {
		vector<vector<int>> V;
		vector<int> Vars;

		for (int i = 0; i < N; i++) {
			vector<int>Padres = GetPadres(i, G);
			Vars = Join(i, Padres);
			V.push_back(Vars);
		}
		return V;
	};

	int Factorial(int Number) {
		int Result = 1;
		if (Number < 0) {
			Result = 1;
		}
		else if (Number == 0) {
			Result = 1;
		}
		else {
			for (int i = 1; i <= Number; i++) {
				Result = Result * i;
			}
		}
		return Result;
	}

	vector<string> GetComand(string comand) {
		vector<string> temp;
		string c;

		if (getline(cin, comand)) {
			stringstream ss(comand);
			for (int j = 0; getline(ss, c, ' '); j++) {
				temp.push_back(c);
			}
		}
		return temp;
	}

	vector<vector<int>>Combinationvars(vector<int>vars) {
		vector<vector<int>>combinatedvars;
		vector<int>temp;
		if (vars.size() == 3) {
			temp.push_back(vars[0]);
			temp.push_back(vars[2]);
			combinatedvars.push_back(temp);
			temp.clear();
			temp.push_back(vars[1]);
			temp.push_back(vars[2]);
			combinatedvars.push_back(temp);
			temp.clear();
			temp.push_back(vars[2]);
			combinatedvars.push_back(temp);
		}
		else {
			for (int i = 0; i < vars.size(); i++) {
				if (i != 1) {
					temp.push_back(vars[i]);
				}
			}
			combinatedvars.push_back(temp);
			temp.clear();
			for (int i = 0; i < vars.size(); i++) {
				if (i != 0) {
					temp.push_back(vars[i]);
				}
			}
			combinatedvars.push_back(temp);
			temp.clear();
			for (int i = 2; i < vars.size(); i++) {
				temp.push_back(vars[i]);
			}
			combinatedvars.push_back(temp);
		}
		return combinatedvars;
	}

	vector<vector<int>>combinationvals(vector<int>vals) {
		vector<vector<int>>combinatedvals;
		vector<int>temp;
		if (vals.size() == 3) {
			temp.push_back(vals[0]);
			temp.push_back(vals[2]);
			combinatedvals.push_back(temp);
			temp.clear();
			temp.push_back(vals[1]);
			temp.push_back(vals[2]);
			combinatedvals.push_back(temp);
			temp.clear();
			temp.push_back(vals[2]);
			combinatedvals.push_back(temp);
		}
		else {
			for (int i = 0; i < vals.size(); i++) {
				if (i != 1) {
					temp.push_back(vals[i]);
				}
			}
			combinatedvals.push_back(temp);
			temp.clear();
			for (int i = 0; i < vals.size(); i++) {
				if (i != 0) {
					temp.push_back(vals[i]);
				}
			}
			combinatedvals.push_back(temp);
			temp.clear();
			for (int i = 2; i < vals.size(); i++) {
				temp.push_back(vals[i]);
			}
			combinatedvals.push_back(temp);
		}
		return combinatedvals;

	}
};
