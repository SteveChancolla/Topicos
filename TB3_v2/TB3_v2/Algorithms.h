#pragma once
#include <iostream>
#include <vector>
#include "Tools.h"
#include "Probabilities.h"
#include <algorithm>
#include "DAG.h"
#include "QMeasure.h"
#include <algorithm>

using namespace std;

class Algorithms
{
private:

public:
	Tools T;
	Probabilities P;
	QMeasure Q;

	void Marginalization(vector<vector<int>> Data, vector<int> Vars, vector<int> Card, vector<string>VarNames, vector<vector<string>> VarVals, int size, vector<int>marg) {
		int w;
		int CardAnt = 1;
		int CardNew = 1;
		int iVar;
		int iVarn;
		int k = 0;
		double newP = 0.0;
		cout << "\n";
		cout << "Factor ";
		cout << "\n";
		cout << "i" << "\t";
		int m = Vars.size();
		int v = marg.size();
		vector<int>vals(m);
		vector<int>vars1;
		vector<int>vals1;

		for (w = 0; w < Vars.size() - 1; w++) {
			cout << VarNames[Vars[w]] << "\t";
		}
		cout << VarNames[Vars[w]] << endl;
		int FactorSize = 1;
		for (int k = 0; k < Vars.size(); k++) {
			FactorSize *= Card[Vars[k]];
		}
		for (int i = 0; i < FactorSize; i++) {
			CardAnt = 1;
			for (int j = 0; j < m; j++) {
				iVar = Vars[j];
				vals[j] = (int)((int)floor(i / CardAnt) % Card[iVar]);
				CardAnt *= Card[iVar];
			}
			double Prob = 0.0;

			Prob = P.JointCalculationWD(Data, Vars, vals, Card, 1.0, false);

			cout << i << "\t";
			for (w = 0; w < vals.size() - 1; w++) {
				cout << VarVals[Vars[w]][vals[w]] + "\t";
			}
			cout << VarVals[Vars[w]][vals[w]] + "\t" << Prob;
			cout << "\n";
			k++;

			for (int m = size - 1; m < Vars.size(); m++) {
				vars1.push_back(Vars[m]);

			}

			for (int p = size - 1; p < vals.size(); p++) {
				vals1.push_back(vals[p]);

			}
			CardNew = 1;
			for (int s = 0; s < v; s++) {
				iVarn = marg[s];
				CardNew *= Card[iVarn];
			}

			newP += Prob;
			if (k%CardNew == 0) {
				cout << "Marginalizacion de ";
				for (int g = 0; g < marg.size(); g++) {
					cout << VarNames[marg[g]];
				}
				cout << " : " << newP << "\n";
				newP = 0.0;
			}
		}
	}

	void MutualInformation(vector<vector<int>>Data, vector<int>Vars, vector<int>Card) {
		int CardAnt = 1;
		int iVar;
		int m = Vars.size();
		double counter = 0.0;
		vector<int>vals(m);

		int FactorSize = 1;
		for (int k = 0; k < Vars.size(); k++) {
			FactorSize *= Card[Vars[k]];
		}
		for (int i = 0; i < FactorSize; i++) {
			CardAnt = 1;
			for (int j = 0; j < m; j++) {
				iVar = Vars[j];
				vals[j] = (int)((int)floor(i / CardAnt) % Card[iVar]);
				CardAnt *= Card[iVar];
			}

			if (Vars.size() <= 2) {
				vector<int>vars1;
				for (int l = 0; l < Vars.size(); l++) {
					if (l != 1) {
						vars1.push_back(Vars[l]);
					}
				}
				vector<int>vals1;
				for (int l = 0; l < vals.size(); l++) {
					if (l != 1) {
						vals1.push_back(vals[l]);
					}
				}

				vector<int>vars2;
				for (int l = 0; l < Vars.size(); l++) {
					if (l != 0) {
						vars2.push_back(Vars[l]);
					}
				}
				vector<int>vals2;
				for (int l = 0; l < vals.size(); l++) {
					if (l != 0) {
						vals2.push_back(vals[l]);
					}
				}

				double Prob = 0.0;

				double X = P.JointCalculationWD(Data, Vars, vals, Card, 1.0, true);
				double Y = X;
				double Z = 0.0;
				double Za = P.JointCalculationWD(Data, vars1, vals1, Card, 1.0, true);
				double Zb = P.JointCalculationWD(Data, vars2, vals2, Card, 1.0, true);
				Z = Za * Zb;
				Prob = X * log2(Y / Z);
				counter += Prob;

			}

			else {
				vector<int>vars1;
				for (int l = 2; l < Vars.size(); l++) {
					vars1.push_back(Vars[l]);
				}
				vector<int>vals1;
				for (int l = 2; l < vals.size(); l++) {
					vals1.push_back(vals[l]);
				}
				/////////////////////////////////////
				vector<int> vars2;
				for (int k = 0; k < Vars.size(); k++) {
					if (k != 1) {
						vars2.push_back(Vars[k]);
					}
				}
				vector<int> vals2;
				for (int k = 0; k < vals.size(); k++) {
					if (k != 1) {
						vals2.push_back(vals[k]);
					}
				}
				vector<int> vars21;
				for (int k = 0; k < vars2.size(); k++) {
					if (k != 0) {
						vars21.push_back(vars2[k]);
					}
				}
				vector<int> vals21;
				for (int k = 0; k < vals2.size(); k++) {
					if (k != 0) {
						vals21.push_back(vals2[k]);
					}
				}
				///////////////////////////////////////////
				vector<int> vars3;
				for (int p = 0; p < Vars.size(); p++) {
					if (p != 0) {
						vars3.push_back(Vars[p]);
					}
				}

				vector<int> vals3;
				for (int p = 0; p < vals.size(); p++) {
					if (p != 0) {
						vals3.push_back(vals[p]);
					}
				}

				vector<int> vars31;
				for (int p = 0; p < vars3.size(); p++) {
					if (p != 0) {
						vars31.push_back(vars3[p]);
					}
				}

				vector<int> vals31;
				for (int p = 0; p < vals3.size(); p++) {
					if (p != 0) {
						vals31.push_back(vals3[p]);
					}
				}


				double Prob = 0.0;

				double X = P.JointCalculationWD(Data, Vars, vals, Card, 1.0, true);
				double Y = 0.0;
				double Ya = X;
				double Yb = P.JointCalculationWD(Data, vars1, vals1, Card, 1.0, true);
				Y = Ya / Yb;
				double Z = 0.0;
				double Za1 = P.JointCalculationWD(Data, vars2, vals2, Card, 1.0, true);
				double Za2 = P.JointCalculationWD(Data, vars21, vals21, Card, 1.0, true);
				double Za = Za1 / Za2;
				double Zb1 = P.JointCalculationWD(Data, vars3, vals3, Card, 1.0, true);
				double Zb2 = P.JointCalculationWD(Data, vars31, vals31, Card, 1.0, true);
				double Zb = Zb1 / Zb2;
				Z = Za * Zb;

				Prob = X * log2(Y / Z);
				counter += Prob;
			}
		}
		std::cout << "Mutual Information score: " << counter << endl;
	}

	void Pearson(vector<vector<int>> Data, vector<int>Vars, vector<int>Card) {
		int w;
		int CardAnt = 1;
		int iVar;

		int m = Vars.size();
		double counter = 0.0;

		vector<int>vals(m);

		int FactorSize = 1;
		for (int k = 0; k < Vars.size(); k++) {
			FactorSize *= Card[Vars[k]];
		}
		for (int i = 0; i < FactorSize; i++) {
			CardAnt = 1;
			for (int j = 0; j < m; j++) {
				iVar = Vars[j];
				vals[j] = (int)((int)floor(i / CardAnt) % Card[iVar]);
				CardAnt *= Card[iVar];
			}
			if (Vars.size() <= 2) {

				double Prob = 0.0;
				double a = P.JointCountTable(Data, Vars, vals);
				double b = 0.0;
				double b1 = P.MarginalCountTable(Data, Vars[0], vals[0]);
				double b2 = P.MarginalCountTable(Data, Vars[1], vals[1]);
				double b3 = Data.size();
				b = (b1*b2) / b3;
				double c = (a - b);
				Prob = (c * c) / b;
				counter += Prob;
			}
			else {
				double Prob = 0.0;
				double a = P.JointCountTable(Data, Vars, vals);
				double b = 0.0;

				vector<int>vars1;
				for (int i = 0; i < Vars.size(); i++) {
					if (i != 1)
					{
						vars1.push_back(Vars[i]);
					}
				}

				vector<int>vals1;
				for (int i = 0; i < vals.size(); i++) {
					if (i != 1)
					{
						vals1.push_back(vals[i]);
					}
				}

				vector<int>vars2;
				for (int i = 0; i < Vars.size(); i++) {
					if (i != 0)
					{
						vars2.push_back(Vars[i]);
					}
				}
				vector<int>vals2;
				for (int i = 0; i < vals.size(); i++) {
					if (i != 1)
					{
						vals2.push_back(vals[i]);
					}
				}

				double b1 = P.JointCountTable(Data, vars1, vals1);
				double b2 = P.JointCountTable(Data, vars2, vals2);
				double b3 = Data.size();
				b = (b1*b2) / b3;
				double c = (a - b);
				Prob = (c*c) / b3;
				if (b == 0) {
					Prob = 0;
				}
				counter += Prob;
			}

		}
		cout << "Pearson	 score: " << counter << endl;

	}

	vector<int> BruteForce(int n, int NumBits) {
		vector<int> Binary(NumBits);

		for (int i = 0; i < NumBits; ++i, n /= 2) {
			switch (n % 2)
			{
			case 0:
				Binary[i] = 0;
				break;
			case 1:
				Binary[i] = 1;
				break;
			}
		}
		return Binary;
	};

	vector<vector<int>> K2(int N) {
		vector<vector<int>> Result(N, vector<int>(N));
		for (int i = 0; i < Result.size(); i++) {
			for (int j = 0; j < Result.size(); j++) {
				if (i > j) {
					Result[i][j] = 1;
					T.PrintMatrix(Result);
					cout << "\n";
				}
			}
		}
		return Result;
	}

	vector<vector<int>> ChowLiu(vector<vector<int>> Matrix, vector<int>Card, vector<vector<int>> Data, double N) {

		vector<vector<int>> edges = T.GetEdge(Matrix);
		vector<int> Vars;
		vector<double> peso;
		vector<vector<double>> Complete;
		vector<vector<int>> Result(N, vector<int>(N));
		DAG *dca;
		for (int i = 0; i < edges.size(); i++) {
			Vars = edges[i];
			int iVar;
			int CardAnt = 1;
			vector<int> Vals(Vars.size());

			int FactorSize = 1;
			for (int v = 0; v < Vars.size(); v++) {
				FactorSize *= Card[Vars[v]];
			}
			vector<double> Prob;

			for (int g = 0; g < FactorSize; g++) {
				CardAnt = 1;
				for (int k = 0; k < Vars.size(); k++) {
					iVar = Vars[k];
					Vals[k] = (int)((int)floor(g / CardAnt) % Card[iVar]);
					CardAnt *= Card[iVar];
				}

				double t = 0.0;
				double r = 0.0;
				double a = P.JointCalculationWD(Data, Vars, Vals, Card, 1.0, true);
				double b = P.MarginalCalculationWD(Data, Vars[0], Vals[0], Card, 1.0);
				double c = P.MarginalCalculationWD(Data, Vars[1], Vals[1], Card, 1.0);
				r = a / (b*c);
				double r1 = log2(r);
				t = r1 * a;
				Prob.push_back(t);
			}
			double pes = 0.0;
			for (int p = 0; p < Prob.size(); p++) {
				pes += Prob[p];
			}
			for (int h = 0; h < Vars.size(); h++) {
				cout << Vars[h] << " ";
			}
			peso.push_back(pes);
			cout << "\t" << pes << "\n";
			vector<double>temp;
			temp.push_back(pes);
			temp.push_back(Vars[0]);
			temp.push_back(Vars[1]);
			Complete.push_back(temp);
		}
		vector<vector<double>>MatrixPeso = T.MakeMatrix2(peso, N);
		//T.PrintMatrix2(MatrixPeso);

		sort(Complete.begin(), Complete.end(), T.sortcol);

		for (int f = 0; f < Complete.size(); f++) {
			int a = Complete[f][1];
			int b = Complete[f][2];

			Result[a][b] = 1;
			int count = 0;

			for (int m = 0; m < Result.size(); m++) {
				if (Result[m][b] == 1) {
					count++;
				}
			}

			if (count > 1) {
				Result[a][b] = 0;
			}
		}
		cout << "\n";	
		T.PrintMatrix(Result);

		return Result;
	}

	void VariableEliminaiton(vector<vector<int>> G, vector<vector<int>> Data, vector<int> Card, int val, vector<int>inference, double N) {

		vector<vector<int>> Dist;
		vector<int> NonObs;
		double Total = 1.0;

		for (int i = 0; i < N; i++) {
			vector<int> Padres = T.GetPadres(i, G);
			vector<int> Vars = T.Join(i, Padres);

			Dist.push_back(Vars);
		}

		for (int p = 0; p < N; p++) {
			if (find(inference.begin(), inference.end(), p) != inference.end()) {
				continue;
			}
			else {
				NonObs.push_back(p);
			}
		}

		int count = 0;

		for (int i = 0; i < NonObs.size(); i++) {
			for (int x = 0; x < Dist.size(); x++) {
				for (int y = 0; y < Dist[x].size(); y++) {
					if (Dist[x][y] == NonObs[i]) {
						Dist[x].erase(Dist[x].begin() + y);
					}
				}
			}
		}

		for (int i = 0; i < Dist.size(); i++) {
			if (Dist[i].size() == 0) {
				Dist.erase(Dist.begin() + i);
				i--;
			}
		}

		for (int n = 0; n < Dist.size(); n++) {
			vector<int> Vars = Dist[n];

			int FactorSize = 1;
			for (int v = 0; v < Vars.size(); v++) {
				FactorSize *= Card[Vars[v]];
			}

			vector<double> Probs(FactorSize);
			if (Vars.size() > 1) {
				Probs = Q.CalculateProbs(Data, Probs, Vars, Card, true);
			}
			else {
				Probs = Q.CalculateProbs(Data, Probs, Vars, Card, false);
			}
			vector<double> selecionados;

			double vf = 0.0;
			if (Vars.size() > 1 && Vars[1] == 1) {
				int cardPadre = Card[Vars[1]];
				for (int p = cardPadre * val; p < cardPadre*val + cardPadre; p++) {

					selecionados.push_back(Probs[p]);
				}
				vf = *max_element(std::begin(selecionados), std::end(selecionados));
			}

			else if (Vars.size() > 1 && Vars[1] != 1) {

				int cardHijo = Card[Vars[0]];
				int cardPadre = Card[Vars[1]];
				for (int p = val; p < Probs.size(); p = p + (cardHijo)) {
					selecionados.push_back(Probs[p]);
				}

				vf = *max_element(std::begin(selecionados), std::end(selecionados));
			}

			else if (Vars.size() == 1 && Vars[0] != 1) {
				vf = *max_element(std::begin(Probs), std::end(Probs));
			}
			else {
				 vf = Probs[val];
			}
			Total *= vf;
		}

		cout << Total;

	}
};
