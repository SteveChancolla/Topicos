#pragma once
#include <iostream>
#include <vector>
#include "Probabilities.h"
#include "Tools.h"

using namespace std;

class QMeasure
{
private:

public:
	Probabilities P;
	Tools T;

	double MarginalEntropy(vector<vector<int>> Data, int Var, int Val, vector<int> Card, double Alpha) {
		double Probability = 0.0;
		double Counter = 0.0;
		if (Data.size() > 0) {
			for (int i = 0; i < Data.size(); i++) {
				if (Data[i][Var] == Val) {
					Probability++;
				}
			}
			Counter = P.MarginalCountTable(Data, Var, Val);
			Probability = (Probability + Alpha) / (Data.size() + Card[Var] * Alpha);
			Probability = Counter * log10(Probability);
		}
		return Probability;
	}

	double JointEntropy(vector<vector<int>> Data, vector<int> Vars, vector<int> Vals, vector<int> Card, double Alpha) {
		double Probalility = 0.0;

		if (Data.size() > 0) {
			double Counter = 0.0;
			bool Boolean = false;
			double FullCard = 1.0;
			for (int i = 0; i < Data.size(); i++) {
				for (int j = 0; j < Vars.size(); j++) {
					if (Data[i][Vars[j]] != Vals[j]) {
						Boolean = false;
						break;
					}
					else {
						Boolean = true;
					}
				}
				if (Boolean) {
					Counter++;
				}
			}
			for (int i = 0; i < Vars.size(); i++) {
				FullCard = FullCard * Card[Vars[i]];
			}

			Probalility = (Counter + Alpha) / ((Data.size()) + (FullCard * Alpha));
			Probalility = Counter * log10(Probalility);
		}
		return Probalility;
	}

	vector<double> ConditionalEntropy(vector<vector<int>> Data, vector<int>Vars, vector<vector<int>> Vals, vector<int> Card) {
		double NumProb = 0.0;
		double DenProb = 0.0;
		int Counter = 0;
		vector<double> Probs(Card[Vars[0]]);
		vector<int> Aux;
		vector<int> Aux2;

		for (int k = 0; k < Vars.size(); k++) {
			Aux.push_back(Vars[k]);
		}
		for (int i = 0; i < Card[Vars[0]]; i++) {
			for (int j = 0; j < Vals[i].size(); j++) {
				Aux2.push_back(Vals[i][j]);
			}
			NumProb = P.JointCalculationWD(Data, Vars, Vals[i], Card, 1.0, false);

			if (Aux.size() > 1) {
				DenProb = P.JointCalculationWD(Data, Aux, Aux2, Card, 1.0, true);

			}
			else {
				DenProb = P.MarginalCalculationWD(Data, Aux[0], Aux2[0], Card, 1.0);
			}
			Probs[i] = NumProb / DenProb;

		}

		double Total = 0.0;
		double Suma = 0.0;
		for (int k = 0; k < Card[Vars[0]]; k++) {
			Total = Probs[k];
			Suma = Suma + Total;
		}
		for (int l = 0; l < Card[Vars[0]]; l++) {
			Probs[l] = Probs[l] / Suma;
		}
		for (int d = 0; d < Probs.size(); d++) {
			Counter = P.JointCountTable(Data, Vars, Vals[d]);
			Probs[d] = Counter * log10(Probs[d]);
		}
		return Probs;
	};

	vector<double> CalculateProbs(vector<vector<int>> Data, vector<double> Probs, vector<int> Vars, vector<int>Card, bool Conditional) {
		int N = Probs.size();
		int M = Vars.size();
		vector<int> Vals(M);
		int CardAnt = 1;
		int iVar;
		if (!Conditional) {
			for (int i = 0; i < N; i++) {
				CardAnt = 1;
				for (int j = 0; j < M; j++) {
					iVar = Vars[j];
					Vals[j] = (int)((int)floor(i / CardAnt) % Card[iVar]);
					CardAnt *= Card[iVar];
				}
				double Prob = 0.0;
				Prob = P.JointCalculationWD(Data, Vars, Vals, Card, 1.0, true);
				Probs[i] = Prob;
			}
		}
		else {
			vector<vector<int>>AuxVals(Card[Vars[0]], vector<int>(M));
			for (int i = 0; i < N; i += Card[Vars[0]]) {
				for (int j = 0; j < Card[Vars[0]]; j++) {
					CardAnt = 1;
					for (int k = 0; k < M; k++) {
						iVar = Vars[k];
						Vals[k] = (int)((int)floor((i + j) / CardAnt) % Card[iVar]);
						CardAnt *= Card[iVar];
					}
					AuxVals[j] = Vals;
				}
				vector<double> Prob;
				Prob = P.ConditionalCalculationWD(Data, Vars, AuxVals, Card, 1.0);

				for (int p = 0; p < Card[Vars[0]]; p++) {
					Probs[i + p] = Prob[p];
				}
			}
		}
		return Probs;
	}


	vector<double> CalculateEntropy(vector<vector<int>> Data, vector<double> Probs, vector<int> Vars, vector<int>Card, bool Conditional) {
		int N = Probs.size();
		int M = Vars.size();
		vector<int> Vals(M);
		int CardAnt = 1;
		int iVar;
		if (!Conditional) {
			for (int i = 0; i < N; i++) {
				CardAnt = 1;
				for (int j = 0; j < M; j++) {
					iVar = Vars[j];
					Vals[j] = (int)((int)floor(i / CardAnt) % Card[iVar]);
					CardAnt *= Card[iVar];
				}
				double Prob = 0.0;
				Prob = MarginalEntropy(Data, Vars[0], Vals[0], Card, 1.0);
				Probs[i] = Prob;
			}
		}
		else {
			vector<vector<int>>AuxVals(Card[Vars[0]], vector<int>(M));
			for (int i = 0; i < N; i += Card[Vars[0]]) {
				for (int j = 0; j < Card[Vars[0]]; j++) {
					CardAnt = 1;
					for (int k = 0; k < M; k++) {
						iVar = Vars[k];
						Vals[k] = (int)((int)floor((i + j) / CardAnt) % Card[iVar]);
						CardAnt *= Card[iVar];
					}
					AuxVals[j] = Vals;
				}
				vector<double> Prob;
				Prob = ConditionalEntropy(Data, Vars, AuxVals, Card);

				for (int p = 0; p < Card[Vars[0]]; p++) {
					Probs[i + p] = Prob[p];
				}
			}
		}
		return Probs;
	}

	vector<double> CalculateK(vector<vector<int>>Data, vector<double>Probs, vector<int> Vars, vector<int>Card, bool Conditional) {
		int N = Probs.size();
		int M = Vars.size();
		vector<int> Vals(M);
		int CardAnt = 1;
		int iVar;
		vector<vector<int>>AuxVals(Card[Vars[0]], vector<int>(M));
		for (int i = 0; i < N; i += Card[Vars[0]]) {
			for (int j = 0; j < Card[Vars[0]]; j++) {
				CardAnt = 1;
				for (int k = 0; k < M; k++) {
					iVar = Vars[k];
					Vals[k] = (int)((int)floor((i + j) / CardAnt) % Card[iVar]);
					CardAnt *= Card[iVar];
				}
				AuxVals[j] = Vals;
			}
			vector<double> Prob;
			int b;
			for (int m = 0; m < AuxVals.size(); m++) {
				int t = P.JointCountTable(Data, Vars, AuxVals[m]);
				double lg = log2(T.Factorial(t));
				Prob.push_back(lg);
			}
			int a = Card[Vars[0]];

		}
	}

	double GetEntropy(vector<vector<int>> G, vector<int> Card, vector<vector<int>> Data, int N) {
		double Entropy = 0.0;

		for (int i = 0; i < N; i++) {
			vector<int> Padres = T.GetPadres(i, G);
			vector<int> Vars = T.Join(i, Padres);
			int FactorSize = 1;
			for (int v = 0; v < Vars.size(); v++) {
				FactorSize *= Card[Vars[v]];
			}
			vector<double> Probs(FactorSize);
			Probs = CalculateEntropy(Data, Probs, Vars, Card, true);

			if (Vars.size() > 1) {
				for (int g = 0; g < Probs.size(); g++) {
					Entropy += Probs[g];
				}
			}
			else {
				for (int k = 0; k < Probs.size(); k++) {
					Entropy += Probs[k];
				}
			}

		}
		return Entropy;
	}

	int NIndependent(vector<vector<int>> Vars, vector<int> Card) {
		int K = 0;
		int CardP;

		for (int i = 0; i < Vars.size(); i++) {
			CardP = 1;
			for (int j = 1; j < Vars[i].size(); j++) {
				CardP *= Card[Vars[i][j]];
			}
			K += (Card[i] - 1)*CardP;
		}
		return K;
	};

	double GetAkaike(double Entropy, int K) {
		double Akaike = 0.0;
		Akaike = Entropy + K;

		return Akaike;
	}

	double GetMDL(double Entropy, double K, int Size) {
		double MDL = 0.0;
		MDL = Entropy + ((K / 2)*log2(Size));

		return MDL;
	}

	double GetK2(vector<vector<int>>G, vector<vector<int>>Data, vector<int>Card, int N, double ndga) {
		double K2 = 0.0;
		int CardAnt = 1;
		int iVar;
		vector<double> kv;
		for (int i = 0; i < N; i++) {
			double x = 0.0;
			vector<int> Padres = T.GetPadres(i, G);
			vector<int> Vars = T.Join(i, Padres);
			int FactorSize = 1;
			for (int v = 0; v < Vars.size(); v++) {
				FactorSize *= Card[Vars[v]];
			}
			vector<int> Vals(Vars.size());

			vector<double> Prob;

			vector<vector<int>>AuxVals(Card[Vars[0]], vector<int>(Vars.size()));
			for (int i = 0; i < FactorSize; i += Card[Vars[0]]) {
				for (int j = 0; j < Card[Vars[0]]; j++) {
					CardAnt = 1;
					for (int k = 0; k < Vars.size(); k++) {
						iVar = Vars[k];
						Vals[k] = (int)((int)floor((i + j) / CardAnt) % Card[iVar]);
						CardAnt *= Card[iVar];
					}
					AuxVals[j] = Vals;
				}
				for (int m = 0; m < AuxVals.size(); m++) {
					int t = P.JointCountTable(Data, Vars, AuxVals[m]);
					double lg = 0.0;
					for (int w = 1; w <= t; w++) {
						lg += log2(w);
					}
					Prob.push_back(lg);
				}
				if (Padres.size() >= 1) {
					vector<vector<int>> valsA;
					for (int f = 0; f < AuxVals.size(); f++) {
						vector<int>temp;
						for (int g = 0; g < AuxVals[f].size(); g++) {
							if (g != 0) {
								temp.push_back(AuxVals[f][g]);
							}
						}
						valsA.push_back(temp);
					}

					double CountP = P.JointCountTable(Data, Padres, valsA[0]);
					double CardVar = Card[Vars[0]];
					double Rest = CardVar - 1;
					double R = 0.0;
					for (int q = 1; q <= CountP; q++) {
						R += log2(1 / (Rest + q));
					}
					Prob.push_back(R);
				}
				else {
					int counP = 0;
					int a = Card[Vars[0]];
					double R = log2(((T.Factorial(a - 1)) / (T.Factorial(a - 1 + counP))));
					Prob.push_back(R);
				}
			}

			for (int t = 0; t < Prob.size(); t++) {
				x += Prob[t];
			}
			kv.push_back(x);
		}

		for (int y = 0; y < kv.size(); y++) {
			K2 += kv[y];
		}
		double pDga = 1 / ndga;
		return K2 + log2(pDga);
	}

};
