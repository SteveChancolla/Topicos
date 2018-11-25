#pragma once
#include<iostream>
#include<vector>

using namespace std;

class Probabilities
{
private:

public:
	int MarginalCountTable(vector<vector<int>> Data, int Var, int Val) {
		int Counter = 0;

		if (Data.size() > 0) {
			for (int i = 0; i < Data.size(); i++) {
				if (Data[i][Var] == Val) {
					Counter++;
				}
			}
		}
		return Counter;
	};

	int JointCountTable(vector<vector<int>> Data, vector<int> Vars, vector<int> Vals) {
		int Counter = 0;

		if (Data.size() > 0) {
			bool Boolean = false;
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
		}
		return Counter++;
	};

	double MarginalCalculation(vector<vector<int>> Data, int Var, int Val) {
		double Probability = 0.0;

		if (Data.size() > 0) {
			for (int i = 0; i < Data.size(); i++) {
				if (Data[i][Var] == Val) {
					Probability++;
				}
			}
			Probability /= Data.size();
		}
		return Probability;
	}

	double JointCalculation(vector<vector<int>> Data, vector<int> Vars, vector<int> Vals, vector<int> Card) {
		double Probability = 0.0;

		if (Data.size() > 0) {
			double Counter = 0;
			bool Boolean = false;
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
			Probability = Counter / Data.size();
		}
		return Probability;
	}

	vector<double> ConditionalCalculation(vector<vector<int>> Data, vector<int> Vars, vector<vector<int>> Vals, vector<int> Card) {
		vector<double> Probability(Card[Vars[0]]);
		vector<int> Aux;
		vector<int> Aux2;
		double Numerator = 0.0;
		double Denominator = 0.0;

		for (int k = 0; k < Vars.size() - 1; k++) {
			Aux.push_back(Vars[k + 1]);
		}
		for (int i = 0; i < Card[Vars[0]]; i++) {
			for (int j = 0; j < Vals[i].size() - 1; j++) {
				Aux2.push_back(Vals[i][j + 1]);
			}
			Numerator = JointCalculation(Data, Vars, Vals[0], Card);
			if (Numerator == 0) {
				Probability[i] = Numerator;
				break;
			}
			if (Aux.size() > 1) {
				Denominator = JointCalculation(Data, Vars, Vals[i], Card);
			}
			else {
				Denominator = MarginalCalculation(Data, Aux[0], Aux2[0]);
			}
			Probability[i] = Numerator / Denominator;
		}
		return Probability;
	}

	double MarginalCalculationWD(vector<vector<int>> Data, int Var, int Val, vector<int> Card, double Alpha) {
		double Probability = 0.0;

		if (Data.size() > 0) {
			for (int i = 0; i < Data.size(); i++) {
				if (Data[i][Var] == Val) {
					Probability++;
				}
			}
			Probability = (Probability + Alpha) / (Data.size() + (Card[Var] * Alpha));
		}
		return Probability;
	};

	double JointCalculationWD(vector<vector<int>> Data, vector<int> Vars, vector<int> Vals, vector<int> Card, double Alpha, bool Conditional) {
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
			if (Conditional && Counter == 0) {
				Counter = (Counter + Alpha) / (Data.size() + (FullCard * Alpha));
				return Counter;
			}
			Probalility = (Counter + Alpha) / ((Data.size()) + (FullCard * Alpha));
		}
		return Probalility;
	};

	vector<double> ConditionalCalculationWD(vector<vector<int>> Data, vector<int> Vars, vector<vector<int>> Vals, vector<int> Card, double alpha) {
		vector<double> Probabilities(Card[Vars[0]]);
		vector<int> Aux;
		vector<int> Aux2;
		double Numerator = 0.0;
		double Denominator = 0.0;

		for (int k = 0; k < Vars.size() - 1; k++) {
			Aux.push_back(Vars[k + 1]);
		}
		for (int i = 0; i < Card[Vars[0]]; i++) {
			for (int j = 0; j < Vals[i].size() - 1; j++) {
				Aux2.push_back(Vals[i][j + 1]);
			}
			Numerator = JointCalculationWD(Data, Vars, Vals[i], Card, alpha, false);
			if (Numerator == 0) {
				Probabilities[i] = Numerator;
				break;
			}
			if (Aux.size() > 1) {
				Denominator = JointCalculationWD(Data, Aux, Aux2, Card, alpha, true);
			}
			else {
				Denominator = MarginalCalculationWD(Data, Aux[0], Aux2[0], Card, alpha);
			}
			Probabilities[i] = Numerator / Denominator;
		}
		//Normalization
		double Sum = 0.0;
		double Total = 1.0;
		int AuxCounter = 0;
		for (int i = 0; i < Probabilities.size(); i++) {
			if (Probabilities[i] == 0) {
				AuxCounter++;
			}
			Sum += Probabilities[i];
		}
		if (Sum != 1.0) {
			Total -= Sum;
			Total /= (AuxCounter);
			for (int i = 0; i < Probabilities.size(); i++) {
				if (Probabilities[i] == 0) {
					Probabilities[i] = Total;
				}
			}
		}
		Sum = 0;
		for (int i = 0; i < Probabilities.size(); i++) {
			Sum += Probabilities[i];
		}
		return Probabilities;
	}

	void ShowPmComandWD(int var, vector<vector<int>>Data, vector<string>VarNames, vector<vector<string>>VarVals, vector<int>Card,int Alpha) {
		int w;
		int CardAnt = 1;
		int iVar;
		cout << "\n";
		cout << "Factor " << VarNames[var];
		cout << "\n";
		cout << "i" << "\t";
		vector<int>Vars;
		Vars.push_back(var);
		int m = Vars.size();
		vector<int>vals(m);
		for (w = 0; w < Vars.size() - 1; w++) {
			cout << VarNames[Vars[w]] << "\t";
		}
		cout << VarNames[Vars[w]] << endl;

		for (int i = 0; i < VarVals[var].size(); i++) {
			CardAnt = 1;
			for (int j = 0; j < m; j++) {
				iVar = Vars[j];
				vals[j] = (int)floor(i / CardAnt) % Card[iVar];
				CardAnt *= Card[iVar];
			}
			double prob = 0.0;

			prob = MarginalCalculationWD(Data, Vars[0], vals[0], Card, Alpha);

			cout << i << "\t";
			for (w = 0; w < vals.size() - 1; w++) {
				cout << VarVals[Vars[w]][vals[w]] + "\t";
			}
			cout << VarVals[Vars[w]][vals[w]] + "\t" << prob;
			cout << "\n";
		}
	}
	
	void ShowPmComand(int var, vector<vector<int>>Data, vector<string>VarNames, vector<vector<string>>VarVals, vector<int>Card) {
		int w;
		int CardAnt = 1;
		int iVar;
		cout << "\n";
		cout << "Factor " << VarNames[var];
		cout << "\n";
		cout << "i" << "\t";
		vector<int>Vars;
		Vars.push_back(var);
		int m = Vars.size();
		vector<int>vals(m);
		for (w = 0; w < Vars.size() - 1; w++) {
			cout << VarNames[Vars[w]] << "\t";
		}
		cout << VarNames[Vars[w]] << endl;

		for (int i = 0; i < VarVals[var].size(); i++) {
			CardAnt = 1;
			for (int j = 0; j < m; j++) {
				iVar = Vars[j];
				vals[j] = (int)floor(i / CardAnt) % Card[iVar];
				CardAnt *= Card[iVar];
			}
			double prob = 0.0;

			prob = MarginalCalculation(Data, Vars[0], vals[0]);

			cout << i << "\t";
			for (w = 0; w < vals.size() - 1; w++) {
				cout << VarVals[Vars[w]][vals[w]] + "\t";
			}
			cout << VarVals[Vars[w]][vals[w]] + "\t" << prob;
			cout << "\n";
		}
	}

	void ShowPconjComandWD(vector<int> Vars, vector<vector<int>>Data, vector<string>VarNames, vector<vector<string>>VarVals, vector<int>Card, int Alpha) {
		int w;
		int CardAnt = 1;
		int iVar;
		cout << "\n";
		cout << "Factor ";
		cout << "\n";
		cout << "i" << "\t";
		int m = Vars.size();
		vector<int>vals(m);
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

			Prob = JointCalculationWD(Data, Vars, vals, Card, Alpha, false);

			cout << i << "\t";
			for (w = 0; w < vals.size() - 1; w++) {
				cout << VarVals[Vars[w]][vals[w]] + "\t";
			}
			cout << VarVals[Vars[w]][vals[w]] + "\t" << Prob;
			cout << "\n";
		}
	}
	
	void ShowPconjComand(vector<int> Vars, vector<vector<int>>Data, vector<string>VarNames, vector<vector<string>>VarVals, vector<int>Card) {
		int w;
		int CardAnt = 1;
		int iVar;
		cout << "\n";
		cout << "Factor ";
		cout << "\n";
		cout << "i" << "\t";
		int m = Vars.size();
		vector<int>vals(m);
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

			Prob = JointCalculation(Data, Vars, vals, Card);

			cout << i << "\t";
			for (w = 0; w < vals.size() - 1; w++) {
				cout << VarVals[Vars[w]][vals[w]] + "\t";
			}
			cout << VarVals[Vars[w]][vals[w]] + "\t" << Prob;
			cout << "\n";
		}
	}

	vector<double> ShowPcondComandWD(vector<int> Vars, vector<vector<int>>Data, vector<double>Probs, vector<string>VarNames, vector<vector<string>>VarVals, vector<int>Card,int Alpha) {
		vector<double>Probabilities;
		int N = Probs.size();
		int w;
		int CardAnt = 1;
		int iVar;
		cout << "\n";
		cout << "Factor ";
		cout << "\n";
		cout << "i" << "\t";
		int m = Vars.size();
		vector<int>vals(m);
		for (w = 0; w < Vars.size() - 1; w++) {
			cout << VarNames[Vars[w]] << "\t";
		}
		cout << VarNames[Vars[w]] << endl;
		int FactorSize = 1;
		for (int k = 0; k < Vars.size(); k++) {
			FactorSize *= Card[Vars[k]];
		}
		vector<vector<int>>AuxVals(Card[Vars[0]], vector<int>(m));
		for (int i = 0; i < N; i += Card[Vars[0]]) {
			for (int k = 0; k < Card[Vars[0]]; k++) {
				CardAnt = 1;
				for (int j = 0; j < m; j++) {
					iVar = Vars[j];
					vals[j] = (int)((int)floor((i + k) / CardAnt) % Card[iVar]);
					CardAnt *= Card[iVar];
				}
				AuxVals[k] = vals;
			}
			vector<double> Prob;
			Prob = ConditionalCalculationWD(Data, Vars, AuxVals, Card,Alpha);
			//double probs;

			for (int p = 0; p < Card[Vars[0]]; p++) {

				double Total = 0.0;
				double Sum = 0.0;
				for (int h = 0; h < Card[Vars[0]]; h++) {
					Total = Prob[h];
					Sum = Sum + Total;
				}
				Probs[p + i] = Prob[p] / Sum;
				cout << (i + p) << "\t";
				for (w = 0; w < vals.size() - 1; w++) {
					cout << VarVals[Vars[w]][AuxVals[p][w]] + "\t";
				}
				cout << VarVals[Vars[w]][AuxVals[p][w]] + "\t" << Prob[p] / Sum;
				cout << "\n";
			}
		}
		return Probs;
	}

	vector<double> ShowPcondComand(vector<int> Vars, vector<vector<int>>Data, vector<double>Probs, vector<string>VarNames, vector<vector<string>>VarVals, vector<int>Card) {
		vector<double>Probabilities;
		int N = Probs.size();
		int w;
		int CardAnt = 1;
		int iVar;
		cout << "\n";
		cout << "Factor ";
		cout << "\n";
		cout << "i" << "\t";
		int m = Vars.size();
		vector<int>vals(m);
		for (w = 0; w < Vars.size() - 1; w++) {
			cout << VarNames[Vars[w]] << "\t";
		}
		cout << VarNames[Vars[w]] << endl;
		int FactorSize = 1;
		for (int k = 0; k < Vars.size(); k++) {
			FactorSize *= Card[Vars[k]];
		}
		vector<vector<int>>AuxVals(Card[Vars[0]], vector<int>(m));
		for (int i = 0; i < N; i += Card[Vars[0]]) {
			for (int k = 0; k < Card[Vars[0]]; k++) {
				CardAnt = 1;
				for (int j = 0; j < m; j++) {
					iVar = Vars[j];
					vals[j] = (int)((int)floor((i + k) / CardAnt) % Card[iVar]);
					CardAnt *= Card[iVar];
				}
				AuxVals[k] = vals;
			}
			vector<double> Prob;
			Prob = ConditionalCalculation(Data, Vars, AuxVals, Card);
			//double probs;

			for (int p = 0; p < Card[Vars[0]]; p++) {

				double Total = 0.0;
				double Sum = 0.0;
				for (int h = 0; h < Card[Vars[0]]; h++) {
					Total = Prob[h];
					Sum = Sum + Total;
				}
				Probs[p + i] = Prob[p] / Sum;
				cout << (i + p) << "\t";
				for (w = 0; w < vals.size() - 1; w++) {
					cout << VarVals[Vars[w]][AuxVals[p][w]] + "\t";
				}
				cout << VarVals[Vars[w]][AuxVals[p][w]] + "\t" << Prob[p] / Sum;
				cout << "\n";
			}
		}
		return Probs;
	}

};
