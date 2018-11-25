#include <iostream>
#include <vector>
#include <math.h>
#include <sstream>
#include <fstream>
#include <stdint.h>
#include <fstream>
#include "DAG.h"
#include "QMeasure.h"
#include "Tools.h"
#include "Algorithms.h"
#include <algorithm>

using namespace std;

int main() {

	double N;
	string input;
	int Alpha = 0;
	int NroDagsFB = 0;
	int NroDagsK2 = 0;

	DAG *dc;
	Tools Tools;
	Algorithms Alg;
	Probabilities P;
	QMeasure QMeasures;

	vector<int> Card;
	vector<string> VarNames;
	vector<vector<string>> Data1;
	vector<vector<int>> Data;
	vector<vector<string>> VarVals;
	vector<vector<int>> ChowGraph;

	do {
		vector<string> Temp = Tools.GetComand(input);
		if (Temp[0] == "data") {
			Data1 = Tools.ReadingData(Temp[1]);
			Tools.PrintMatrix1(Data1);
			cout << "\tNumero de instancias: " << Data1.size() << endl;
			cout << "\tNumero de variables: " << Data1[0].size() << endl;
			N = Data1[0].size();
			cout << "\n";
		}
		else if (Temp[0] == "vars") {
			for (int i = 0; i < Temp.size() - 1; i++) {
				VarNames.push_back(Temp[i + 1]);
			}
			for (int j = 0; j < VarNames.size(); j++) {
				cout << "[" << j + 1 << "] : " << VarNames[j] << endl;
			}
			cout << "\n";
		}
		else if (Temp[0] == "vals") {
			vector<string>aux;
			for (int i = 0; i < Temp.size() - 2; i++) {
				aux.push_back(Temp[i + 2]);
			}
			VarVals.push_back(aux);
			Card.push_back(aux.size());

			cout << "\n";
		}
		else if (Temp[0] == "card") {
			for (int i = 0; i < Card.size(); i++) {
				cout << "\t" << Card[i] << endl;
			}
			cout << "\n";
		}
		else if (Temp[0] == "vals.show") {
			for (int i = 0; i < VarVals.size(); i++) {
				cout << "[" << VarNames[i] << "]";
				for (int j = 0; j < VarVals[i].size(); j++) {
					cout << VarVals[i][j] << " ";
				}
				cout << "\n";
			}
			cout << "\n";
			Data = Tools.Dataint(Data1, VarVals);
			Tools.PrintMatrix(Data);
			cout << "\n";
		}
		else if (Temp[0] == "alpha") {
			Alpha = stoi(Temp[1]);
			cout << "Alpha = " << Alpha << endl;
			cout << "\n";
		}
		else if (Temp[0] == "pm") {
			int pos1 = find(VarNames.begin(), VarNames.end(), Temp[1]) - VarNames.begin();
			if (Alpha > 0) {
				P.ShowPmComandWD(pos1, Data, VarNames, VarVals, Card, Alpha);
			}
			else {
				P.ShowPmComand(pos1, Data, VarNames, VarVals, Card);
			}
			cout << "\n";
		}
		else if (Temp[0] == "pconj") {
			vector<int>vars;

			for (int i = 1; i < Temp.size(); i++) {
				int pos1 = find(VarNames.begin(), VarNames.end(), Temp[i]) - VarNames.begin();
				vars.push_back(pos1);
			}
			if (Alpha > 0) {
				P.ShowPconjComandWD(vars, Data, VarNames, VarVals, Card, Alpha);
			}
			else {
				P.ShowPconjComand(vars, Data, VarNames, VarVals, Card);
			}
			cout << "\n";
		}
		else if (Temp[0] == "pcond") {
			vector<int>vars;
			for (int i = 1; i < Temp.size(); i++) {
				if (i != 2) {
					int pos1 = find(VarNames.begin(), VarNames.end(), Temp[i]) - VarNames.begin();
					vars.push_back(pos1);
				}
			}
			int TamFactor = 1;
			for (int j = 0; j < vars.size(); j++) {
				TamFactor *= Card[vars[j]];
			}
			vector<double> probs(TamFactor);
			if (Alpha > 0) {
				probs = P.ShowPcondComandWD(vars, Data, probs, VarNames, VarVals, Card, Alpha);
			}
			else {
				probs = P.ShowPcondComand(vars, Data, probs, VarNames, VarVals, Card);
			}
			cout << "\n";
		}
		else if (Temp[0] == "chow") {
			vector<vector<int>> Matrix = { { 1,1,1,1 },{ 1,1,1,1 },{ 1,1,1,1 },{ 1,1,1,1 } };

			ChowGraph = Alg.ChowLiu(Matrix, Card, Data, N);

			cout << "\n";
		}
		else if (Temp[0] == "ve") {

			vector<int> inferencia;
			for (int i = 1; i < Temp.size() - 1; i++) {
				int pos1 = find(VarNames.begin(), VarNames.end(), Temp[i]) - VarNames.begin();
				inferencia.push_back(pos1);
			}
			int valor = stoi(Temp[3]);


			Alg.VariableEliminaiton(ChowGraph, Data, Card, valor, inferencia, N);

		}
		else if (Temp[0] == "mi") {
			vector<int>vars;
			for (int i = 1; i < Temp.size(); i++) {
				if (i != 3) {
					int pos1 = find(VarNames.begin(), VarNames.end(), Temp[i]) - VarNames.begin();
					vars.push_back(pos1);
				}
			}
			Alg.MutualInformation(Data, vars, Card);
			cout << "\n";
		}
		else if (Temp[0] == "pearson") {
			vector<int>vars;
			for (int i = 1; i < Temp.size(); i++) {
				if (i != 3) {
					int pos1 = find(VarNames.begin(), VarNames.end(), Temp[i]) - VarNames.begin();
					vars.push_back(pos1);
				}
			}
			Alg.Pearson(Data, vars, Card);
		}
		else if (Temp[0] == "marginalizacion")
		{
			vector<int>vars;


			vector<string>::iterator itr = find(Temp.begin(), Temp.end(), "m");

			int posm = distance(Temp.begin(), itr);

			vector<int>marg;
			for (int i = posm + 1; i < Temp.size(); i++) {
				int pos1 = find(VarNames.begin(), VarNames.end(), Temp[i]) - VarNames.begin();
				marg.push_back(pos1);
			}
			int size = marg.size();


			for (int i = 1; i < posm; i++) {
				int pos1 = find(VarNames.begin(), VarNames.end(), Temp[i]) - VarNames.begin();
				vars.push_back(pos1);
			}

			Alg.Marginalization(Data, vars, Card, VarNames, VarVals, size, marg);
		}
		else if (Temp[0] == "dagsFB") {
			int n = (int)pow(2, (pow(N, 2) - N));
			int nAristas = pow(N, 2) - N;
			int count = 0;
			int count2 = 0;
			int k = 0;
			vector<int> Aux;
			vector<vector<int>> Result;
			for (int i = 0; i < n; i++) {
				Aux = Alg.BruteForce(i, nAristas);

				count++;
				Result = Tools.MakeMatrix(Aux, N);

				dc = new DAG(Result);

				if (dc->HasCycle()) {
					continue;
				}
				else {
					count2++;
				}
			}
			NroDagsFB = count2;
			cout << "Numero de Grafos por Fuerza Bruta : " << count << endl;
			cout << "Numero de Dags por Fuerza Bruta : " << count2 << endl;
			cout << "\n";
		}
		else if (Temp[0] == "entropiaFB") {
			int n = (int)pow(2, (pow(N, 2) - N));
			int nAristas = pow(N, 2) - N;
			int count = 0;
			int count2 = 0;
			int k = 0;
			ofstream file;
			vector<int> Aux;
			vector<vector<int>> Result;
			file.open("EntropiaFB.txt");

			for (int i = 0; i < n; i++) {
				Aux = Alg.BruteForce(i, nAristas);

				count++;
				Result = Tools.MakeMatrix(Aux, N);

				dc = new DAG(Result);

				if (dc->HasCycle()) {
					continue;
				}
				else {
					cout << "\n";
					//Tools.PrintMatrix(Result);
					double Entropy = QMeasures.GetEntropy(Result, Card, Data, N);
					cout << "Score Entropia del grafo " << i << " : " << Entropy;
					for (int i = 0; i < Result.size(); i++) {
						for (int j = 0; j < Result.size(); j++) {
							file << Result[i][j];
						}
						file << "\n";
					}
					file << "\n";
					file << "Entropia: " << Entropy;
				}
			}
			file.close();

			cout << "\n";
		}
		else if (Temp[0] == "akaikeFB") {
			int n = (int)pow(2, (pow(N, 2) - N));
			int nAristas = pow(N, 2) - N;
			vector<int> Aux;
			vector<vector<int>> Result;
			ofstream file2;
			file2.open("AkaikeFB.txt");

			for (int i = 0; i < n; i++) {
				Aux = Alg.BruteForce(i, nAristas);
				Result = Tools.MakeMatrix(Aux, N);

				dc = new DAG(Result);

				if (dc->HasCycle()) {
					continue;
				}
				else {
					cout << "\n";
					//Tools.PrintMatrix(Result);
					vector<vector<int>> Vars = Tools.GetVars(N, Result);
					int k = QMeasures.NIndependent(Vars, Card);
					double Entropy = QMeasures.GetEntropy(Result, Card, Data, N);
					double Akaike = QMeasures.GetAkaike(Entropy, k);
					cout << "Score Akaike del grafo " << i << " : " << Akaike;
					for (int i = 0; i < Result.size(); i++) {
						for (int j = 0; j < Result.size(); j++) {
							file2 << Result[i][j];
						}
						file2 << "\n";
					}
					file2 << "\n";
					file2 << "Akaike: " << Akaike << endl;
				}
			}
			cout << "\n";
		}
		else if (Temp[0] == "mdlFB") {
			int n = (int)pow(2, (pow(N, 2) - N));
			int nAristas = pow(N, 2) - N;
			int k = 0;
			vector<int> Aux;
			vector<vector<int>> Result;
			ofstream file3;
			file3.open("DMLFB.txt");

			for (int i = 0; i < n; i++) {
				Aux = Alg.BruteForce(i, nAristas);
				Result = Tools.MakeMatrix(Aux, N);

				dc = new DAG(Result);

				if (dc->HasCycle()) {
					continue;
				}
				else {
					cout << "\n";
					//Tools.PrintMatrix(Result);
					vector<vector<int>> Vars = Tools.GetVars(N, Result);
					double Entropy = QMeasures.GetEntropy(Result, Card, Data, N);
					int k = QMeasures.NIndependent(Vars, Card);
					double MDL = QMeasures.GetMDL(Entropy, k, Data.size());
					cout << "Score MDL del Grafo " << i << " : " << MDL;
					for (int i = 0; i < Result.size(); i++) {
						for (int j = 0; j < Result.size(); j++) {
							file3 << Result[i][j];
						}
						file3 << "\n";
					}
					file3 << "\n";
					file3 << "MDL: " << MDL << endl;
				}
			}
			cout << "\n";
		}
		else if (Temp[0] == "k2FB") {
			int n = (int)pow(2, (pow(N, 2) - N));
			int nAristas = pow(N, 2) - N;
			int count = 0;
			int count2 = 0;
			int k = 0;
			vector<int> Aux;
			vector<vector<int>> Result;
			ofstream file4;
			file4.open("K2FB.txt");

			for (int i = 0; i < n; i++) {
				Aux = Alg.BruteForce(i, nAristas);

				count++;
				Result = Tools.MakeMatrix(Aux, N);

				dc = new DAG(Result);

				if (dc->HasCycle()) {
					continue;
				}
				else {
					cout << "\n";
					//Tools.PrintMatrix(Result);
					double K2 = QMeasures.GetK2(Result, Data, Card, N, NroDagsFB);
					cout << "Score K2 del grafo " << i << " : " << K2;
					for (int i = 0; i < Result.size(); i++) {
						for (int j = 0; j < Result.size(); j++) {
							file4 << Result[i][j];
						}
						file4 << "\n";
					}
					file4 << "\n";
					file4 << "K2: " << K2 << endl;
				}
			}
			std::cout << "\n";
		}
		else if (Temp[0] == "entropiaK2") {
			ofstream file5;
			file5.open("entropiaK2.txt");
			int padres = stoi(Temp[1]);

			vector<vector<int>> Result(N, vector<int>(N));
			vector<vector<int>> graph = Result;
			vector<vector<int>> graph2 = Result;
			int Contador = 0;
			int Contador2 = 0;
			int esOnoes = 0;

			double entropiaini = QMeasures.GetEntropy(Result, Card, Data, N);
			std::cout << "Score Entropia Inicial " << entropiaini << endl;

			for (int i = 0; i < graph.size(); i++) {
				for (int j = 0; j < graph.size(); j++) {
					file5 << Result[i][j];
				}
				file5 << "\n";
			}
			file5 << entropiaini << endl;

			do {
				esOnoes = 0;
				for (int i = 0; i < graph.size(); i++) {
					Contador2 = 0;
					for (int j = 0; j < graph.size(); j++) {
						if (graph[i][j] == 0 && !(j >= i)) {
							Contador2++;

							if (Contador2 > padres) {
								continue;
							}

							graph[i][j] = 1;

							dc = new DAG(graph);

							if (dc->HasCycle()) {
								continue;
							}
							else {

								for (int i = 0; i < graph.size(); i++) {
									for (int j = 0; j < graph.size(); j++) {
										file5 << graph[i][j];
									}
									file5 << "\n";
								}
								Contador++;
								double e = 0.0;
								double entropia2 = QMeasures.GetEntropy(graph, Card, Data, N);
								e = entropia2;
								std::cout << "Score Entropia del grafo " << Contador << " : " << entropia2 << endl;

								file5 << entropia2 << endl;
								if (e < entropiaini) {
									entropiaini = e;
									graph2 = graph;
									esOnoes++;
									graph[i][j] = 0;
								}
								else {
									graph[i][j] = 0;
								}
							}
						}
					}
				}
				esOnoes = 0;
				graph = graph2;
			} while (esOnoes > 0);

			std::cout << "\n";
		}
		else if (Temp[0] == "akaikeK2") {
			ofstream file6;
			file6.open("Akaikek2.txt");

			int padres = stoi(Temp[1]);

			vector<vector<int>> Result(N, vector<int>(N));
			vector<vector<int>> graph = Result;
			vector<vector<int>> graph2 = Result;
			int Contador = 0;
			int Contador2 = 0;
			int esOnoes = 0;


			double entropia = QMeasures.GetEntropy(Result, Card, Data, N);
			vector<vector<int>> Vars = Tools.GetVars(N, Result);
			int k = QMeasures.NIndependent(Vars, Card);
			double akaikeini = QMeasures.GetAkaike(entropia, k);
			cout << "Score Akaike incial: " << akaikeini << endl;

			for (int i = 0; i < graph.size(); i++) {
				for (int j = 0; j < graph.size(); j++) {
					file6 << Result[i][j];
				}
				file6 << "\n";
			}
			file6 << akaikeini << endl;

			do {
				esOnoes = 0;
				for (int i = 0; i < graph.size(); i++) {
					Contador2 = 0;

					for (int j = 0; j < graph.size(); j++) {
						if (graph[i][j] == 0 && !(j >= i)) {

							Contador2++;

							if (Contador2 > padres) {
								continue;
							}

							graph[i][j] = 1;

							dc = new DAG(graph);

							if (dc->HasCycle()) {
								continue;
							}
							else {
								Contador++;

								double entropia2 = QMeasures.GetEntropy(graph, Card, Data, N);
								vector<vector<int>> vars2 = Tools.GetVars(N, graph);
								int ki2 = QMeasures.NIndependent(vars2, Card);
								double akaike = QMeasures.GetAkaike(entropia2, ki2);

								for (int i = 0; i < graph.size(); i++) {
									for (int j = 0; j < graph.size(); j++) {
										file6 << graph[i][j];
									}
									file6 << "\n";
								}
								std::cout << "Score Akaike del grafo " << Contador << " : " << akaike << endl;
								file6 << entropia2 << endl;

								if (akaike < akaikeini) {
									akaikeini = akaike;
									graph2 = graph;
									esOnoes++;
									graph[i][j] = 0;
								}
								else {
									graph[i][j] = 0;
								}
							}
						}
					}
				}
				esOnoes = 0;
				graph = graph2;
			} while (esOnoes > 0);

			std::cout << "\n";
		}
		else if (Temp[0] == "mdlK2") {
			ofstream file7;

			file7.open("MDLK2.txt");
			int padres = stoi(Temp[1]);

			vector<vector<int>> Result(N, vector<int>(N));
			vector<vector<int>> graph = Result;
			vector<vector<int>> graph2 = Result;
			int Contador = 0;
			int Contador2 = 0;
			int esOnoes = 0;



			double entropia = QMeasures.GetEntropy(Result, Card, Data, N);
			vector<vector<int>> Vars = Tools.GetVars(N, Result);
			int k = QMeasures.NIndependent(Vars, Card);
			double MDLinicial = QMeasures.GetMDL(entropia, k, Data.size());
			cout << "Score MDL inicial: " << MDLinicial << endl;

			for (int i = 0; i < graph.size(); i++) {
				for (int j = 0; j < graph.size(); j++) {
					file7 << Result[i][j];
				}
				file7 << "\n";
			}

			file7 << MDLinicial << endl;
			do {
				esOnoes = 0;
				for (int i =0; i < graph.size(); i++) {
					Contador2 = 0;
					for (int j = 0; j < graph.size(); j++) {
						if (graph[i][j] == 0 && !(j >= i)) {
							Contador2++;

							if (Contador2 > padres) {
								continue;
							}

							graph[i][j] = 1;

							dc = new DAG(graph);

							if (dc->HasCycle()) {
								continue;
							}
							else {
								Contador++;

								double entropia2 = QMeasures.GetEntropy(graph, Card, Data, N);
								vector<vector<int>> vars2 = Tools.GetVars(N, graph);
								int ki2 = QMeasures.NIndependent(vars2, Card);
								double mdl2 = QMeasures.GetMDL(entropia2, ki2, Data.size());
								std::cout << "Score MDL del grafo " << Contador << " : " << mdl2 << endl;

								for (int i = 0; i < graph.size(); i++) {
									for (int j = 0; j < graph.size(); j++) {
										file7 << graph[i][j];
									}
									file7 << "\n";
								}
								file7 << mdl2 << endl;


								if (mdl2 < MDLinicial) {
									MDLinicial = mdl2;
									graph2 = graph;
									esOnoes++;
									graph[i][j] = 0;
								}
								else {
									graph[i][j] = 0;
								}
							}
						}
					}
				}
				esOnoes = 0;
				graph = graph2;
			} while (esOnoes > 0);

			std::cout << "\n";
		}
		else if (Temp[0] == "k2K2") {
			ofstream file8;
			file8.open("k2K2.txt");

			int padres = stoi(Temp[1]);

			vector<vector<int>> Result(N, vector<int>(N));
			vector<vector<int>> graph = Result;
			vector<vector<int>> graph2 = Result;
			int Contador = 0;
			int Contador2 = 0;
			int esOnoes = 0;

			double eInicial = 0.0;
			double eCopia = 0.0;
			double k2Inicial = 0.0;
			k2Inicial = QMeasures.GetK2(Result, Data, Card, N, 64);
			std::cout <<"Score K2 inicial: " <<k2Inicial << endl;

			for (int i = 0; i < graph.size(); i++) {
				for (int j = 0; j < graph.size(); j++) {
					file8 << Result[i][j];
				}
				file8 << "\n";
			}
			file8 << k2Inicial << endl;

			do {
				esOnoes = 0;
				for (int i = 0; i < graph.size(); i++) {
					Contador2 = 0;

					for (int j = 0; j < graph.size(); j++) {
						if (graph[i][j] == 0 && !(j >= i)) {
							Contador2++;

							if (Contador2 > padres) {
								continue;
							}

							graph[i][j] = 1;

							dc = new DAG(graph);

							if (dc->HasCycle()) {
								continue;
							}
							else {
								Contador++;

								double k2 = 0.0;
								k2 = QMeasures.GetK2(graph, Data, Card, N, 64);
								std::cout <<"Score K2 del grafo "<< Contador <<" : "<<k2<<endl;

								for (int i = 0; i < graph.size(); i++) {
									for (int j = 0; j < graph.size(); j++) {
										file8 << graph[i][j];
									}
									file8 << "\n";
								}
								file8 << k2 << endl;

								if (k2 > k2Inicial) {
									k2Inicial = k2;
									graph2 = graph;
									esOnoes++;
									graph[i][j] = 0;
								}
								else {
									graph[i][j] = 0;
								}
							}
						}
					}
				}
				esOnoes = 0;
				graph = graph2;
			} while (esOnoes > 0);
			cout << "\n";
		}
		else if (Temp[0] == "demo")
		{

			vector<vector<string>> Datax = Tools.ReadingData("Demo.csv");
			int N = Datax[0].size();

			vector<string> VarNames = { "Conocimiento","Estudia","Asite","Aprueba" };
			vector<vector<string>> VarVals = { {"bajo","medio","alto"},
			{"bajo","medio","alto"},
			{"no","si"},
			{"no","si"}
			};

			vector<int>Card = { 3,3,2,2 };

			vector<vector<int>> DataR = Tools.Dataint(Datax, VarVals);
			int NroDagsFB = 0;
			int NroDagsK2 = 0;

			int n = (int)pow(2, (pow(N, 2) - N));
			int nAristas = pow(N, 2) - N;
			int count = 0;
			int count2 = 0;
			int k = 0;
			vector<int> Aux;
			vector<vector<int>> Result2;
			for (int i = 0; i < n; i++) {
				Aux = Alg.BruteForce(i, nAristas);

				count++;
				Result2 = Tools.MakeMatrix(Aux, N);

				dc = new DAG(Result2);

				if (dc->HasCycle()) {
					continue;
				}
				else {
					count2++;
				}
			}
			NroDagsFB = count2;
			cout << "Numero de Grafos por Fuerza Bruta : " << count << endl;
			cout << "Numero de Dags por Fuerza Bruta : " << count2 << endl;
			cout << "\n";

			int a = pow(2, N);
			NroDagsK2 = a;
			cout << "Numero de Grafos por K2 : " << a << endl;
			cout << "Numero de Dags por K2 : " << a << endl;
			cout << "\n";

			vector<vector<int>> Result1;

			for (int i = 0; i < n; i++) {
				Aux = Alg.BruteForce(i, nAristas);
				Result1 = Tools.MakeMatrix(Aux, N);

				dc = new DAG(Result1);

				if (dc->HasCycle()) {
					continue;
				}
				else {
					cout << "\n";
					Tools.PrintMatrix(Result1);
					vector<vector<int>> Vars = Tools.GetVars(N, Result1);
					int k = QMeasures.NIndependent(Vars, Card);
					double Entropy = QMeasures.GetEntropy(Result1, Card, DataR, N);
					double Akaike = QMeasures.GetAkaike(Entropy, k);
					cout << "Akaike FB : " << Akaike << endl;

				}
			}
			cout << "\n";
			cout << "////////////////////////////////////";

			vector<vector<int>> Result(N, vector<int>(N));
			int Contador = 0;
			int esOnoes = 0;
			bool esK2 = false;
			vector<vector<int>> graph = Result;
			vector<vector<int>> graph2 = Result;
			vector<double> k2Entropia;
			double eInicial = 0.0;
			double eCopia = 0.0;
			double k2Inicial = 0.0;
			k2Inicial = QMeasures.GetK2(Result, DataR, Card, N, NroDagsK2);
			Tools.PrintMatrix(Result);
			std::cout << k2Inicial << endl;

			do {
				esOnoes = 0;
				for (int i = 0; i < graph.size(); i++) {
					for (int j = 0; j < graph.size(); j++) {
						if (graph[i][j] == 0 && !(j >= i)) {
							graph[i][j] = 1;

							dc = new DAG(graph);

							if (dc->HasCycle()) {
								continue;
							}
							else {
								Contador++;
								double k2 = 0.0;

								Tools.PrintMatrix(graph);

								k2 = QMeasures.GetK2(graph, DataR, Card, N, NroDagsK2);
								std::cout << "K2 K2" << k2;

								if (k2 > k2Inicial) {
									k2Inicial = k2;
									graph2 = graph;
									esOnoes++;
									graph[i][j] = 0;
								}
								else {
									graph[i][j] = 0;
								}
							}
						}
					}
				}
				graph = graph2;
			} while (esOnoes > 0);

			vector<vector<int>> Matrix = { { 1,1,1,1 },{ 1,1,1,1 },{ 1,1,1,1 },{ 1,1,1,1 } };

			ChowGraph = Alg.ChowLiu(Matrix, Card, Data, N);

			Tools.PrintMatrix(ChowGraph);
			cout << "\n";

			double t = QMeasures.GetEntropy(ChowGraph, Card, Data, 4);
			vector<vector<int>> Vars = Tools.GetVars(N, ChowGraph);

			int m = QMeasures.NIndependent(Vars, Card);
			double Entropy = QMeasures.GetMDL(t, m, 4);
			cout << "Calcular MDL con Chow Liu" << Entropy << endl;

			double k2 = QMeasures.GetK2(ChowGraph, Data, Card, 4, 1);
			cout << "Calculo de K2 con Chow liu" << k2 << endl;

			double aic = QMeasures.GetAkaike(a, k);
			cout << "Calculo de aic con chow liu" << aic << endl;


		}
	} while (input != "salir");


	getchar();
	getchar();
	return 0;
}