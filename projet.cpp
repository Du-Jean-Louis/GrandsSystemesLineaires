#include "fredholm.hpp"
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include "vectorops.hpp"

using namespace std;

int main(void){
	
	// Question 1
	
	cout << endl << "Test des fonctions de la classe FredholmMatrix (Question 1)" << endl << endl;
	
	FredholmMatrix a(3);
	
	a.insert(2,2,-1);
	a.insert(2,2,3);
	
	FredholmMatrix abis = a;
	cout << "Matrice abis issue de l'operateur =" << endl << abis << endl;
	
	vector<double> c(3,0.1);
	
	a.insert(c, c);
	a.insert(c, c);
	
	cout << "Matrice a" << endl << a << endl;
	
	vector<double> d = {0.12, 4.12, 0.12};
	cout << "Vecteur d" << endl << d << endl << endl;
	
	cout << "Resolution de aX = d par MinResSolve" << endl;
	vector<double> res = MinresSolve(a, d);
	cout << "Vecteur X = " << res <<endl;
	
	// Question 2/3
	
	cout << endl << "Test de CrossApproximation et des FrobeniusNorm (Question 2 et 3)" << endl << endl;
	
	DenseMatrix b(3,3);
	
	b(0,0) = 60;
	b(0,1) = 79;
	b(0,2) = 33;
	b(1,0) = 42;
	b(1,1) = 68;
	b(1,2) = 19;
	b(2,0) = 13;
	b(2,1) = 58;
	b(2,2) = 41;
	
	FredholmMatrix appro = CrossApproximation(b, 3);
	
	cout << "CrossApproximation de b" << endl << appro << endl;
	
	double norm = FrobeniusNorm(b);
	
	cout << "Frobenius Norme de b =\t" << norm << endl;
	
	double norm2 = FrobeniusNorm(appro);
	
	cout << "Frobenius Norme de son approximation =\t" << norm2 <<endl<<endl;
	
	// Question 4 
	
	cout << "Test de la fonction qui cree la matrice de la question 4 et de PlotGraph (Question 4)" << endl;
	
	DenseMatrix b4 = create_mat_q4(3);
	
	cout << "Matrice B pour n = 3" << endl << b4 << endl ;
	
	cout << "Graphe de la question 4 pour n = 100" << endl;
	PlotGraph(100);
	
	// Question 5/6
	
	cout << endl << "Test de la fonction LoadDenseMatrix et PlotGraph (Question 5 et 6)" << endl;
	
	PlotGraph("matrix.txt");
	
	return 0;
}