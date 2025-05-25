#ifndef MAP_MATRIX
#define MAP_MATRIX
#include <iostream>
#include <vector>
#include <map>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>

using namespace std;
template <typename value_t> 

class MapMatrix{
	
	int nr, nc; // Nombre de row, nombre de column
	
	typedef std::tuple<int,int> NxN;  // On crée un tuple de 2 int
	
	map<NxN, value_t> data;  // Système clés / valeurs avec clés uniques
	
public : 
	
	// Constructor
	
	MapMatrix(const int& nr0=0, const int& nc0=0):
		nr(nr0), nc(nc0) {}; // On fait rien de data
	
	// Copy constructor
	
	MapMatrix(const MapMatrix& m):
		nr(m.nr), nc(m.nc), data(m.data) {};
		
	// Operator = qui copie
	
	MapMatrix& operator=(const MapMatrix& m){
		nr = m.nr;
		nc = m.nc;
		data = m.data;
		return *this;
	}
	
	// Fonction insert qui insère un nouvel élément dans data (l'ajoute au precedent)

	void insert(const int& j, const int& k, const value_t& v){
		tuple<int, int> t = make_tuple(j,k);  // NxN t = {j,k}
		data[t] += v;
	}
	
	// Fonction insert qui insère un nouvel élément dans data le remplace

	void insert2(const int& j, const int& k, const value_t& v){
		tuple<int, int> t = make_tuple(j,k);  // NxN t = {j,k}
		data[t] = v;
	}
	
	// Fonction pour regarder les valeurs de la map
	
	double operator()(const int& j, const int& k) const {
		assert(j<nr && k<nc);
		tuple<int, int> t = {j,k};
		
		if(auto search = data.find(t); search != data.end()){  // Si on trouve t, et qu'on est pas arrivé à la fin
			return search->second;  // On prend la valeur (pas la clé)
		} else {
			return 0.;
		}	
	}

	// Operator + entre deux MapMatrix
	
	MapMatrix operator+(const MapMatrix& m){
		assert(m.nc == nc && m.nr == nr);
		MapMatrix& n = *this;
		MapMatrix c(nr, nc);
		for (int j = 0; j<nr; j++){
			for (int k = 0; k<nc; k++){
				NxN t = {j,k};
				c.data[t] = n(j,k)+m(j,k);
			}
		}
		return c;
	}
	
	// Operator += entre deux MapMatrix
	
	MapMatrix& operator+=(const MapMatrix& m){
		assert(m.nc == nc && m.nr == nr);
		for (int j = 0; j<nr; j++){
			for (int k = 0; k<nc; k++){
				NxN t = {j,k};
				data[t] += m(j,k);
			}
		}
		return *this; 
	}
	
	// Operator - entre deux MapMatrix
	
	MapMatrix operator-(const MapMatrix& m){
		assert(m.nc == nc && m.nr == nr);
		MapMatrix& n = *this;
		MapMatrix c(nr, nc);
		for (int j = 0; j<nr; j++){
			for (int k = 0; k<nc; k++){
				NxN t = {j,k};
				c.data[t] = n(j,k)-m(j,k);
			}
		}
		return c;
	}
	
	// Operator -= entre deux MapMatrix
	
	MapMatrix& operator-=(const MapMatrix& m){
		assert(m.nc == nc && m.nr == nr);
		for (int j = 0; j<nr; j++){
			for (int k = 0; k<nc; k++){
				NxN t = {j,k};
				data[t] -= m(j,k);
			}
		}
		return *this; 
	}
	
	// Operator * entre deux MapMatrix
	
	MapMatrix operator*(const MapMatrix& m){
		assert(nc == m.nr);
		MapMatrix& n = *this;
		MapMatrix c(nr, m.nc);

		for (int i = 0; i<nr; i++){
			for (int j = 0; j<m.nc; j++){
				for(int k = 0; k<nc; k++){
					NxN t = {i,j};
					c.data[t] += n(i,k)*m(k,j);
				}
			}
		}
		return c;
	
	}
	
	// Operator * entre une MapMatrix et un vecteur
	
	vector<value_t> operator*(const vector<value_t>& v){ // Const veut dire que l'on peut utiliser la fonction même si c'est constant
		
		assert(v.size() == nc);
		MapMatrix& n = *this;
		vector<value_t> w(nr);
		
		for (int i = 0; i<nr; i++){
			for (int k = 0; k<nc; k++){
				w[i]+=n(i,k)*v[k];
			}
		}
		return w;
	}
	
	// Operator * entre une MapMatrix et un vecteur
	
	vector<value_t> operator*(const vector<value_t>& v) const{ // Const veut dire que l'on peut utiliser la fonction même si c'est constant
		
		assert(v.size() == nc);
		const MapMatrix& n = *this;
		vector<value_t> w(nr);
		
		for (int i = 0; i<nr; i++){
			for (int k = 0; k<nc; k++){
				w[i]+=n(i,k)*v[k];
			}
		}
		return w;
	}
	
	// Operator *= entre deux MapMatrix
	
	MapMatrix& operator*=(const MapMatrix& m){
		
		MapMatrix& n = *this;
		n = n*m;
		return n;
	}
	
	// Operator *= entre une matrice et un double
	
	MapMatrix& operator*=(const value_t& a){
		
		MapMatrix& n = *this;
		for (auto it=data.begin(); it!=data.end(); it++){
			it->second*=a;
		}
		return n;
	}
	
	// Fonction AMI qui donne le nombre de lignes
	friend const int& NbRow(const MapMatrix& m){
		return m.nr;
	}
	
	// Fonction AMI qui donne le nombre de colonnes
	friend const int& NbCol(const MapMatrix& m){
		return m.nc;
	}
	
	// Operator * entre un double et une DenseMatrix

	friend MapMatrix operator*(const value_t& a, const MapMatrix& m){
	
		MapMatrix c = m;
		
		for(auto it = c.data.begin(); it!=c.data.end(); it++){
			it->second *= a;
		}
		return c;
	}

	// Operateur << pour afficher
	
	friend ostream& operator<<(std::ostream& o,const MapMatrix& m){
		for(int j=0; j<NbRow(m); ++j){
			for(int k=0; k<NbCol(m); ++k){
				o << m(j,k) << "\t";
			}
			o << "\n";
		}
		return o;
	}
	
	friend void Write(const string&, const MapMatrix<double>&);

};	

// Fonction qui transcrit un fichier en MapMatrix

MapMatrix<double> LoadMapMatrix(const string& filename){
	
	ifstream fic;
	fic.open(filename);
	
	if(fic)
	{
		int nr, nc;
		fic >> nr;
		fic >> nc;
		
		MapMatrix<double> m(nr,nc);
		
		// cout << "Nombre de ligne : " << nr << endl;
		// cout << "Nombre de colonne : " << nc << endl;
		
		string ligne;
		int j, k;
		double v;
		
		while(getline(fic, ligne))
		{
			fic >> j >> k >> v;
			m.insert(j,k,v);
			// cout << j << "\t" << k << "\t" << v << endl;
		
		}
		fic.close();
		return(m);
		
	} else {
		cout << "Erreur" << endl;
		MapMatrix<double> m;
		fic.close();
		return(m);
	}
}

// Fonction qui transcrit une MapMatrix en fichier (avec un nom prédéfinit)

void Write(const string& filename, const MapMatrix<double>& m){
	
	ofstream fic;
	fic.open(filename);
	
	if(fic)
	{
		fic << NbRow(m) << "\t" << NbCol(m) << endl;
		
		for(auto it = m.data.begin(); it!=m.data.end(); it++){
			fic << get<0>(it->first) << "\t" << get<1>(it->first) << "\t" << it->second << endl;
		}
		
		fic.close();
		
	} else {
		cout << "Erreurr" << endl;
	}
	
}


#endif