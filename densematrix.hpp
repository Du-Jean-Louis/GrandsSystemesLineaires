#ifndef DENSEMATRIX_OPS
#define DENSEMATRIX_OPS
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include "vectorops.hpp"

using namespace std;

class DenseMatrix{
	
	int nr, nc;
	vector<double> data; // On va les stocker dans un vecteur lignes par lignes 
	
	public:
	
	// Constructor qui initialise la DenseMatrix
	DenseMatrix(const int& nr0 = 0, const int& nc0 = 0): 
		nr(nr0), nc(nc0), data(nc0*nr0, 0.) {};
	
	// Copy constructor
	DenseMatrix(const DenseMatrix& m):
		nr(m.nr), nc(m.nc), data(m.data) {};
	
	// Operator = qui copie la matrice
	DenseMatrix& 
	operator=(const DenseMatrix& m){
		nr = m.nr;
		nc = m.nc;
		data = m.data;
		return *this;     // On est au sein de la class DenseMatrix, du coup nc, nr et data existe déjà, on a pas à les nommer
	}                     // this désigne l'adresse de la DenseMatrix, *this désigne donc la Densematrix elle-même

	// Operator , qui retrouve le coeff à la jème ligne et kème colonne, version non const, on peut la modifier
	double& operator()(const int& j, const int& k){
		assert(j<nr && k<nc);  // On prend j qui démarre à la 0ème ligne et k à la 0ème colonne
		return data[j*nc+k];
	}
	
	// Operator , qui retrouve le coeff à la jème ligne et kème colonne, version const, accessible en mode lecture seulement
	double operator()(const int& j, const int& k) const {
		assert(j<nr && k<nc);  // On prend j qui démarre à la 0ème ligne et k à la 0ème colonne
		return data[j*nc+k];
	}
	
	// Operator +=
	DenseMatrix& operator+=(const DenseMatrix& m){
		assert(nc == m.nc && nr == m.nr);
		DenseMatrix& n = *this; 
		for (int i = 0; i<nr; i++){
			for (int j = 0; j<nc; j++){
				n(i,j) += m(i,j);
			}
		}
		return n;
	}
	
	// Operator -=
	DenseMatrix& operator-=(const DenseMatrix& m){
		assert(nc == m.nc && nr == m.nr);
		DenseMatrix& n = *this; 
		for (int i = 0; i<nr; i++){
			for (int j = 0; j<nc; j++){
				n(i,j) -= m(i,j);
			}
		}
		return n;
	}
	
	// Operator +
	DenseMatrix operator+(const DenseMatrix& m){
		assert(nc == m.nc && nr == m.nr);
		DenseMatrix c(m.nr,m.nc);
		for (int i = 0; i<nr; i++){
			for (int j = 0; j<nc; j++){
				c.data[i*nc+j] = data[i*nc+j] + m.data[i*nc+j];
			}
		}
		return c;
	}
	
	// Operator -
	DenseMatrix operator-(const DenseMatrix& m) const{
		assert(nc == m.nc && nr == m.nr);
		DenseMatrix c(m.nr,m.nc);
		for (int i = 0; i<nr; i++){
			for (int j = 0; j<nc; j++){
				c.data[i*nc+j] = data[i*nc+j] - m.data[i*nc+j];
			}
		}
		return c;
	}
	
	// Operator -
	DenseMatrix operator-(const DenseMatrix& m) {
		assert(nc == m.nc && nr == m.nr);
		DenseMatrix c(m.nr,m.nc);
		for (int i = 0; i<nr; i++){
			for (int j = 0; j<nc; j++){
				c.data[i*nc+j] = data[i*nc+j] - m.data[i*nc+j];
			}
		}
		return c;
	}
	
	// Operator * entre DenseMatrix
	DenseMatrix operator*(const DenseMatrix& m){
		assert(nc == m.nr);
		DenseMatrix& n = *this;
		DenseMatrix c(nr, m.nc);

		for (int i = 0; i<nr; i++){
			for (int j = 0; j<m.nc; j++){
				for(int k = 0; k<nc; k++){
					c(i,j) += n(i,k)*m(k,j);
				}
			}
		}
		return c;
	}
	
	// Operator * entre une DenseMatrix et un vector
	
	vector<double> operator*(const vector<double>& vect){
		
		assert(vect.size() == nc);
		
		DenseMatrix& n = *this;
		vector<double> res(nr, 0.);
		
		for (int j = 0; j<nr ; j++){
			for (int i = 0; i<vect.size(); i++){
				res[j] +=  n(j,i)*vect[i];
			}
		}
		return res;
	}
	
	// Operator *= qui multiplie une DenseMatrix par un double
	
	DenseMatrix& operator*=(const double& a){
		DenseMatrix& n = *this;
		
		for (int i = 0; i<nr; i++){
			for (int j = 0; j<nc; j++){
				n(i,j) *= a;
			}
		}
		return n;
	}
	
	// Operator *= qui multiplie une Densematrix par une autre
	
	DenseMatrix& operator*=(DenseMatrix& m){
		
		DenseMatrix& n = *this;
		n = n*m;
		return n;
	}
	
	// Fonction AMI qui donne le nombre de lignes
	friend const int& NbRow(const DenseMatrix& m){
		return m.nr;
	}
	
	// Fonction AMI qui donne le nombre de colonnes
	friend const int& NbCol(const DenseMatrix& m){
		return m.nc;
	}

};

// Operator * entre un double et une DenseMatrix

DenseMatrix operator*(const double& a, DenseMatrix& m){
	
	DenseMatrix c(NbRow(m), NbCol(m));
	
	for (int i = 0; i<NbRow(m); i++){
		for (int j = 0; j<NbCol(m); j++){
			c(i,j) = a*m(i,j);
		}
	}
	return c;
}

// Operator output stream <<

std::ostream& operator<<(std::ostream& o, const DenseMatrix& m){
	for(std::size_t j=0; j<NbRow(m); ++j){
		for(std::size_t k=0; k<NbCol(m); ++k){
			o << m(j,k) << "\t";
		}
		o << "\n";
	}
  return o;
}

// Donner la factorisation LU de la densematrix 

void LUFactorize(DenseMatrix& mat){
	
	for (int i = 0; i<NbRow(mat)-1; i++){  // Boucle pour tout les pivots
		
		assert(mat(i,i) != 0);
		
		for (int j = i+1; j<NbRow(mat); j++){  // Boucle sur les lignes en dessous
			
			mat(j,i) = mat(j,i)/mat(i,i); // On donne la valeur dans L
			
			for (int k = i+1; k<NbCol(mat); k++){   // Boucle sur les colonnes à partir de i
				
				mat(j,k) = mat(j,k) - mat(j,i)*mat(i,k);
				
			}				
		}
	}
	
}

// Donner la factorisation LU de la densematrix (Version on renvoie la matrice)

DenseMatrix LUFactorizebis(DenseMatrix& mat){
	
	for (int i = 0; i<NbRow(mat)-1; i++){  // Boucle pour tout les pivots
		
		assert(mat(i,i) != 0);
		
		for (int j = i+1; j<NbRow(mat); j++){  // Boucle sur les lignes en dessous
			
			mat(j,i) = mat(j,i)/mat(i,i); // On donne la valeur dans L
			
			for (int k = i+1; k<NbCol(mat); k++){   // Boucle sur les colonnes à partir de i
				
				mat(j,k) = mat(j,k) - mat(j,i)*mat(i,k);
				
			}				
		}
	}
	return(mat);
}

// Donner la LU Factorisation avec pivotement

DenseMatrix LUFactorizePivotRow(DenseMatrix& mat, vector<int>& pivot){
	
	DenseMatrix copie = mat;
	
	for (int i = 0; i<NbRow(mat)-1; i++){  // Boucle pour tout les pivots
		
		// On parcourt la colonne pour chercher le plus gros coeff
		
		int pv = i;  // On initialise à pivot i
		for (int p = i+1; p<NbRow(mat); p++){
			if(abs(mat(pivot[p],i)) > abs(mat(pivot[pv],i))){  // Si la valeur est plus grand, on la garde en mémoire
				pv = p;
			}
		}
		// A la fin on prend le plus grand
		int temp = pivot[i];
		pivot[i] = pivot[pv];
		pivot[pv] = temp;
		
		assert(mat(pivot[i],i) != 0); // On verifie que le plus grand n'est pas zeros
		
		for (int j = i+1; j<NbRow(mat); j++){  // Boucle sur les lignes en dessous
			
			mat(pivot[j],i) = mat(pivot[j],i)/mat(pivot[i],i); // On donne la valeur dans L
			
			for (int k = i+1; k<NbCol(mat); k++){   // Boucle sur les colonnes à partir de i
				
				mat(pivot[j],k) = mat(pivot[j],k) - mat(pivot[j],i)*mat(pivot[i],k);
				
			}
		}
	}
	
	for (int i = 0;  i<NbRow(copie); i++){
		for (int j = 0; j<NbCol(copie); j++){
			copie(i, j) = mat(pivot[i], j);
		}
	}
	
	mat = copie;
	return(copie);
}

DenseMatrix LUFactorizePivotCol(DenseMatrix& mat, vector<int>& pivot){
	
	DenseMatrix copie = mat;
	
	for (int i = 0; i<NbCol(mat)-1; i++){  // Boucle pour tout les pivots
		
		// On parcourt la colonne pour chercher le plus gros coeff
		
		int pv = i;  // On initialise à pivot i
		for (int p = i+1; p<NbCol(mat); p++){
			if(abs(mat(i,pivot[p])) > abs(mat(i,pivot[pv]))){  // Si la valeur est plus grand, on la garde en mémoire
				pv = p;
			}
		}
		// A la fin on prend le plus grand
		int temp = pivot[i];
		pivot[i] = pivot[pv];
		pivot[pv] = temp;
		
		assert(mat(i,pivot[i]) != 0); // On verifie que le plus grand n'est pas zeros
		
		for (int j = i+1; j<NbRow(mat); j++){  // Boucle sur les lignes en dessous
			
			mat(j,pivot[i]) = mat(j,pivot[i])/mat(i,pivot[i]); // On donne la valeur dans L
			
			for (int k = i+1; k<NbCol(mat); k++){   // Boucle sur les colonnes à partir de i
				
				mat(j,pivot[k]) = mat(j,pivot[k]) - mat(j,pivot[i])*mat(i,pivot[k]);
				
			}
		}
	}
	
	for (int i = 0;  i<NbRow(copie); i++){
		for (int j = 0; j<NbCol(copie); j++){
			copie(i, j) = mat(i, pivot [j]);
		}
	}
	
	mat = copie;
	return(copie);
}


/*
void LUFactorizePivot(DenseMatrix& mat, vector<int>& perm){
	
	int N = NbRow(mat);
	
	for (int p = 0; p<N-1; p++){
		
		// On parcourt la colonne pour chercher le plus gros coeff
		
		double pivot = mat(perm[p],p);
		for (int j = p+1; p<N; p++){
			if(abs(mat(perm[j],p)) > abs(pivot)){  // Si la valeur est plus grand, on la garde en mémoire
				pivot = mat(perm[j],p);
				swap(perm[j],perm[p]);
			}
		}
		
		assert(mat(perm[p],p) != 0); 
		cout << "o" << endl;
		for (int j = p+1; j<N; j++){  
		
			mat(perm[j],p) /= pivot; 
			
			for (int k = p+1; k<N; k++){   // Boucle sur les colonnes à partir de i
				
				mat(perm[j],k) -= mat(perm[j],p)*mat(perm[p],k);
				
			}
			
		cout << mat << endl ;
		
		}
	}
	
}
*/


#endif