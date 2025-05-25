#ifndef FREDHOLM
#define FREDHOLM

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include "mapmatrix.hpp"
#include "vectorops.hpp"
#include "densematrix.hpp"

using namespace std;

class FredholmMatrix{
	
	int n;     // size of the matrix
	MapMatrix<double> diag;  // MapMatrix qui représente D
	vector < vector <double> > lru;   // représente les uj
	vector < vector <double> > lrv;   // représente les vj
	
	public :
	
	// Constructor qui initialise avec les données

	FredholmMatrix(const int& n0 = 0):
		n(n0), diag(n0,n0) {};
		
	// Copy constructor
	
	FredholmMatrix(const FredholmMatrix& a):
		n(a.n), diag(a.diag), lru(a.lru), lrv(a.lrv) {};
		
	// Operator = 
	
	FredholmMatrix operator=(const FredholmMatrix& a){
		
		n = a.n;
		diag = a.diag;
		lru = a.lru;
		lrv = a.lrv;
		return *this;
		
	}
		
	// Member function insert to increment D
	
	void insert(const int& j, const int& k, const double& val){  
		
		assert(j <= n && k <= n && j>0 && k > 0); // j, k compris entre 1 et n
		
		diag.insert(j-1,k-1,val);  // La fonction insert des MapMatrix commence le compte d'indice a 0
		
	}
	
	// Member function insert to increment lru and lrv
	
	void insert(const vector<double>& u, const vector<double>& v){
		
		assert(u.size() == n && v.size() == n);
		
		lru.push_back(u);   // On rajoute un vecteur au vecteur lru, meme chose pour lrv
		lrv.push_back(v);
		
	}
	
	// Operator *
	
	vector<double> operator*(const vector<double>& vect){
		
		assert(vect.size() == n);
		
		vector<double> res = diag*vect;  // Par linearite du produit on peut decomposer le produit sommes
	                                     // Ainsi on ajoute tout d'abord le produit de diag par vect
		
		for (int j = 0; j<lru.size(); j++){
			double coeff = (lrv[j],vect); // Ensuite au lieu de creer une matrice pour ensuite faire le produit matrice vecteur
			res += coeff*lru[j];          // On voit que c'est plus simple de faire vj scalaire vect qui donne un double qu'on peut ensuite multiplie a uj
		}
		return res;
	}
	
	// Operator * version const
	
	vector<double> operator*(const vector<double>& vect) const{
		
		assert(vect.size() == n);
		
		vector<double> res = diag*vect;
		for (int j = 0; j<lru.size(); j++){
			double coeff = (lrv[j],vect);  // , correspond au produit scalaire ici
			res += coeff*lru[j];
		}
		return res;
	}
	
	// Operator () to acess the value
	
	double operator()(const int& j, const int& k){
		
		assert(j <= n && k <= n && j>0 && k > 0);
		
		FredholmMatrix& m = *this;
		
		vector<double> ek(n);  // Pour lire la valeur d'indice (j,k) on remarque qu'on peut regarder la j-eme composante
		ek[k-1] = 1;           // du produit par le vecteur de la base canonique ek, on creer donc ce vecteur ek
		
		double res = (m*ek)[j-1];
		
		return res;
		
	}
	
	// Version const
	
	double operator()(const int& j, const int& k) const{
		
		assert(j <= n && k <= n && j>0 && k > 0);
		
		const FredholmMatrix& m = *this;
		FredholmMatrix copie(m);
		
		vector<double> ek(n);
		ek[k-1] = 1;
		
		double res = (copie*ek)[j-1];
		
		return res;
		
	}
	
	// Fonction friend pour recuperer la taille n d'une FredholmMatrix
	
	friend const int& Recup_n(const FredholmMatrix& m){
		return m.n;
	}
	
	// Fonction friend pour print
	
	friend ostream& operator<<(ostream& o,const FredholmMatrix& m){
		
		int n_m = Recup_n(m);
		
		for(int j=1; j<=n_m; ++j){
			for(int k=1; k<=n_m; ++k){
				o << m(j,k) << "\t";
			}
			o << "\n";
		}
		
		return o;
	}
	
	// Fonction friend qui resoud Ax = b
	
	friend vector<double> MinresSolve(const FredholmMatrix& A, const vector<double>& b){
		
		assert(Recup_n(A) == b.size());
		
		vector<double> x(Recup_n(A));
		vector<double> r;
		double alpha;
		int iterator = 0;
		const int maxit = 100;
		const double tol = 1e-6;
		
		while (iterator < maxit){
			
			r = b - A*x;
			alpha = (A*r,r)/(Norm(A*r)*Norm(A*r));
			x += alpha*r;
			iterator++;

			if (Norm(r) < tol) {
				cout << "Solution convenable trouvee" << endl;
				break;
			}
		}
		return x;
	}
};

typedef tuple<int,int> NxN;

// Fonction qui renvoie l'indice (j,k) pour lequel on atteint le max des valeurs en valeur absolue

NxN argmax(const DenseMatrix& B){
	
	int j_max = 0;
	int k_max = 0;
	
	for (int j = 0; j < NbRow(B); j++){
		for (int k = 0; k < NbCol(B); k++){
			if (abs(B(j_max,k_max)) < abs(B(j,k))){  // On cherche l'indice ou est atteint la plus grande valeur de B en valeur absolue par comparaison succesive
				j_max = j;
				k_max = k;
			}	
		}
	}
	
	NxN j_k_max = make_tuple(j_max+1, k_max+1);
	
	return j_k_max;
}

// Fonction qui transforme la j-eme ligne en vecteur colonne

vector<double> Row(const DenseMatrix& B, const int& j){
	
	assert(j-1 < NbRow(B) && j > 0);
	
	vector<double> row(NbCol(B));
	for (int i = 0; i < NbCol(B); i++){
		row[i] = B(j-1,i);  // On recopie la j-eme ligne dans un vecteur colonne
	}
	return row;
}

// Fonction qui transforme la k-eme colonne en vecteur colonne

vector<double> Col(const DenseMatrix& B, const int& k){
	
	assert(k-1 < NbCol(B) && k > 0);
	
	vector<double> col(NbRow(B));
	for (int i = 0; i < NbRow(B); i++){
		col[i] = B(i,k-1);  // On recopie la k-eme colonne dans un vecteur colonne
	}
	return col;
}

// Fonction qui donne la matrice forme par le produit de deux vecteur de bonnes tailles

DenseMatrix mat_prod(const vector<double> u, const vector<double> v){
	
	assert(u.size() == v.size());
	int n = u.size();
	DenseMatrix mat_prod(n,n);
	
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			mat_prod(i,j) = u[i]*v[j]; 
		}
	}
	return mat_prod;
}

// Fonction qui effectue l'algorithme de CrossApproximation

FredholmMatrix CrossApproximation(const DenseMatrix& B, const int& r){
	
	assert(NbRow(B) == NbCol(B));
	FredholmMatrix A(NbRow(B));
	
	DenseMatrix R = B;
	DenseMatrix Btilde(NbRow(B), NbCol(B));
	NxN j_k_star = make_tuple(0,0);
	
	for (int p = 0; p < r; p++){
		
		j_k_star = argmax(R);  // On trouve le couple (j,k) donnant la plus grande valeur en valeur absolue
		int jstar = get<0>(j_k_star);   // On les depack
		int kstar = get<1>(j_k_star); 
		
		if (R(jstar-1, kstar-1) == 0){  // Si toute les valeurs sont nul, on arrete l'algorithme
			break;
		}
		
		vector<double> up = (1/(R(jstar-1,kstar-1)))*Col(R,kstar);  // On calcule up et vp
		vector<double> vp = Row(R,jstar);
		//cout << "up\n" << up << endl << "vp\n" << vp << endl;
		A.insert(up, vp);                                           // Puis on les rentres dans lru et lrv de la FredholmMatrix
		DenseMatrix add = mat_prod(up,vp);
		R -= add;
		Btilde += add;	
		//cout << "add pour " << p << endl << add << endl;
		//cout << "R pour " << p << endl << R << endl;
		//cout << "Btilda pour " << p << endl << Btilda << endl << endl;
	}
	
	cout << "R final " << endl << R << endl;
	cout << "Btilde final" << endl << Btilde << endl << endl;
	
	return A;
}

// Fonction qui calcule la frobenius norm d'une DenseMatrix

double FrobeniusNorm(const DenseMatrix& B){
	
	assert(NbRow(B) == NbCol(B));
	
	double res = 0;
	
	for (int j = 0; j < NbRow(B); j++){
		for(int k = 0; k < NbCol(B); k++){
			res += B(j,k)*B(j,k);
		}
	}
	return sqrt(res);
}

// Fonction qui calcule la frobenius norm d'une FredholmMatrix 

double FrobeniusNorm(const FredholmMatrix& A){
	
	double res = 0;
	
	for (int j = 1; j <= Recup_n(A); j++){
		for(int k = 1; k <= Recup_n(A); k++){
			res += A(j, k)*A(j,k);
		}
	}
	return sqrt(res);
}

// Fonction qui produit la matrice de la question 4 etant donne une taille n

DenseMatrix create_mat_q4(const size_t& n){
	
	DenseMatrix B(n,n);
	
	for (int j = 1; j<=n; j++){
		for (int k = 1; k<=n; k++){
			
			double inter = -(j-k)*(j-k);   // On a separer la formule car l'ordinateur n'arrive pas a interpreter l'entierete
			inter /= n*n;
			B(j-1,k-1) = exp(inter);
			
		}
	}
	
	return B;
	
}

// Fonction qui écrit les datas dans un fichier texte bis

void create_datatxt4(const size_t& n, const DenseMatrix& B){   // On reprend l'algorithme de CrossApproximation et on la modifie pour calculer la norme du reste a chaque etape
	
	ofstream fic("dataq4.txt");
	
	if(fic){   // Et on rentre les valeurs dans un fichier texte accompagne de la taille de lru et lrv correspondante
		
		double val;
		
		DenseMatrix R = B;
		DenseMatrix Btilda(NbRow(B), NbCol(B));
		NxN j_k_star = make_tuple(0,0);
		
		for (int p = 0; p < n; p++){
			
			j_k_star = argmax(R);
			int jstar = get<0>(j_k_star); 
			int kstar = get<1>(j_k_star); 
			
			if (R(jstar-1, kstar-1) == 0){
				break;
			}
			
			vector<double> up = (1/(R(jstar-1,kstar-1)))*Col(R,kstar);
			vector<double> vp = Row(R,jstar);
			
			DenseMatrix add = mat_prod(up,vp);
			R -= add;   // C'est le reste qui nous interesse ici
			Btilda += add;	
			
			val = FrobeniusNorm(R);
			// cout << p+1 << "\t" << val << endl;
			fic << p+1 << "\t" << val << endl;
		}
		fic.close();
	}
}

// Fonction qui plot un graphe, il produit un pdf

void PlotGraph(const size_t& n){
	
	DenseMatrix B = create_mat_q4(n);  // On cree la matrice demande
	
	create_datatxt4(n, B);   // On cree le fichier texte demande
	
	// Commande pour générer le pdf
    const char* gnuplotCommand = "gnuplot -e \"set term pdf; set output 'graph4.pdf'; \
									set logscale xy; plot 'dataq4.txt' with lines\"";

    int exec = system(gnuplotCommand);  // On execute la commande, si exec = 0 alors ca a bien ete executé

    // On verifie que la commande a bien ete execute
    if (exec != 0) {
        cout << "Erreur" << endl;
    } else {
        cout << "PDF cree!" << std::endl;
		cout << "Vous pouvez consulter la courbe en ouvrant graph4.pdf" << endl;
    }
	
}

// Fonction qui lit un fichier texte et le transforme en DenseMatrix

DenseMatrix LoadDenseMatrix(const string& filename){
	
	ifstream fic;
	fic.open(filename);
	
	if(fic)
	{
		string ligne;
		getline(fic, ligne);  // On passe la ligne avec #SIZE
		
		int nr, nc;
		fic >> nr;
		fic >> nc;
		
		DenseMatrix m(nr,nc);
		
		// cout << "Nombre de ligne : " << nr << endl;
		// cout << "Nombre de colonne : " << nc << endl;
		
		getline(fic, ligne);  // On passe la ligne apres avoir lu nr et nc
		getline(fic, ligne);  // On passe la ligne avec #DATA
		
		double value;
		int j, k;
		int compteur = 0;
		
		for (int i = 0; i < nr*nc; i++){
			
			j = i/nc;  // Numero de la ligne
			k = i%nc;  // Numero de la colonne
			fic >> value;
			m(j,k) = value;  // On rentre la valeur lu dans la matrice
			compteur ++;
			// cout << j << "\t" << k << "\t" << value << endl;
			if (compteur == 10){  // Puisque les lignes ont 10 valeurs, des qu'on arrive a la 10eme valeur on passe a la ligne suivante
				getline(fic, ligne);
			}
		}
		fic.close();
		return m;
	
	} else {
		cout << "Erreur" << endl;
		DenseMatrix m;
		fic.close();
		return(m);
	}
}

// Fonction qui écrit les datas dans un fichier texte bis

void create_datatxt6(const size_t& n, const DenseMatrix& B){  // Meme chose que l'autre fonction create_data, juste nom du fichier texte genere different
	
	ofstream fic("dataq6.txt");
	
	if(fic){
		
		double val;
		
		DenseMatrix R = B;
		DenseMatrix Btilda(NbRow(B), NbCol(B));
		NxN j_k_star = make_tuple(0,0);
		
		for (int p = 0; p < n; p++){
			
			j_k_star = argmax(R);
			int jstar = get<0>(j_k_star); 
			int kstar = get<1>(j_k_star); 
			
			if (R(jstar-1, kstar-1) == 0){
				break;
			}
			
			vector<double> up = (1/(R(jstar-1,kstar-1)))*Col(R,kstar);
			vector<double> vp = Row(R,jstar);
			
			DenseMatrix add = mat_prod(up,vp);
			R -= add;
			Btilda += add;	
			
			val = FrobeniusNorm(R);
			// cout << p+1 << "\t" << val << endl;
			fic << p+1 << "\t" << val << endl;
		}
		fic.close();
	}
}

// Fonction finale

void PlotGraph(const string& filename){
	
	DenseMatrix B = LoadDenseMatrix(filename);
	
	int n = NbRow(B); // Matrice carree ici
	
	create_datatxt6(n, B);
	
	// Command pour générer le pdf
    const char* gnuplotCommand = "gnuplot -e \"set term pdf; set output 'graph6.pdf'; \
									set logscale xy; plot 'dataq6.txt' with lines\"";

    // On l'execute
    int exec = system(gnuplotCommand);  // Si exec = 0, bien executé

    // On verifie que la commande a bien ete execute
    if (exec != 0) {
        cerr << "Erreur" << endl;
    } else {
        cout << "PDF cree!" << std::endl;
		cout << "Vous pouvez consulter la courbe en ouvrant graph6.pdf" << endl <<endl;
    }
}	


#endif