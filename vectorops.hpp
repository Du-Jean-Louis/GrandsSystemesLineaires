#ifndef VECTOR_OPS
#define VECTOR_OPS

#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>

using namespace std;

ostream& operator<<(ostream& o, const vector<double>& u){
	for(size_t j=0; j<u.size(); j++){
		o << u[j] << "\t";
	}
	return o;
}

void print(vector<double>& vect){
	for (size_t j=0; j<vect.size(); j++){
		cout << vect[j] << endl;
	}		
}

void print(vector<int>& vect){
	for (size_t j=0; j<vect.size(); j++){
		cout << vect[j] << endl;
	}		
} 

vector<double>& operator+=(vector<double>& lhs, const vector<double>& rhs){
	
	assert(lhs.size()==rhs.size());
	
	for (size_t j = 0; j<lhs.size(); j++){
		lhs[j] += rhs[j];
	}
	return lhs;
}


vector<double>& operator-=(vector<double>& lhs, const vector<double>& rhs){
	
	assert(lhs.size()==rhs.size());
	
	for (size_t j = 0; j<lhs.size(); j++){
		lhs[j] -= rhs[j];
	}
	return lhs;
}


vector<double> operator+(const vector<double>& lhs, const vector<double>& rhs){

	assert(lhs.size()==rhs.size());
	vector<double> res(lhs.size(),0.);
	
	for (size_t j = 0; j<lhs.size();j++){
		res[j] = lhs[j]+rhs[j];
	}
	return res;
}

vector<double> operator-(const vector<double>& lhs, const vector<double>& rhs){

	assert(lhs.size()==rhs.size());
	vector<double> res(lhs.size(),0.);
	
	for (size_t j = 0; j<lhs.size();j++){
		res[j] = lhs[j]-rhs[j];
	}
	return res;
}		

double operator,(const vector<double>& lhs, const vector<double>& rhs){

	assert(lhs.size()==rhs.size());
	double res = 0;
	
	for (size_t j = 0; j<lhs.size();j++){
		res += lhs[j]*rhs[j];
	}
	return res;
}		

double Norm(const vector<double>& vect){

	double norm = 0.0;
	
	for (size_t j = 0; j<vect.size();j++){
	norm += vect[j]*vect[j];            // Comment on fait un carré      
	}
	return sqrt(norm);
}		

vector<double>& operator*=(vector<double>& vect, double multiplicateur){
	
	for (size_t j=0; j<vect.size(); j++){
	vect[j]*=multiplicateur;
	}
	return vect;
}

// vector<double>& operator*(double multiplicateur, vector<double>& vect){
	
	// for (size_t j=0; j<vect.size(); j++){
	// vect[j]*=multiplicateur;
	// }
	// cout << "used" << endl;
	// return vect;
// }

// La meme chose avec des const pour le TP4

// vector<double> operator*(const double multiplicateur, const vector<double>& vect){
	
	// vector<double> res(vect.size());
	// for (size_t j=0; j<vect.size(); j++){
	// res[j] = vect[j]*multiplicateur;
	// }
	// return res;	
// }

// Version pour le projet

vector<double> operator*(double multiplicateur, const vector<double>& vect){
	
	vector<double> res(vect.size());
	for (size_t j=0; j<vect.size(); j++){
		res[j] = vect[j]*multiplicateur;
	}
	return res;
}

double norme2(const vector<double>& vect){
	
	double sum = 0;
	
	for (size_t j = 0; j<vect.size(); j++){
		sum += vect[j]*vect[j];
	}		
	return sqrt(sum);
}

#endif