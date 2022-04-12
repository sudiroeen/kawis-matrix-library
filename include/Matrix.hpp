/*
	Copyright:
				Sudiro
						[at] SudiroEEN@gmail.com

*/
#pragma once

#include <iostream>
#include <sstream>
#include <initializer_list>
#include <unordered_map>

#include <cmath>
#include <limits>

#include <random>

#define EPSILON 1e-5

using namespace std;

typedef std::pair<int, int> size_m;

class Matrix{
public:
	int row;
	int col;
	int n_rank_of_this_matrix = 0;
private:
	double**  entry;
public:
	////////////////  CONSTRUCTOR
	Matrix();
	// ~Matrix(){                     // DESTRUCTOR -> SUMBER MASALAH ERROR HADEH
	// 	for(int r=0; r<row; r++)
	// 		delete[] entry[r];
	// 	delete[] entry;
	// }

	Matrix(int row_, int col_);
	Matrix(size_m sm);
	Matrix(const initializer_list< initializer_list <double> >& iif);
	Matrix(const Matrix& M);

	////////////////  MEMBER FUNCTION BIASA
	Matrix t();

	double& set(int r, int c){
		return entry[r][c];
	}
	double& set(int idx){
		int r = int(idx/col);
		int c = idx % col;
		return entry[r][c];
	}
	size_m size(){
		return std::make_pair(row, col);
	}

	Matrix getDiag();
	void setDiag(const Matrix& rowmat);
	void setDiag(initializer_list< initializer_list<double> > iif){
		Matrix rowmat = iif;
		setDiag(rowmat);
	}

	Matrix block(int ri, int ci, int rf, int cf);
	void setBlock(int ri, int ci, int rf, int cf, const Matrix& mblock);

	void setBlock(int ri, int ci, int rf, int cf, initializer_list< initializer_list<double> > iif){
		Matrix miif = iif;
		setBlock(ri, ci, rf, cf, miif);
	}

	// single row
	Matrix getRow(int r);
	void setRow(int r, const Matrix& mrow);
	void setRow(int r, initializer_list< initializer_list<double> > iif);
	// multiple rows
	Matrix getRow(int ri, int rf);
	void setRow(int ri, int rf, const Matrix& mrow);
	void setRow(int ri, int rf, initializer_list< initializer_list<double> > iif);
	void deleteRow(int ri, int rf);
	void deleteRow(int r);

	// single coloumn
	Matrix getCol(int c);
	void setCol(int c, const Matrix& mcol);
	void setCol(int c, initializer_list < initializer_list<double> > iif);

	// multiple coloumns
	Matrix getCol(int ci, int cf);
	void setCol(int ci, int cf, const Matrix& mcol);
	void setCol(int ci, int cf, initializer_list < initializer_list<double> > iif);
	void deleteCol(int ci, int cf);
	void deleteCol(int c);

	double mul_entry(){
		double val_all(1.);
		for(int r=0; r<row; r++)
			for(int c=0; c<col; c++)
				val_all *= (*this)(r,c);
		return val_all;
	}

	double sum_entry(){
		double sum_all(0.);
		for(int r=0; r<row; r++)
			for(int c=0; c<col; c++)
				sum_all += (*this)(r,c);
		return sum_all;
	}

	Matrix pow_entry(double valpow){
		Matrix respowentry(row, col);
		for(int r=0; r<row; r++)
			for(int c=0; c<col; c++)
				respowentry.set(r,c) = pow(((*this)(r,c)), valpow);
		return respowentry;
	}

	int n_entry(){
		return row * col;
	}

	double det(int method = 0);

	int n_diag(){
		return row<col ? row : col;
	}

	int rank(){
		LUdecomp();
		return n_rank_of_this_matrix;
	}

	std::unordered_map<std::string, Matrix> LUdecomp();
	std::unordered_map<std::string, Matrix> rref();

	std::unordered_map<std::string, Matrix> null_space();

	std::unordered_map<std::string, Matrix> qr(bool simple=false);

	std::unordered_map<std::string, Matrix> qr_iter(double tol=EPSILON, int nloop= -1, bool isprint=false);

	std::unordered_map<std::string, Matrix> eig(bool compute_eigvec=true);
	// std::unordered_map<std::string, Matrix> SchurDecomp();

	////////////////  FRIEND FUNCTION
	// must be friend function,  otherwise not be called
	friend Matrix operator*(double scalar, Matrix M);
	friend Matrix operator/(double scalar, Matrix M);

	////////////////  OPERATOR
	double operator()(int r, int c){
		return entry[r][c];
	}
	double operator()(int idx){
		int r = int(idx/col);
		int c = idx % col;
		return entry[r][c];
	}

	const Matrix& operator=(const Matrix& m);

	bool isZero(){
		for(int r=0; r<row; r++)
			for(int c=0; c<col; c++)
				if(fabs((*this)(r,c)) > EPSILON)
					return false;
		return true;
	}

	bool isUpper(double tol=EPSILON){
		for(int c=0; c<col-1; c++)
			for(int r=c+1; r<row; r++)
				if(fabs((*this)(r,c)) > tol)
					return false;

		if(row > col)
			if(!getRow(col,-1).isZero())
				return false;
		return true;
	}

	bool isLower(double tol=EPSILON){
		for(int c=col-1; c>0; c--)
			for(int r=c-1; r>=0; r--)
				if(fabs((*this)(r,c)) > tol)
					return false;

		return true;
	}

	bool isDiagonal(double tol=EPSILON){
		n_rank_of_this_matrix = rank();
		Matrix rmat = block(0,0,n_rank_of_this_matrix-1,n_rank_of_this_matrix-1);
		if(rmat.isUpper(tol) && rmat.isLower(tol))
			return true;
		return false;
	}

	bool isSame(const Matrix& M, double tol=EPSILON){
		for(int r=0; r<row; r++)
			for(int c=0; c<col; c++)
				if(fabs((*this)(r,c) - M.entry[r][c]) > tol)
					return false;
		return true;
	}

	int numOfNonZero(double tol=EPSILON){
		int nz(0);
		for(int r=0; r<row; r++)
			for(int c=0; c<col; c++)
				if(fabs((*this)(r,c)) > tol) nz++;
		return nz;
	}

	double norm(){
		double val_norm(0.);

		for(int r=0; r<row; r++){
			for(int c=0; c<col; c++){
				val_norm += pow((*this)(r,c),2.);
			}
		}

		val_norm = pow(val_norm, 0.5);

		return val_norm;
	}

	Matrix unit(){
		double val_norm = norm();
		return val_norm>EPSILON ? (*this)/val_norm : Matrix::Zeros(size());
	}

	void setAsUnit(){
		*this = unit();
	}

	Matrix operator-=(const Matrix& M){
		return ((*this) - M);
	}

	Matrix operator+=(const Matrix& M){
		return (*this) + M;
	}

	Matrix operator*(double a);

	Matrix operator*(const Matrix& M);
	Matrix operator+(const Matrix& M);
	Matrix operator-(const Matrix& M);

	Matrix operator-(){
		return (*this) * double(-1.);
	}

	Matrix operator/(double a){
		return (*this) * double(1/a);
	}

	////////////////  STATIC
	static Matrix MatrixNaN(int row_, int col_);
	static Matrix MatrixINF(int row_, int col_);
	static Matrix RandomUniformReal(int row_, int col_, Matrix bound);
	static Matrix Zeros(size_m sm);
	static Matrix Zeros(int row, int col);

	static Matrix Ones(size_m sm);
	static Matrix Ones(int row, int col);

	static Matrix Identity(int ndiag);

	////////////////  OSTREAM
	void print(std::string _name_mat = "", bool print_round=true, int mult_round=1000);
	friend std::ostream& operator<<(std::ostream& os, const Matrix& M);
};

std::ostream& operator<<(std::ostream& os, std::pair<int, int> sm);
std::ostream& operator<<(std::ostream& os, std::unordered_map<std::string, Matrix> um);

Matrix concat(const std::initializer_list< std::initializer_list <Matrix> >& iim);
Matrix createDiag(initializer_list< initializer_list<double> > iif);
Matrix createDiag(const Matrix& mdiag);

std::unordered_map<std::string, Matrix > solve(Matrix A, Matrix b);

double mul_entry(initializer_list< initializer_list <double> > iif);
double dot(Matrix v1, Matrix v2);
Matrix round(Matrix M, int mult=1000);