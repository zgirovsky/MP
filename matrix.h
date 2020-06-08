#ifndef MATRIX_PROGONKA_MATRIX_H
#define MATRIX_PROGONKA_MATRIX_H

#include <vector>
#include <fstream>
using namespace std;

// should coincide with BLOCK from main.cpp
constexpr int BS = 4;

// All generates, prints, transform from vector to pointer functions.

// copy data 'from' pointer into 'to' pointer of 'size'.
void cp(double* from, double* to, int size);
// copy data 'from' pointer into 'to' pointer of size BS*BS.
void make_block_copy_ptr(double* matrix, double* result);
// copy matrix in vectors to pointer.
void matrix_to_pointer(const vector<vector<double>>& matrix, double* result, int L);

int get_rand();

vector<vector<double>> generate_vector(int N);
vector<double> generate_zero_matrix();
vector<double> generate_diagonal_matrix();
vector<double> generate_tridiagonal_matrix(vector<double>* A, vector<double>* B);
vector<double> generate_one_matrix();
void generate_matrices(vector<vector<double>>& A, vector<vector<double>>& B, vector<vector<double>>& C, int N);

void print_matrix(vector<double> A, ostream& out);
void print_matrix_ptr(double* A, ostream& out);
void print_vector(vector<double> A, int n, ostream& out);
void print_vector_ptr(double* A, int n, ostream& out);

#endif