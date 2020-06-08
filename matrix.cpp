#include <vector>
#include <fstream>
#include "matrix.h"


using namespace std;


// copy data 'from' pointer into 'to' pointer of 'size'.
void cp(double* from, double* to, int size) {
    for (int i = 0; i < size; ++i) {
        to[i] = from[i];
    }
}

// copy data 'from' pointer into 'to' pointer of size BS*BS.
void make_block_copy_ptr(double* matrix, double* result) {
    for (int i = 0; i < BS*BS; i++) {
        result[i] = matrix[i];
    }
}

// copy matrix in vectors to pointer.
void matrix_to_pointer(const vector<vector<double>>& matrix, double* result, int L) {
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < BS*BS; ++j) {
            result[i*BS*BS+j] = matrix[i][j];
        }
    }
}

int get_rand()
{
    int res = 0;
    while (res == 0)
    {
        res = rand() % 8 - 4;
    }
    return res;
}

// BLOCK
vector<vector<double>> generate_vector(int N)
{
    vector<vector<double>> v(N);
    for (int i = 0; i < N; i++)
    {
        v[i].resize(BS);
        for (int j = 0; j < BS; j++) {
            v[i][j] = rand() % 20 - 10;
        }
    }
    return v;
}
vector<double> generate_zero_matrix()
{
    return vector<double>(BS * BS);
}

vector<double> generate_diagonal_matrix()
{
    vector<double> matrix(BS * BS, 0);
    for (int i = 0; i < BS; i++)
    {
        matrix[i * BS + i] = get_rand();
    }
    return matrix;
}

vector<double> generate_tridiagonal_matrix(vector<double>* A, vector<double>* B)
{
    vector<double> matrix(BS * BS);
    for (int i = 0; i < BS; i++) {
        int d_ind = i * BS + i;
        matrix[d_ind] = 41 + get_rand(); // diagonal dominance is required
        if (d_ind - 1 > 0) {
            matrix[d_ind - 1] = get_rand();
        }
        if (d_ind + 1 < BS * BS) {
            matrix[d_ind + 1] = get_rand();
        }
    }
    return matrix;
}

vector<double> generate_one_matrix()
{
    vector<double> matrix(BS * BS, 0);
    for (int i = 0; i < BS; i++)
    {
        matrix[i * BS + i] = 1;
    }
    return matrix;
}

void generate_matrices(vector<vector<double>>& A, vector<vector<double>>& B, vector<vector<double>>& C, int N)
{
    for (int i = 0; i < N; i++)
    {
        A[i] = generate_diagonal_matrix();
        B[i] = generate_diagonal_matrix();
    }
    A[0] = generate_zero_matrix();
    B[N - 1] = generate_zero_matrix();

    C[0] = generate_tridiagonal_matrix(nullptr, &B[0]);
    C[N - 1] = generate_tridiagonal_matrix(&A[N - 2], nullptr);
    for (int i = 1; i < N - 1; i++)
    {
        C[i] = generate_tridiagonal_matrix(&A[i], &B[i]);
    }
}

void print_matrix(vector<double> A, ostream& out)
{
    out.precision(4);
    for (int i = 0; i < BS; i++)
    {
        for (int j = 0; j < BS; j++)
        {
            out << fixed << A[i * BS + j] << "\t";
        }
        out << "\n";
    }
    out << "\n\n";
}
void print_matrix_ptr(double* A, ostream& out)
{
    out.precision(4);
    for (int i = 0; i < BS; i++)
    {
        for (int j = 0; j < BS; j++)
        {
            out << fixed << A[i * BS + j] << "\t";
        }
        out << "\n";
    }
    out << "\n\n";
}

void print_vector(vector<double> A, int n, ostream& out)
{
    out.precision(8);
    for (int i = 0; i < n; i++)
    {
        out << fixed << A[i] << "\t" << endl;
    }
}

void print_vector_ptr(double* A, int n, ostream& out)
{
    out.precision(8);
    for (int i = 0; i < n; i++)
    {
        out << fixed << A[i] << "\t" << endl;
    }
}
