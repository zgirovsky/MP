#include <vector>
#include <fstream>
#include "matrix_vect.h"


void swap_rows_vect(vector<double> &matrix, int i, int j) {
    for (int l = 0; l < BS; l++) {
        swap(matrix[i * BS + l], matrix[j * BS + l]);
    }
}

vector<double> negate_vector_vect(vector<double> const &A) {
    vector<double> result(BS);
    for (int i = 0; i < BS; i++) {
        result[i] = -A[i];
    }
    return result;
}


vector<double> negate_matrix_vect(vector<double> const &A) {
    vector<double> result(BS * BS);
    for (int i = 0; i < BS; i++) {
        for (int j = 0; j < BS; j++) {
            result[i * BS + j] = -A[i * BS + j];
        }
    }
    return result;
}

vector<double> matrix_add_vect(vector<double> const &A, vector<double> const &B) {
    vector<double> result(BS * BS);
    for (int i = 0; i < BS; i++) {
        for (int j = 0; j < BS; j++) {
            result[i * BS + j] = A[i * BS + j] + B[i * BS + j];
        }
    }
    return result;
}


vector<double> vector_add_vect(vector<double> const &A, vector<double> const &B) {
    vector<double> result(BS);
    for (int i = 0; i < BS; i++) {
        result[i] = A[i] + B[i];
    }
    return result;
}

vector<double> matrix_vector_mul_vect(vector<double> const &A, vector<double> const &v, int n) {
    vector<double> result(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += A[i * n + j] * v[j];
        }
    }
    return result;
}

vector<double> matrix_mul_vect(vector<double> const &A, vector<double> const &B) {
    vector<double> result(BS * BS);
    for (int i = 0; i < BS; i++) {
        for (int j = 0; j < BS; j++) {
            for (int k = 0; k < BS; k++) {
                result[i * BS + j] += A[i * BS + k] * B[k * BS + j];
            }
        }
    }
    return result;
}

vector<double> inverse_matrix_vect(vector<double> matrix) {
    vector<double> inversal = generate_one_matrix();
    for (int64_t l = BS - 1; l >= 0; l--) {
        if (matrix[l * BS + l] == 0) {
            for (int64_t j = l - 1; j >= 0; j--) {
                if (matrix[j * BS + l] != 0) {
                    swap_rows_vect(matrix, l, j);
                    swap_rows_vect(inversal, l, j);
                    break;
                }
            }

            if (matrix[l * BS + l] == 0) {
                //print_matrix(matrix, fout);
                //  cout << "IMPOSSIBLE EVALUATE INVERSE MATRIX\n";
                //  system("pause");
                return inversal;
            }
        }

        for (int64_t i = l - 1; i >= 0; i--) {
            if (matrix[i * BS + l] != 0) {
                double multiplier = matrix[l * BS + l] / matrix[i * BS + l];
                for (int64_t j = l; j >= 0; j--) {
                    matrix[i * BS + j] *= multiplier;
                    matrix[i * BS + j] -= matrix[l * BS + j];
                }

                for (int64_t j = BS - 1; j >= 0; j--) {
                    inversal[i * BS + j] *= multiplier;
                    inversal[i * BS + j] -= inversal[l * BS + j];
                }
            }
        }
    }

    for (int64_t l = 0; l < BS; l++) {
        for (int64_t j = 0; j < BS; j++) {
            inversal[l * BS + j] /= matrix[l * BS + l];
        }
        matrix[l * BS + l] = 1;

        for (int64_t i = l + 1; i < BS; i++) {
            if (matrix[i * BS + l] != 0) {
                for (int64_t j = 0; j < BS; j++) {
                    inversal[i * BS + j] /= matrix[i * BS + l];
                    inversal[i * BS + j] -= inversal[l * BS + j];
                }

                for (int64_t j = l + 1; j < i + 1; j++) {
                    matrix[i * BS + j] /= matrix[i * BS + l];
                }

                matrix[i * BS + l] = 0;
            }
        }
    }
    return inversal;
}