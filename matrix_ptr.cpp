#include <vector>
#include <fstream>
#include "matrix_ptr.h"

using namespace std;

// change result
void vector_add_ptr(double *A, double *B, double *result) {
    for (int i = 0; i < BS; i++) {
        result[i] = A[i] + B[i];
    }
}

// change A
void negate_vector_ptr(double *A) {
    for (int i = 0; i < BS; i++) {
        A[i] = -A[i];
    }
}

// change result
void matrix_vector_mul_ptr(double *A, double *v, double *result, int n) {
    for (int i = 0; i < n; i++) {
        result[i] = 0;
        for (int j = 0; j < n; j++) {
            result[i] += A[i * n + j] * v[j];
        }
    }
}

//change result
void matrix_add_ptr(double *A, double *B, double *result) {
    for (int i = 0; i < BS; i++) {
        for (int j = 0; j < BS; j++) {
            result[i * BS + j] = A[i * BS + j] + B[i * BS + j];
        }
    }
}

//change A
void negate_matrix_ptr(double *A) {
    for (int i = 0; i < BS; i++) {
        for (int j = 0; j < BS; j++) {
            A[i * BS + j] = -A[i * BS + j];
        }
    }
}

// change result
void matrix_mul_ptr(double *A, double *B, double *result) {
    for (int i = 0; i < BS; i++) {
        for (int j = 0; j < BS; j++) {
            result[i * BS + j] = 0;
            for (int k = 0; k < BS; k++) {
                result[i * BS + j] += A[i * BS + k] * B[k * BS + j];
            }
        }
    }
}

//change matrix
void swap_rows_ptr(double *matrix, int i, int j) {
    for (int l = 0; l < BS; l++) {
        swap(matrix[i * BS + l], matrix[j * BS + l]);
    }
}

void generate_one_matrix_ptr(double *matrix) {
    for (int i = 0; i < BS; i++) {
        for (int j = 0; j < BS; ++j) {
            matrix[i * BS + j] = 0;
        }
        matrix[i * BS + i] = 1;
    }
}

// change matrix, inversal.
void inverse_matrix_ptr(double *matrix, double *inversal) {
    generate_one_matrix_ptr(inversal);
    for (int64_t l = BS - 1; l >= 0; l--) {
        if (matrix[l * BS + l] == 0) {
            for (int64_t j = l - 1; j >= 0; j--) {
                if (matrix[j * BS + l] != 0) {
                    swap_rows_ptr(matrix, l, j);
                    swap_rows_ptr(inversal, l, j);
                    break;
                }
            }

            if (matrix[l * BS + l] == 0) {
                //print_matrix(matrix, fout);
                //  cout << "IMPOSSIBLE EVALUATE INVERSE MATRIX\n";
                // system("pause");
                return;
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
}
