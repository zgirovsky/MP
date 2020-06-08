#ifndef MATRIX_PROGONKA_MATRIX_PTR_H
#define MATRIX_PROGONKA_MATRIX_PTR_H

#include <vector>
#include <fstream>
#include "matrix.h"

using namespace std;


// change result
void vector_add_ptr(double *A, double *B, double *result);

//change result
void matrix_add_ptr(double *A, double *B, double *result);

// change A
void negate_vector_ptr(double *A);

//change A
void negate_matrix_ptr(double *A);

// change result
void matrix_vector_mul_ptr(double *A, double *v, double *result, int n);

// change result
void matrix_mul_ptr(double *A, double *B, double *result);

//change matrix
void swap_rows_ptr(double *matrix, int i, int j);

void generate_one_matrix_ptr(double *matrix);

// change matrix, inversal.
void inverse_matrix_ptr(double *matrix, double *inversal);

#endif