//
// Created by hanna on 6.01.20.
//

#ifndef LIHODED_MATRIX_VECT_H
#define LIHODED_MATRIX_VECT_H

#include <vector>
#include <fstream>
#include "matrix.h"


void swap_rows_vect(vector<double> &matrix, int i, int j);

vector<double> negate_vector_vect(vector<double> const &A);

vector<double> negate_matrix_vect(vector<double> const &A);

vector<double> matrix_add_vect(vector<double> const &A, vector<double> const &B);

vector<double> vector_add_vect(vector<double> const &A, vector<double> const &B);

vector<double> matrix_vector_mul_vect(vector<double> const &A, vector<double> const &v, int n);

vector<double> matrix_mul_vect(vector<double> const &A, vector<double> const &B);

vector<double> inverse_matrix_vect(vector<double> matrix);


#endif //LIHODED_MATRIX_VECT_H
