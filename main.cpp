#include <iostream>
#include <mpi.h>
#include <fstream>
#include <vector>
#include "matrix.h"
#include "matrix_ptr.h"
#include "matrix_vect.h"

using namespace std;

constexpr int N = 256000; // number of equations
constexpr int K = 256; // count of blocks
constexpr int BLOCK = 4; // size of block-coefficients, should coincide with BLOCK_MATRIX_SIZE from matrix.cpp
constexpr int L = N / K; // size of each block/part
constexpr int BLOCK_SQ = BLOCK * BLOCK;
constexpr int matrix_size = N * BLOCK_SQ;
constexpr int main_eq_size = K * BLOCK_SQ;
constexpr int block_size = L * BLOCK * BLOCK;

/* Solve block-tridiagonal system sequentially ("progonka")
 * -A[i] * X[i - 1] + C[i] * X[i] - B[i] * X[i + 1] = F[i]
 */


vector<vector<double>> right_matrix_tridiagonal(vector<vector<double>> &A, vector<vector<double>> &B,
                                                vector<vector<double>> &C, vector<vector<double>> &F, int size) {
    vector<vector<double>> Y(size, vector<double>(BLOCK, 0));
    vector<vector<double>> U(size - 1);
    vector<vector<double>> G(size);
    vector<double> inverse_C;

    inverse_C = inverse_matrix_vect(C[0]);
    U[0] = matrix_mul_vect(inverse_C, B[0]);
    G[0] = matrix_vector_mul_vect(inverse_C, F[0], BLOCK);

    vector<double> multiplier;
    vector<double> multiplier_F;
    for (int i = 1; i < size - 1; i++) {
        multiplier = matrix_mul_vect(A[i], U[i - 1]);
        multiplier = inverse_matrix_vect(matrix_add_vect(C[i], negate_matrix_vect(multiplier)));
        U[i] = matrix_mul_vect(multiplier, B[i]);

        multiplier_F = matrix_vector_mul_vect(A[i], G[i - 1], BLOCK);
        multiplier_F = vector_add_vect(F[i], negate_vector_vect(multiplier_F));
        G[i] = matrix_vector_mul_vect(multiplier, multiplier_F, BLOCK);
    }
    multiplier = matrix_mul_vect(A[size - 1], U[size - 2]);
    multiplier = matrix_add_vect(C[size - 1], negate_matrix_vect(multiplier));
    multiplier_F = matrix_vector_mul_vect(A[size - 1], G[size - 2], BLOCK);
    multiplier_F = vector_add_vect(F[size - 1], negate_vector_vect(multiplier_F));
    G[size - 1] = matrix_vector_mul_vect(inverse_matrix_vect(multiplier), multiplier_F, BLOCK);

    Y[size - 1] = G[size - 1];
    vector<double> temp;
    for (int i = size - 2; i >= 0; i--) {
        temp = matrix_vector_mul_vect(U[i], Y[i + 1], BLOCK);
        Y[i] = vector_add_vect(G[i], negate_vector_vect(temp));
    }
    return Y;
}

void right_matrix_tridiagonal_ptr(double *A, double *B,
                                  double *C, double *F, double *result, int size) {
    // vector<vector<double>> Y(size); Y is result and shoud be size * BLOCK
    // vector<vector<double>> U(size - 1);
//    vector<vector<double>> G(size);
    auto *inverse_C = new double[BLOCK_SQ];
    auto *C_cp = new double[BLOCK_SQ];
    auto *U = new double[(size - 1) * BLOCK_SQ];
    auto *G = new double[size * BLOCK];

    cp(C, C_cp, BLOCK_SQ);
    inverse_matrix_ptr(&C_cp[0], inverse_C);
    matrix_mul_ptr(inverse_C, B, U);
    matrix_vector_mul_ptr(inverse_C, F, G, BLOCK);

    auto *multiplier = new double[BLOCK_SQ];
    auto *tmp = new double[BLOCK_SQ];
    auto *multiplier_F = new double[BLOCK];
    for (int i = 1; i < size - 1; i++) {
        matrix_mul_ptr(&A[i * BLOCK_SQ], &U[(i - 1) * BLOCK_SQ], multiplier);

        negate_matrix_ptr(multiplier);
        matrix_add_ptr(&C[i * BLOCK_SQ], multiplier, tmp);
        inverse_matrix_ptr(tmp, multiplier);

        matrix_mul_ptr(multiplier, &B[i * BLOCK_SQ], &U[i * BLOCK_SQ]);

        matrix_vector_mul_ptr(&A[i * BLOCK_SQ], &G[(i - 1) * BLOCK], multiplier_F, BLOCK);

        negate_vector_ptr(multiplier_F);
        vector_add_ptr(&F[i * BLOCK], multiplier_F, multiplier_F);

        matrix_vector_mul_ptr(multiplier, multiplier_F, &G[i * BLOCK], BLOCK);
    }

    matrix_mul_ptr(&A[(size - 1) * BLOCK_SQ], &U[(size - 2) * BLOCK_SQ], multiplier);

    negate_matrix_ptr(multiplier);
    matrix_add_ptr(&C[(size - 1) * BLOCK_SQ], multiplier, multiplier);

    matrix_vector_mul_ptr(&A[(size - 1) * BLOCK_SQ], &G[(size - 2) * BLOCK], multiplier_F, BLOCK);

    negate_vector_ptr(multiplier_F);
    vector_add_ptr(&F[(size - 1) * BLOCK], multiplier_F, multiplier_F);

    inverse_matrix_ptr(multiplier, tmp);
    matrix_vector_mul_ptr(tmp, multiplier_F, &G[(size - 1) * BLOCK], BLOCK);

    cp(&G[(size - 1) * BLOCK], &result[(size - 1) * BLOCK], BLOCK);

    for (int i = size - 2; i >= 0; i--) {
        matrix_vector_mul_ptr(&U[i * BLOCK_SQ], &result[(i + 1) * BLOCK], tmp, BLOCK);
        negate_vector_ptr(tmp);
        vector_add_ptr(&G[i * BLOCK], tmp, &result[i * BLOCK]);
    }


    delete[] inverse_C;
    delete[] C_cp;
    delete[] U;
    delete[] G;
    delete[] multiplier;
    delete[] tmp;
    delete[] multiplier_F;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int processes_count;
    MPI_Status status;
    int myrank;

    MPI_Comm_size(MPI_COMM_WORLD, &processes_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    srand(time(nullptr));

    clock_t tp1s = 0, tp1e = 0, tp2s = 0, tp2e = 0;

    // firstly use classic sequential algorithm to get correct results
    vector<vector<double>> X_TRUE;

    // Main matrix in algorithm
    vector<vector<double>> A;
    vector<vector<double>> B;
    vector<vector<double>> C;
    vector<vector<double>> F;

    // Phase 1

    // Allocate local matrixes.
    auto *L_A = new double[matrix_size];
    auto *L_B = new double[matrix_size];
    auto *L_C = new double[matrix_size];
    auto *L_F = new double[N * BLOCK];
    auto *U_A = new double[main_eq_size];
    auto *U_B = new double[main_eq_size];
    auto *U_C = new double[main_eq_size];
    auto *U_F = new double[K * BLOCK];

    // Reduced matrix.
    double *A1;
    double *B1;
    double *C1;
    double *F1;
    A1 = new double[block_size + BLOCK_SQ];
    B1 = new double[block_size + BLOCK_SQ];
    C1 = new double[block_size + BLOCK_SQ];
    F1 = new double[L * BLOCK + BLOCK];

    // Copy of matrix for phase 3.
    auto *part_A = new double[block_size + BLOCK_SQ];
    auto *part_B = new double[block_size + BLOCK_SQ];
    auto *part_C = new double[block_size + BLOCK_SQ];
    auto *part_F = new double[L * BLOCK + BLOCK];

    // Reduce and map matrix.
    if (myrank == 0) {

        // generate matrix
        A.resize(N);
        B.resize(N);
        C.resize(N);
        generate_matrices(A, B, C, N);
        F = generate_vector(N);

        // transform vectors to pointers.
        matrix_to_pointer(A, L_A, N);
        matrix_to_pointer(B, L_B, N);
        matrix_to_pointer(C, L_C, N);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < BLOCK; ++j) {
                L_F[i * BLOCK + j] = F[i][j];
            }
        }

        // make local copy for process 1.
        matrix_to_pointer(A, A1, L);
        matrix_to_pointer(B, B1, L);
        matrix_to_pointer(C, C1, L);
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < BLOCK; ++j) {
                F1[i * BLOCK + j] = F[i][j];
            }
        }

        tp1s = clock();
    }

    MPI_Scatter(&L_A[0], block_size, MPI_DOUBLE, &A1[0], block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(&L_B[0], block_size, MPI_DOUBLE, &B1[0], block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(&L_C[0], block_size, MPI_DOUBLE, &C1[0], block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(&L_F[0], L * BLOCK, MPI_DOUBLE, &F1[0], L * BLOCK, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // make copy for phase 3.
    cp(A1, part_A, block_size);
    cp(B1, part_B, block_size);
    cp(C1, part_C, block_size);
    cp(F1, part_F, L * BLOCK);


    // Phase 1.
    // All matrices work in 0 to L-1 range.
    int k = 0;

    // Get bottom coefficients.
    cp(&A1[BLOCK_SQ], L_A, BLOCK_SQ);
    cp(&B1[BLOCK_SQ], L_B, BLOCK_SQ);
    cp(&C1[BLOCK_SQ], L_C, BLOCK_SQ);
    cp(&F1[BLOCK], L_F, BLOCK);
    auto *temp = new double[BLOCK_SQ];
    auto *temp_block = new double[BLOCK_SQ];
    auto *tmp = new double[BLOCK_SQ];
    for (int r = 2; r <= L - 1; r++) {
        make_block_copy_ptr(&A1[r * BLOCK_SQ], tmp);
        inverse_matrix_ptr(tmp, temp);
        matrix_mul_ptr(&L_C[k], temp, temp_block);

        matrix_mul_ptr(temp_block, &C1[r * BLOCK_SQ], temp);
        negate_matrix_ptr(temp);
        matrix_add_ptr(&L_B[k], temp, &L_C[k]);

        matrix_mul_ptr(temp_block, &B1[r * BLOCK_SQ], &L_B[k]);
        negate_matrix_ptr(&L_B[k]);

        matrix_vector_mul_ptr(temp_block, &F1[r * BLOCK], temp, BLOCK);
        negate_vector_ptr(temp);
        vector_add_ptr(&L_F[k], temp, &L_F[k]);
    }

    // Get top coefficients.
    cp(&A1[BLOCK_SQ * (L - 2)], U_A, BLOCK_SQ);
    cp(&B1[BLOCK_SQ * (L - 2)], U_B, BLOCK_SQ);
    cp(&C1[BLOCK_SQ * (L - 2)], U_C, BLOCK_SQ);
    cp(&F1[BLOCK * (L - 2)], U_F, BLOCK);
    k = 0;
    for (int r = L - 3; r >= 0; r--) {
        make_block_copy_ptr(&B1[r * BLOCK_SQ], tmp);
        inverse_matrix_ptr(tmp, temp);
        matrix_mul_ptr(&U_C[k], temp, temp_block);

        matrix_mul_ptr(temp_block, &A1[r * BLOCK_SQ], &U_A[k]);
        negate_matrix_ptr(&U_A[0]);

        matrix_mul_ptr(temp_block, &C1[r * BLOCK_SQ], temp);
        negate_matrix_ptr(temp);
        matrix_add_ptr(&U_A[k], temp, &U_C[k]);

        matrix_vector_mul_ptr(temp_block, &F1[r * BLOCK], temp, BLOCK);
        negate_vector_ptr(temp);
        vector_add_ptr(&U_F[k], temp, &U_F[k]);
    }

    // Take all matrices together.
    MPI_Gather(&U_A[0], BLOCK_SQ, MPI_DOUBLE, &U_A[0], BLOCK_SQ, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&U_B[0], BLOCK_SQ, MPI_DOUBLE, &U_B[0], BLOCK_SQ, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&U_C[0], BLOCK_SQ, MPI_DOUBLE, &U_C[0], BLOCK_SQ, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&U_F[0], BLOCK, MPI_DOUBLE, &U_F[0], BLOCK, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&L_A[0], BLOCK_SQ, MPI_DOUBLE, &L_A[0], BLOCK_SQ, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&L_B[0], BLOCK_SQ, MPI_DOUBLE, &L_B[0], BLOCK_SQ, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&L_C[0], BLOCK_SQ, MPI_DOUBLE, &L_C[0], BLOCK_SQ, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&L_F[0], BLOCK, MPI_DOUBLE, &L_F[0], BLOCK, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myrank == 0) {
        tp1e = clock();
    }

    // Allocate memory for reduced equations from phase 2.
    auto *Y_k = new double[BLOCK];
    auto *Y_k1 = new double[BLOCK];
    for (int i = 0; i < BLOCK; ++i) {
        Y_k[i] = 0;
        Y_k1[i] = 0;
    }
    auto *Y_res = new double[L * BLOCK];
    double *Y = nullptr;

    MPI_Datatype YS;
    MPI_Type_vector(1, BLOCK, 2 * BLOCK, MPI_DOUBLE, &YS);

    if (myrank == 0) {
        // Phase 2
        auto *reduced_A = new double[2 * K * BLOCK_SQ];
        auto *reduced_B = new double[2 * K * BLOCK_SQ];
        auto *reduced_C = new double[2 * K * BLOCK_SQ];
        auto reduced_F = new double[2 * K * BLOCK];

        for (int i = 0; i < K; i++) {
            for (int j = 0; j < BLOCK_SQ; ++j) {
                reduced_A[2 * i * BLOCK_SQ + j] = U_A[i * BLOCK_SQ + j];
                reduced_A[(2 * i + 1) * BLOCK_SQ + j] = L_A[i * BLOCK_SQ + j];

                reduced_B[2 * i * BLOCK_SQ + j] = U_B[i * BLOCK_SQ + j];
                reduced_B[(2 * i + 1) * BLOCK_SQ + j] = L_B[i * BLOCK_SQ + j];

                reduced_C[2 * i * BLOCK_SQ + j] = U_C[i * BLOCK_SQ + j];
                reduced_C[(2 * i + 1) * BLOCK_SQ + j] = L_C[i * BLOCK_SQ + j];
            }
            for (int j = 0; j < BLOCK; ++j) {
                reduced_F[2 * i * BLOCK + j] = U_F[i * BLOCK + j];
                reduced_F[(2 * i + 1) * BLOCK + j] = L_F[i * BLOCK + j];
            }
        }

        Y = new double[2 * K * BLOCK];
        right_matrix_tridiagonal_ptr(reduced_A, reduced_B, reduced_C, reduced_F, Y, 2 * K);

        // send reduced solutions of Y.
//        for (int k = 1; k < processes_count; k++) {
//            MPI_Send(&Y[k * 2 * BLOCK], BLOCK, MPI_DOUBLE, k, 9, MPI_COMM_WORLD);
//            if (k + 1 < K) {
//                MPI_Send(&Y[(2 * k + 1) * BLOCK], BLOCK, MPI_DOUBLE, k, 10, MPI_COMM_WORLD);
//            }
//        }
//        cp(Y, Y_k, BLOCK);
//        cp(&Y[BLOCK], Y_k1, BLOCK);
    }
//    else {
//        // recieve Ys
//        MPI_Recv(&Y_k[0], BLOCK, MPI_DOUBLE, 0, 9, MPI_COMM_WORLD, &status);
//        if (myrank + 1 < K) {
//            MPI_Recv(&Y_k1[0], BLOCK, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD, &status);
//        }
//    }

    // send reduced solutions of Y.
    MPI_Scatter(&Y[0], 1, YS, &Y_k[0], BLOCK, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(&Y[1], 1, YS, &Y_k1[0], BLOCK, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Phase 3

    if (myrank == 0) {
        tp2s = clock();
    }

    int system_size = L - 1;
    matrix_vector_mul_ptr(part_A, Y_k, tmp, BLOCK);
    negate_matrix_ptr(tmp);
    vector_add_ptr(&part_F[0], tmp, &part_F[0]);

    if (myrank + 1 < processes_count) {
        matrix_vector_mul_ptr(&part_B[(system_size - 1) * BLOCK_SQ], Y_k1, tmp, BLOCK);
        vector_add_ptr(&part_F[(system_size - 1) * BLOCK], tmp, &part_F[(system_size - 1) * BLOCK]);
    }

    right_matrix_tridiagonal_ptr(part_A, part_B, part_C, part_F, Y_res, system_size);

    auto *result_X1 = new double[K * L * BLOCK];

    MPI_Gather(&Y_res[0], BLOCK * L, MPI_DOUBLE, &result_X1[0], BLOCK * L, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myrank == 0) {
        tp2e = clock();
    }

    // Final results
    vector<vector<double>> result_X(N);
    if (myrank == 0) {
//        for (int k = 0; k < K; ++k) {
//            start = k * L;
//            end = start + L - 1;
//            for (int i = start + 1; i < end; i++) {
//                result_X[i] = vector<double>(&result_X1[i * BLOCK], &result_X1[(i + 1) * BLOCK]);
//            }
//            result_X[start] = vector<double>(&Y[2 * k * BLOCK], &result_X1[(2 * k + 1) * BLOCK]);
//            result_X[end] = vector<double>(&Y[(2 * k + 1) * BLOCK], &result_X1[(2 * k + 2) * BLOCK]);
//        }
        X_TRUE = right_matrix_tridiagonal(A, B, C, F, N);

        for (int i = 0; i < N; i++) {
            vector<double> T = matrix_vector_mul_vect(A[i], (i - 1 >= 0 ? X_TRUE[i - 1] : vector<double>(BLOCK, 0)),
                                                      BLOCK);
            T = vector_add_vect(T, matrix_vector_mul_vect(C[i], X_TRUE[i], BLOCK));
            T = vector_add_vect(T,
                                matrix_vector_mul_vect(B[i], (i + 1 < N - 1 ? X_TRUE[i + 1] : vector<double>(BLOCK, 0)),
                                                       BLOCK));
            T = vector_add_vect((T), negate_vector_vect(F[i]));
        }

//        ofstream fout("report.txt", ios_base::trunc);
//        fout << "======================================\n\n";
//        fout << "Diffs between modified algo result and sequential\n";
//        fout << "======================================\n\n";
//        for (int i = 0; i < N; i++) {
//            fout << "***** " << i << " *****\n\n";
//            result_X[i] = vector_add_vect(result_X[i], negate_vector_vect(X_TRUE[i]));
//            print_vector(result_X[i], BLOCK, fout);
//        }

        ofstream f_time_out("time.txt", ios_base::app);
        f_time_out << tp1e - tp1s + tp2e - tp2s << "\n";
    }

    delete[] part_A;
    delete[] part_B;
    delete[] part_C;
    delete[] part_F;
    delete[] L_A;
    delete[] L_B;
    delete[] L_C;
    delete[] L_F;
    delete[] U_A;
    delete[] U_B;
    delete[] U_C;
    delete[] U_F;
    delete[] A1;
    delete[] B1;
    delete[] C1;
    delete[] F1;

    delete[] tmp;
    delete[] temp;
    delete[] temp_block;

    delete[] Y_k;
    delete[] Y_k1;
    delete[] Y_res;

    MPI_Finalize();

    return 0;
}
