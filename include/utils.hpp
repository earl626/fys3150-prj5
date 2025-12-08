
#if !defined(UTILS_H)

#include <armadillo>
#include <vector>
#include <string>
#include <iostream>

/**
 * @brief Print usage information for the Schrodinger simulation program.
 *
 * This function prints instructions on how to run the program, including
 * required arguments, optional parameters, and an example command-line call.
 *
 * @param filename Name of the executable (typically argv[0]).
 */
void print_usage(const char* filename);

/**
 * @brief Print the sparsity pattern of a complex sparse matrix.
 *
 * This function visualizes the non-zero structure of a given Armadillo
 * sparse complex matrix by printing a grid where each non-zero entry
 * is marked with a dot. It also reports the matrix dimensions and
 * total number of non-zero elements.
 *
 * @param A The sparse complex matrix whose structure is to be printed.
 */
void print_sp_matrix_structure(const arma::sp_cx_mat& A);

#define UTILS_H
#endif
