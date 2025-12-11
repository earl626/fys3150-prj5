#if !defined(SOLVER_H)

#include <armadillo>
#include <complex>

// TODO: Convert this setup and simulation logic into a dedicated Solver class.
//       The class should own parameters like M, m, K, Nt, h, dt, V, A, B, U, u, etc.,
//       and provide member functions for initializing the grid, potential, wavepacket,
//       building CN matrices, running the time evolution, and storing results.
//       This will remove the clutter from main() and ensure consistent handling
//       of simulation state and resources.
//
//       Not implemented now due to deadline constraints and other exams.

/**
 * @brief Convert 2D grid indices to a 1D vector index.
 *
 * This function maps a pair of integer indices (i, j) from a
 * two-dimensional grid into the corresponding single index in a
 * one-dimensional vector using column-major ordering.
 *
 * @param i  Index along the x-direction (0 <= i < m)
 * @param j  Index along the y-direction (0 <= j < m)
 * @param m  Number of interior grid points per dimension.
 * @return   The corresponding 1D index.
 */
int idx(int i, int j, int m);

/**
 * @brief Construct an n-slit potential barrier on a 2D grid.
 *
 * This function creates a matrix V(x,y) representing a vertical barrier
 * centered at x = 0.5 with a given thickness. The barrier contains a
 * user-specified number of open slits, all aligned symmetrically around
 * the domain midpoint. Grid points inside wall regions are assigned the
 * potential height v0, while slit openings remain zero.
 *
 * @param M        Number of grid points along each dimension.
 * @param h        Spatial step size (domain is [0,1] x [0,1]).
 * @param v0       Potential height for barrier regions.
 * @param n_slits  Number of slits to embed in the barrier.
 *
 * @return A matrix V(M x M) containing the constructed potential.
 */
arma::mat init_n_slit_potential(const int M, const double h, const double v0,
                                const int n_slits);

/**
 * @brief Initialize a normalized 2D Gaussian wave packet on the interior grid.
 *
 * This function constructs the complex-valued initial state u(x,y,0) for the
 * Schrodinger simulation. The packet is centered at (x_c, y_c), has Gaussian
 * widths sigma_x and sigma_y, and carries momenta p_x and p_y.
 * Only interior grid points are included, and the resulting vector is normalized.
 *
 * @param M        Total number of grid points including boundaries.
 * @param h        Spatial step size.
 * @param x_c      Initial center of the packet in x.
 * @param y_c      Initial center of the packet in y.
 * @param sigma_x  Gaussian width in the x-direction.
 * @param sigma_y  Gaussian width in the y-direction.
 * @param p_x      Momentum component in the x-direction.
 * @param p_y      Momentum component in the y-direction.
 *
 * @return A normalized complex vector containing the initial wave packet
 *         values on the (M-2)^2 interior points.
 */
arma::cx_vec init_wavepacket(const int M, const double h,
                             const double x_c, const double y_c,
                             const double sigma_x, const double sigma_y,
                             const double p_x, const double p_y);

/**
 * @brief Assemble the Crank–Nicolson matrices A and B for the 2D Schrödinger equation.
 *
 * This function constructs the sparse complex matrices A and B appearing in the
 * Crank–Nicolson update equation
 *
 *      A * u^{n+1} = B * u^{n},
 *
 * using a five-point Laplacian stencil and Dirichlet boundary conditions.
 * Only interior grid points are included, giving matrices of size (M-2)^2 × (M-2)^2.
 * The potential values are sampled from the full M × M grid.
 *
 * @param M   Total number of grid points per spatial dimension (including boundaries).
 * @param h   Spatial step size.
 * @param dt  Time-step size.
 * @param V   Real-valued potential matrix on the full grid.
 * @param A   Output sparse complex matrix for the left-hand Crank–Nicolson operator.
 * @param B   Output sparse complex matrix for the right-hand Crank–Nicolson operator.
 */
void build_CN_AB(const int M, const double h, const double dt,
                 const arma::mat& V,
                 arma::sp_cx_mat& A,
                 arma::sp_cx_mat& B);

/**
 * @brief Store a vectorized wavefunction state into a 3D cube.
 *
 * This function inserts the interior values of a state vector u_vec into
 * the corresponding slice of a 3D Armadillo cube U. Boundary values are
 * not written, as they are assumed to already be zero due to Dirichlet
 * conditions. The mapping from 1D index to 2D grid position follows the
 * same ordering used when assembling the Crank–Nicolson matrices.
 *
 * @param u_vec        Vector containing the (M-2)^2 interior wavefunction values.
 * @param U            Cube storing the wavefunction over time.
 * @param slice_index  Index of the time slice into which the state is stored.
 * @param M            Total number of grid points per spatial dimension.
 */
void store_state_in_cube(const arma::cx_vec& u_vec,
                         arma::cx_cube& U,
                         const int slice_index,
                         const int M);

/**
 * @brief Compute the total probability from a vectorized wavefunction.
 *
 * This function sums |u_k|^2 over all elements of the complex state vector
 * representing the interior grid points of the wavefunction.
 *
 * @param u_vec  Complex vector containing the interior wavefunction values.
 * @return       The total probability (sum of squared magnitudes).
 */
double total_probability(const arma::cx_vec& u_vec);

#define SOLVER_H
#endif
