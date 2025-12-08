
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>

#include "utils.hpp"

using namespace std;
using namespace arma;

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
inline int idx(int i, int j, int m)
{
    return i + m * j;
}

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
mat init_n_slit_potential(const int M, const double h, const double v0,
                          const int n_slits)
{
    // Initialise potential matrix
    mat V = zeros<mat>(M, M);
    
    // x-position of the center of the vertical wall, and its thickness
    const double x_wall     = 0.5;
    const double wall_thick = 0.02;
    
    // Geometric parameters for the slits and the gaps between them (in y)
    const double slit_height = 0.05; // height of each opening
    const double gap_y       = 0.05; // vertical wall segment between neighbouring slits
    const double y_center    = 0.5;  // center of n-slit structure
    
    // Total vertical extent occupied by all slits and interior wall segments
    double total_height = n_slits * slit_height + (n_slits - 1) * gap_y;
    
    // y-coordinate of the bottom edge of the lowest slit
    double y_bottom = y_center - 0.5 * total_height;
    
    // Build a list of [y_min, y_max] intervals for all slit openings
    std::vector<std::pair<double,double>> slits;
    double current = y_bottom;
    
    for (int s = 0; s < n_slits; s++) {
        double y_min = current;               // lower edge of slit s
        double y_max = current + slit_height; // upper edge of slit s
        slits.emplace_back(y_min, y_max);
        
        // Move current position up past the wall segment to the start of next slit
        current = y_max + gap_y;
    }
    
    // Convert wall centre and thickness in x to left/right x-coordinates
    const double x_left  = x_wall - wall_thick / 2.0;
    const double x_right = x_wall + wall_thick / 2.0;
    
    // Loop over all grid points
    for (int j = 0; j < M; j++) {
        double y = j * h;
        
        for (int i = 0; i < M; i++) {
            double x = i * h;
            
            // Only modify points that lie within the vertical wall band in x
            if (x >= x_left && x <= x_right) {
                
                // Check if this (x,y) is inside any of the slit openings
                bool in_slit = false;
                for (const auto& interval : slits) {
                    // If y lies inside one of the [y_min, y_max] intervals,
                    // this point is part of a slit (opening), not a barrier.
                    if (y >= interval.first && y <= interval.second) {
                        in_slit = true;
                        break;
                    }
                }
                
                // If the point is in the wall band and not in a slit opening,
                // assign the barrier potential v0
                if (!in_slit) {
                    V(i, j) = v0;
                }
            }
        }
    }
    
    // Return the completed potential matrix
    return V;
}

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
cx_vec init_wavepacket(const int M, const double h,
                       const double x_c, const double y_c,
                       const double sigma_x, const double sigma_y,
                       const double p_x, const double p_y)
{
    const complex<double> I(0.0, 1.0);
    
    const int m = M - 2;      // interior points per dimension
    const int K = m * m;      // total interior points
    
    // Initial wave function at interior points
    cx_vec u0(K, fill::zeros);
    
    // Accumulator for the squared norm (for normalization)
    double norm2 = 0.0;
    
    // Loop over interior grid
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < m; i++) {
            
            // Map to full-grid coordinates (skip boundaries)
            double x = (i + 1) * h;
            double y = (j + 1) * h;
            
            // Displacement from the packet center
            double dx = x - x_c;
            double dy = y - y_c;
            
            double gaussian =
                exp(- dx*dx / (2*sigma_x*sigma_x)
                    - dy*dy / (2*sigma_y*sigma_y));
            
            complex<double> phase = exp(I * (p_x*x + p_y*y));
            complex<double> value = gaussian * phase;
            
            int k = idx(i, j, m);
            u0(k) = value;
            
            norm2 += norm(value);
        }
    }
    
    // Normalize
    double N = sqrt(norm2);
    if (N > 0)
        u0 /= N;
    
    return u0;
}

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
                 const mat& V,
                 sp_cx_mat& A,
                 sp_cx_mat& B)
{
    const complex<double> I(0.0, 1.0);
    
    const int m = M - 2; // internal points per dimension
    const int K = m * m; // total internal points
    
    A.zeros(K, K);
    B.zeros(K, K);
    
    const complex<double> r = I * dt / (2.0 * h * h);
    
    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < m; i++)
        {
            int k = idx(i, j, m);
            
            // Map internal (i,j) -> full grid (i+1, j+1) because of Dirichlet boundaries
            double vij_real = V(i + 1, j + 1);
            complex<double> Vij(vij_real, 0.0);
            
            complex<double> ak = 1.0 + 4.0 * r + 0.5 * I * dt * Vij;
            complex<double> bk = 1.0 - 4.0 * r - 0.5 * I * dt * Vij;
            
            // Diagonal entries
            A(k, k) = ak;
            B(k, k) = bk;
            
            // Left neighbour (i-1, j)
            if (i > 0)
            {
                int kL = idx(i - 1, j, m);
                A(k, kL) = -r;
                B(k, kL) =  r;
            }
            
            // Right neighbour (i+1, j)
            if (i < m - 1)
            {
                int kR = idx(i + 1, j, m);
                A(k, kR) = -r;
                B(k, kR) =  r;
            }
            
            // Down neighbour (i, j-1)
            if (j > 0)
            {
                int kD = idx(i, j - 1, m);
                A(k, kD) = -r;
                B(k, kD) =  r;
            }
            
            // Up neighbour (i, j+1)
            if (j < m - 1)
            {
                int kU = idx(i, j + 1, m);
                A(k, kU) = -r;
                B(k, kU) =  r;
            }
        }
    }
}

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
inline void
store_state_in_cube(const cx_vec& u_vec,
                    cx_cube& U,
                    const int slice_index,
                    const int M)
{
    const int m = M - 2;
    
    // Boundaries are already zero from U(M, M, Nt, fill::zeros),
    // so we only fill interior points.
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < m; i++) {
            int k = idx(i, j, m);
            U(i + 1, j + 1, slice_index) = u_vec(k);
        }
    }
}

/**
 * @brief Compute the total probability from a vectorized wavefunction.
 *
 * This function sums |u_k|^2 over all elements of the complex state vector
 * representing the interior grid points of the wavefunction.
 *
 * @param u_vec  Complex vector containing the interior wavefunction values.
 * @return       The total probability (sum of squared magnitudes).
 */
double total_probability(const cx_vec& u_vec)
{
    double P = 0.0;
    
    // Sum |u_k|^2 over all interior points
    for (arma::uword k = 0; k < u_vec.n_elem; k++) {
        P += std::norm(u_vec(k));
    }
    return P;
}

int main(int argc, char** argv)
{
    //
    // Read CMD arguments
    //
    
    int min_argument_count = 13;
    if (argc < min_argument_count) {
        print_usage(argv[0]); // print usage if not enough arguments
        return 1;
    }
    
    // Required command line arguments
    const string file_name_prefix = argv[1];         // prefix for output files
    double h                      = atof(argv[2]);   // spatial step size
    const double dt               = atof(argv[3]);   // time step
    const double T_desired        = atof(argv[4]);   // desired total simulation time
    const double x_c              = atof(argv[5]);   // initial packet center in x
    const double sigma_x          = atof(argv[6]);   // wave packet width in x
    const double p_x              = atof(argv[7]);   // momentum in x
    const double y_c              = atof(argv[8]);   // initial packet center in y
    const double sigma_y          = atof(argv[9]);   // wave packet width in y
    const double p_y              = atof(argv[10]);  // momentum in y
    const double v0               = atof(argv[11]);  // barrier height
    const int n_slits             = atoi(argv[12]);  // number of slits in barrier
    
    // Optional flags
    int M_override = -1; // if > 0, overrides M computed from h
    
    // Handle optional arguments
    for (int i = min_argument_count; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--M") {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return 1;
            }
            M_override = atoi(argv[++i]);
        }
        else {
            print_usage(argv[0]);
            return 1;
        }
    }
    
    //
    // Determining the output folder
    //
    
    string folder_path = filesystem::exists("../output/") ? "../output/" : "output/";
    while (!filesystem::exists(folder_path)) {
        cout << "Enter path to output folder: ";
        getline(cin, folder_path);
        
        if (!folder_path.empty())
        {
            if (folder_path.back() != '/') {
                folder_path.push_back('/');
            }
        }
    }
    
    //
    // Set up grid + derived parameters
    //
    
    int M;
    if (M_override > 0) {
        M = M_override;
        h = 1.0 / (M - 1); // adjust h to fit [0,1]
    } else {
        M = static_cast<int>(round(1.0 / h)) + 1;
        h = 1.0 / (M - 1); // adjust h to be consistent
    }
    
    int m  = M - 2;  // interior points per dimension
    int K  = m * m;  // total interior points
    int Nt = static_cast<int>(T_desired / dt) + 1; // number of time steps
    // double T_actual = (Nt - 1) * dt; // actual simulated time
    
    // Potential on full M x M grid (including boundaries)
    mat V = init_n_slit_potential(M, h, v0, n_slits);
    
    // Initial wave function on interior grid as a flattened vector
    cx_vec u = init_wavepacket(M, h, x_c, y_c, sigma_x, sigma_y, p_x, p_y);
    
    // Crank–Nicolson system matrices A (LHS) and B (RHS), sparse complex
    sp_cx_mat A, B;
    build_CN_AB(M, h, dt, V, A, B);
    
    // Storage for full wave function U(x,y,t): M x M spatial grid, Nt time steps
    cx_cube U(M, M, Nt, fill::zeros);
    
    // Storage for time values and total probability at each time step
    vec times(Nt, fill::none);
    vec total_prob(Nt, fill::none);
    
    // Store initial state at time t = 0
    store_state_in_cube(u, U, 0, M);
    times(0) = 0.0;
    total_prob(0)  = total_probability(u);
    
    // Preallocate memory for RHS vector b = B * u^n
    cx_vec b(K, fill::none);
    
    // TODO:
    //   Precompute the LU factorization of A
    //   Upgrade to Armadillo version 12.2 and newer
    //   as it has the functionality we need
    //   (see: spsolve_factoriser)
    
    //
    // Time loop
    //
    
    for (int n = 0; n < Nt - 1; n++)
    {
        b = B * u; // Right-hand side: b = B * u^n
        
        // Solve A * u^{n+1} = b
        bool ok = spsolve(u, A, b, "superlu"); // overwrites u with u^{n+1}
        if (!ok) {
            cerr << "spsolve failed at time step n = " << n << endl;
            return 1;
        }
        
        // Current time index is n+1
        double t_n1 = (n + 1) * dt;
        times(n + 1) = t_n1;
        total_prob(n + 1)  = total_probability(u);
        
        // Store full grid state (with boundaries) at time step n+1
        store_state_in_cube(u, U, n + 1, M);
    }
    
    //
    // Save to files
    //
    
    // Construct output file names based on prefix and folder
    string wavefile = folder_path + file_name_prefix + "_wave.bin";
    string potfile  = folder_path + file_name_prefix + "_potential.bin";
    string probfile = folder_path + file_name_prefix + "_probability.csv";
    
    // Save wave function and potential in raw binary format
    U.save(wavefile, raw_binary);
    V.save(potfile,  raw_binary);
    
    // Save probability deviation: t, (P(t) - 1.0)
    mat out(Nt, 2);
    out.col(0) = times;
    out.col(1) = total_prob - 1.0; // deviation from 1
    out.save(probfile, csv_ascii);
    
    return 0;
}
