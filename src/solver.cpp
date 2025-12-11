
#include "solver.hpp"

using namespace std;
using namespace arma;

int idx(int i, int j, int m)
{
    return i + m * j;
}

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

void store_state_in_cube(const cx_vec& u_vec,
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

double total_probability(const cx_vec& u_vec)
{
    double P = 0.0;
    
    // Sum |u_k|^2 over all interior points
    for (arma::uword k = 0; k < u_vec.n_elem; k++) {
        P += std::norm(u_vec(k));
    }
    return P;
}
