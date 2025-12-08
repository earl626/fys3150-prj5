
#include "utils.hpp"

using namespace std;
using namespace arma;

void print_usage(const char* filename) {
    cerr << "Usage: " << filename << "\n"
         << "  <file_name_prefix> <h> <dt> <T> "
         << "<x_c> <sigma_x> <p_x> <y_c> <sigma_y> <p_y> <v0> [options]\n\n"
         << "Arguments:\n"
         << "  file_name_prefix:\tPrefix for output files\n"
         << "  h:\t\t\tSpatial step size (domain is [0,1] x [0,1])\n"
         << "  dt:\t\t\tTime step\n"
         << "  T:\t\t\tDesired total simulation time\n"
         << "  x_c:\t\t\tInitial wave packet center in x\n"
         << "  sigma_x:\t\tWave packet width in x\n"
         << "  p_x:\t\t\tInitial momentum in x\n"
         << "  y_c:\t\t\tInitial wave packet center in y\n"
         << "  sigma_y:\t\tWave packet width in y\n"
         << "  p_y:\t\t\tInitial momentum in y\n"
         << "  v0:\t\t\tBarrier height for double-slit potential\n"
         << "  n_slits:\t\tNumber of slits in the barrier (e.g. 1 = single-slit, 2 = double-slit, 3 = triple-slit)\n\n"
         << "Options:\n"
         << "  --M <int>\t\tOverride grid size M (number of points per dim).\n"
         << "\nExample:\n"
         << "  " << filename
         << "  run1 0.005 2.5e-5 0.002 0.25 0.05 200.0 0.5 0.2 0.0 1e10\n\n"
         << "  This runs the double-slit simulation on a grid with spacing h=0.005,\n"
         << "  time step dt=2.5e-5, total time T=0.002, a Gaussian wave packet\n"
         << "  centered at (0.25, 0.5) with widths sigma_x=0.05, sigma_y=0.2,\n"
         << "  momentum (p_x, p_y) = (200, 0), and barrier height v0 = 1e10.\n\n";
}

void print_sp_matrix_structure(const sp_cx_mat& A)
{
    // Use dynamic 2D vector instead of non-standard VLA
    vector<vector<string>> S(A.n_rows, vector<string>(A.n_cols, " "));
    
    // Iterate through non-zero elements
    int nnz = 0;
    for (auto it = A.begin(); it != A.end(); ++it)
    {
        S[it.row()][it.col()] = "•";
        nnz++;
    }
    
    // Print structure
    cout << "\n";
    for (uword i = 0; i < A.n_rows; i++)
    {
        cout << "| ";
        for (uword j = 0; j < A.n_cols; j++)
        {
            cout << S[i][j] << " ";
        }
        cout << "|\n";
    }
    
    cout << "\nmatrix size: " << A.n_rows << "x" << A.n_cols << "\n";
    cout << "non-zero elements: " << nnz << "\n\n";
}
