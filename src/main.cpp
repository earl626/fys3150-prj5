
#include <filesystem>
#include <fstream>
#include <iomanip>

#include "utils.hpp"
#include "solver.hpp"

using namespace std;
using namespace arma;

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
