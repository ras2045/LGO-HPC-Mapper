// lgo_mapper.cpp - Law of Geometric Order (LGO) Continuous Prime Gap Mapper
//
// Role: High-performance C++ core. This file contains the proprietary geometric formula
// and is the engine for high-throughput data generation.
//
// **LICENSING NOTICE:** This file is covered by the **PROPRIETARY USE LICENSE** found in 
// PROPRIETARY_LICENSE.txt. Commercial use requires a separate written agreement.
//
// COMPILATION (GCC/Clang): g++ lgo_mapper.cpp -o lgo_mapper -std=c++17 -O3 -fopenmp
// Run Example: ./lgo_mapper 1000000 (Maps 1 million primes)

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>
#include <limits>
#include <iomanip>
#include <omp.h>
#include <fstream>
#include <csignal>
#include <algorithm>

using PrimeType = long long;

// --- GLOBAL STATE & FILES ---
PrimeType global_last_prime_found = 0; // Tracks the last successfully processed Pn
const std::string CHECKPOINT_FILE = "lgo_checkpoint.txt";
const std::string LGO_LOG_FILE = "lgo_run_log.csv";
const std::string VOLATILE_FILE = "lgo_volatile_sequence.txt"; // Channel for high-volatility primes

// --- LGO GEOMETRIC CONSTANTS (v1.1 Certified Defaults - Trade Secrets) ---
// Note: ROOT_SCALING_CONSTANT can be overridden by the LGO_CUSTOM_ROOT_SCALING environment variable.
constexpr double ADDITIVE_FACTOR = 0.01;
constexpr double THRESHOLD = 0.0025;
constexpr double ROOT_SCALING_CONSTANT_DEFAULT = 2.730; 
constexpr double LOW_DENSITY_CORRECTION = -0.08;
constexpr double HIGH_DENSITY_CORRECTION = 0.40;
constexpr PrimeType DEFAULT_START_PRIME = 401;

// --- DIGITAL WATERMARK (DO NOT REMOVE) ---
// This unique string is embedded in the binary to prove authorship and copyright ownership.
// It is intentionally complex and unused to survive code compilation/optimization.
const char* INTERNAL_ASSET_TAG = "LGO-HPC-V1.1-PROPRIETARY-CORE-ASSET-ID-20251203X00A";

// --- Data Structure for LGO Output ---
struct LgoResult {
    int predicted_gap;
    double psi_magnitude;
};

// --- CONFIGURATION FUNCTIONS ---
double get_active_root_scaling() {
    double active_constant = ROOT_SCALING_CONSTANT_DEFAULT;
    const char* custom_val = std::getenv("LGO_CUSTOM_ROOT_SCALING");
    
    if (custom_val) {
        try {
            active_constant = std::stod(custom_val);
            std::cout << "[TUNING ON] Using custom ROOT_SCALING_CONSTANT: " << std::fixed << std::setprecision(8) << active_constant << "\n";
            return active_constant;
        } catch (const std::exception& e) {
            std::cerr << "Warning: LGO_CUSTOM_ROOT_SCALING invalid. Using default " << ROOT_SCALING_CONSTANT_DEFAULT << ".\n";
        }
    }
    std::cout << "[TUNING OFF] Using default certified ROOT_SCALING_CONSTANT: " << std::fixed << std::setprecision(8) << ROOT_SCALING_CONSTANT_DEFAULT << "\n";
    return active_constant;
}

PrimeType get_end_prime_limit() {
    const char* limit_val = std::getenv("LGO_END_PRIME");
    if (limit_val) {
        try {
            return std::stoll(limit_val);
        } catch (const std::exception& e) {
            std::cerr << "Warning: LGO_END_PRIME invalid. Ignoring limit.\n";
        }
    }
    return std::numeric_limits<PrimeType>::max(); // No effective limit
}

// --- LGO CORE CALCULATION ---
LgoResult calculate_lgo_gap(PrimeType Pn, double active_root_scaling) {
    if (Pn < 2) return {0, 0.0};
        
    double log_Pn = std::log(static_cast<double>(Pn));
    double Phi_n = log_Pn * log_Pn;
    double Omega_n = 1.0 / static_cast<double>(Pn);
    bool is_sieve_field = Omega_n > THRESHOLD;
    
    double Psi_n_final = 0.0;
    
    // The core piecewise formula: Sieve Field vs. Entropy Field
    if (is_sieve_field) {
        // Sieve Field: Log-squared magnitude + Additive factor + Low Density Correction
        Psi_n_final = Phi_n * (1.0 + ADDITIVE_FACTOR) + LOW_DENSITY_CORRECTION;
    } else {
        // Entropy Field: Root-scaled magnitude + High Density Correction
        Psi_n_final = active_root_scaling * log_Pn + HIGH_DENSITY_CORRECTION;
    }
    
    int predicted_gap = static_cast<int>(std::floor(Psi_n_final));
    
    return {predicted_gap, Psi_n_final};
}

// --- Primality Test (Optimized Trial Division) ---
bool is_prime(PrimeType n) {
    if (n <= 1) return false;
    if (n <= 3) return true; 
    if (n % 2 == 0 || n % 3 == 0) return false;
    
    for (PrimeType i = 5; i * i <= n; i = i + 6) {
        if (i * i > n) break; // Optimization break
        if (n % i == 0 || n % (i + 2) == 0) {
            return false;
        }
    }
    return true;
}

// --- CHECKPOINTING AND RESUMPTION LOGIC (SIGTERM Handling) ---
void save_checkpoint(PrimeType Pn) {
    std::ofstream outfile(CHECKPOINT_FILE);
    if (outfile.is_open()) {
        outfile << Pn;
        outfile.close();
    } else {
        std::cerr << "\n[CRITICAL ERROR] Failed to write checkpoint file.\n";
    }
}

void sigterm_handler(int signum) {
    std::cerr << "\n--- CAUGHT SIGNAL " << signum << " (Graceful Exit) ---\n";
    if (global_last_prime_found > 0) {
        save_checkpoint(global_last_prime_found);
    }
    std::exit(signum); 
}

PrimeType get_start_prime() {
    PrimeType start_prime = DEFAULT_START_PRIME;
    const char* env_val = std::getenv("LGO_START_PRIME");
    
    if (env_val) {
        try {
            start_prime = std::stoll(env_val);
            std::cout << "RESUMING from ENVIRONMENT variable LGO_START_PRIME: " << start_prime << "\n";
            return start_prime;
        } catch (const std::exception& e) {
            std::cerr << "Warning: LGO_START_PRIME invalid. Checking checkpoint file.\n";
        }
    }

    std::ifstream infile(CHECKPOINT_FILE);
    if (infile.is_open()) {
        if (infile >> start_prime) {
            std::cout << "RESUMING from CHECKPOINT file: " << start_prime << "\n";
            infile.close();
            std::remove(CHECKPOINT_FILE.c_str()); // Remove checkpoint after loading
            return start_prime;
        }
        infile.close();
    }

    std::cout << "Starting from default Pn: " << DEFAULT_START_PRIME << ".\n";
    return DEFAULT_START_PRIME;
}

// --- Main Application Logic ---
int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    // Register signal handlers for graceful exit (Ctrl+C or kill signal)
    std::signal(SIGINT, sigterm_handler);
    std::signal(SIGTERM, sigterm_handler);

    const double active_root_scaling = get_active_root_scaling();
    const PrimeType end_prime_limit = get_end_prime_limit();
    const bool parallelism_enabled = std::getenv("LGO_DISABLE_PARALLELISM") == nullptr;

    // --- Data Logging Setup ---
    std::ofstream logfile(LGO_LOG_FILE, std::ios_base::app); 
    std::ofstream volatile_file(VOLATILE_FILE, std::ios_base::app);
    if (!logfile.is_open() || !volatile_file.is_open()) {
        std::cerr << "[CRITICAL ERROR] Failed to open one or more log files. Exiting.\n";
        return 1;
    }

    // Write CSV header only if the log file is new/empty
    if (logfile.tellp() == 0) {
        logfile << "Pn,LgoMagnitude,PredictedGap,ActualGap,NextPrime,ChecksPerformed,RootScalingConstantUsed\n";
    }

    PrimeType current_prime = get_start_prime();
    int max_iterations = (argc > 1) ? std::stoi(argv[1]) : 50; 
    int count = 0;
    
    std::cout << "-----------------------------------------------------------------------------------------------------------------\n";
    std::cout << "| " << std::left << std::setw(18) << "Pn"
              << " | " << std::left << std::setw(18) << "LGO Mag (Psi_n)"
              << " | " << std::left << std::setw(15) << "Predicted G_n"
              << " | " << std::left << std::setw(10) << "Actual Gap" 
              << " | " << std::left << std::setw(20) << "Next Prime (Pn+1)" 
              << " | " << "Checks" << " |\n";
    std::cout << "-----------------------------------------------------------------------------------------------------------------\n";
    
    while (count < max_iterations && current_prime < end_prime_limit) {
        PrimeType Pn_start = current_prime; 
        
        LgoResult result = calculate_lgo_gap(Pn_start, active_root_scaling);
        PrimeType lgo_gap_pred = result.predicted_gap;
        double psi_magnitude = result.psi_magnitude;
        
        PrimeType start_search = Pn_start + 1;
        PrimeType end_search = Pn_start + lgo_gap_pred;
        
        PrimeType true_next_prime = 0;
        PrimeType checks_performed = 0;
        
        // --- HPC Parallel Search Block ---
        if (parallelism_enabled) {
            #pragma omp parallel for shared(true_next_prime) 
            for (PrimeType search_num_i = start_search; search_num_i <= end_search; ++search_num_i) {
                if (true_next_prime != 0) continue; 
                if (is_prime(search_num_i)) {
                    #pragma omp critical
                    {
                        if (true_next_prime == 0) { 
                            true_next_prime = search_num_i;
                        }
                    }
                    #pragma omp flush(true_next_prime)
                }
            }
        } else {
            // --- Serial Search Block ---
            for (PrimeType search_num_i = start_search; search_num_i <= end_search; ++search_num_i) {
                 if (is_prime(search_num_i)) {
                    true_next_prime = search_num_i;
                    break;
                }
            }
        }
        
        // --- Volatility Search (If Prediction Failed) ---
        if (true_next_prime == 0) {
             // LGO prediction failed (G_actual > G_pred). Search continues indefinitely.
             PrimeType search_beyond = end_search + 1;
             while (true_next_prime == 0) {
                 if (is_prime(search_beyond)) {
                    true_next_prime = search_beyond;
                    break;
                 }
                 search_beyond++;
             }
        }
        
        PrimeType actual_gap = true_next_prime - Pn_start;
        checks_performed = actual_gap / 3; // Approximation of checks

        // --- VOLATILE SEQUENCE CHANNEL LOGIC ---
        // Log Pn if the LGO prediction was exceeded (Actual > Predicted).
        if (actual_gap > lgo_gap_pred) {
            volatile_file << Pn_start << "," << actual_gap << "," << lgo_gap_pred << "\n";
        }

        global_last_prime_found = true_next_prime;

        // --- Console and Log Output ---
        std::cout << "| " << std::left << std::setw(18) << Pn_start
                  << " | " << std::left << std::setw(18) << std::fixed << std::setprecision(8) << psi_magnitude
                  << " | " << std::left << std::setw(15) << lgo_gap_pred
                  << " | " << std::left << std::setw(10) << actual_gap 
                  << " | " << std::left << std::setw(20) << true_next_prime 
                  << " | " << checks_performed << " |\n";

        logfile << Pn_start << ","
                << std::fixed << std::setprecision(10) << psi_magnitude << "," 
                << lgo_gap_pred << ","
                << actual_gap << ","
                << true_next_prime << ","
                << checks_performed << ","
                << active_root_scaling << "\n";
        
        current_prime = true_next_prime;
        count++;
    }
    
    save_checkpoint(current_prime);
    
    std::cout << "-----------------------------------------------------------------------------------------------------------------\n";
    if (current_prime >= end_prime_limit) {
        std::cout << "SESSION STOPPED. Reached end prime limit: " << end_prime_limit << "\n";
    } else {
        std::cout << "SESSION COMPLETE. Final Pn: " << current_prime << " (Checkpoint saved).\n";
    }
    std::cout << "-----------------------------------------------------------------------------------------------------------------\n";
    
    logfile.close();
    volatile_file.close();

    return 0;
}
