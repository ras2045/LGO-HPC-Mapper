# lgo_recalibrate.py - LGO Model Analysis and Recalibration Tool
# License: MIT License (Permissive)
# Role: Allows users and collaborators to analyze the LGO model's behavior 
# and test recalibration strategies against the volatile sequence data.

import math

# --- LGO GEOMETRIC CONSTANTS (Production v1.1 - for analysis consistency) ---
# NOTE: These constants are the finalized values used in the C++ core (lgo_mapper.cpp).
ADDITIVE_FACTOR = 0.01
THRESHOLD = 0.0025  # Pn=400 structural break point
ROOT_SCALING_CONSTANT = 2.730 
LOW_DENSITY_CORRECTION = -0.08
HIGH_DENSITY_CORRECTION = 0.40

# --- LGO CORE CALCULATION ---

def calculate_lgo_magnitude(Pn: int, root_scaling: float = ROOT_SCALING_CONSTANT) -> tuple[int, float, str]:
    """
    Calculates the LGO predicted gap and State Magnitude (Psi_n) 
    using the piecewise geometric formula with production constants.

    Args:
        Pn (int): The prime number being analyzed.
        root_scaling (float): Allows overriding the default root scaling constant 
                              for recalibration testing.

    Returns:
        tuple[int, float, str]: (Predicted Gap, Raw Psi_n Magnitude, Field Classification)
    """
    if Pn < 2:
        return 0, 0.0, "Invalid"
        
    # 1. Base Logarithmic Magnitude (Phi_n)
    Phi_n = math.log(Pn) ** 2
    
    # 2. Volatility Magnitude (Omega_n)
    Omega_n = 1.0 / Pn
    
    # 3. Field Separation Logic
    is_sieve_field = Omega_n > THRESHOLD
    
    if is_sieve_field:
        # Sieve Field (Pn <= 400): High volatility
        # Formula: Phi_n * (1 + C_add) + Correction
        Psi_n_final = Phi_n * (1.0 + ADDITIVE_FACTOR) + LOW_DENSITY_CORRECTION
        field_name = "Sieve Field"
    else:
        # Entropy Field (Pn > 400): Low volatility
        # Formula: C_root * ln(Pn) + Correction
        Psi_n_final = root_scaling * math.log(Pn) + HIGH_DENSITY_CORRECTION
        field_name = "Entropy Field"
    
    predicted_gap = math.floor(Psi_n_final)
    
    return int(predicted_gap), Psi_n_final, field_name

def run_test_cases():
    """Runs a series of tests to show how the LGO analysis module functions."""
    
    # Test primes include Sieve Field (7), Transition Point (401), and Entropy Field (1709, 50021)
    test_primes = [7, 401, 1709, 50021]
    
    print("-" * 70)
    print("LGO Model Analysis Test Cases (Using Production Constants)")
    print("-" * 70)
    print(f"{'Pn':<10} | {'Gap (Pred)':<12} | {'Psi_n Magnitude':<20} | {'Field'}")
    print("-" * 70)
    
    for Pn in test_primes:
        predicted_gap, psi_magnitude, field = calculate_lgo_magnitude(Pn)
        
        marker = " <--- Transition" if Pn == 401 else ""
        
        print(f"{Pn:<10} | {predicted_gap:<12} | {psi_magnitude:.16f:<20} | {field}{marker}")
        
    print("-" * 70)
    print("\nAnalysis test complete. Run 'calculate_lgo_magnitude(Pn, custom_root_scaling)' to test recalibration.")

if __name__ == "__main__":
    run_test_cases()
