import numpy as np
from scipy.optimize import fsolve

def compute_T_active(T_irr, Sigma_gas, alpha, OmegaK, kappa=1.0):
    # Constantes
    sigma_sb = 5.67e-5
    k_B = 1.38e-16
    m_p = 1.67e-24
    mu = 2.3
    
    # T^4 = T_irr^4 + T_visc^4
    # T_visc^4 = (3 * tau / 8 / sigma_sb) * (9/4 * nu * Sigma_gas * OmegaK^2)
    # tau = kappa * Sigma_gas / 2
    # nu = alpha * c_s^2 / OmegaK
    # c_s^2 = k_B * T / (mu * m_p)
    #
    # T^4 = T_irr^4 + (3/8/sigma_sb) * (kappa * Sigma_gas / 2) * (9/4) * alpha * (k_B * T / (mu * m_p)) * Sigma_gas * OmegaK
    # T^4 = T_irr^4 + C * T
    
    C = (3.0 / 8.0 / sigma_sb) * (kappa * Sigma_gas / 2.0) * (9.0 / 4.0) * alpha * (k_B / (mu * m_p)) * Sigma_gas * OmegaK
    
    def eq(T):
        return T**4 - T_irr**4 - C * T
    
    T_guess = T_irr + 10.0
    T_sol = fsolve(eq, T_guess)
    return T_sol

print("Testing T_active...")
