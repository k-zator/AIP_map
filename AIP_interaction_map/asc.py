import numpy as np
import numpy.polynomial.polynomial as polynom
from AIP_interaction_map.constants import POLY_COEFFS

def calculate_association_energy(epsilon_i, epsilon_j, solvent, temperature=298.0):
    """Calculate free energy of binding in exchage for AIP-solvent interactions
       for AIPs with epsilon_i and _j values in given solvent."""
    t = POLY_COEFFS[solvent][0]
    poly = POLY_COEFFS[solvent][1]
    EvdW = -5.6  # kJmol^-1.
    GAS_CONSTANT = 8.3144598  # Jmol^-1.
    K_vdW = 0.5 * np.exp(-EvdW * np.power(10, 3) / (GAS_CONSTANT * temperature))

    #including repulsion
    interaction_energy = epsilon_i * epsilon_j # if epsilon_i * epsilon_j < 0.0 else 0.0
    K = 0.5 * np.exp(-(interaction_energy + EvdW) * np.power(10, 3) / (GAS_CONSTANT * temperature))

    association_constant = K + K_vdW
    frac_num = -1.0 + np.sqrt(1.0 + 4 * association_constant * t)
    frac_den = 2 * association_constant * t
    binding_energy = 2 * GAS_CONSTANT * temperature * np.log(frac_num / frac_den) / np.power(10, 3)

    c_frac_num = -1.0 + np.sqrt(1.0 + 8 * t)
    c_frac_den = 4 * t
    confinement_energy = - 2 * GAS_CONSTANT * temperature * np.log(c_frac_num / c_frac_den) / np.power(10, 3)

    solvation_energy_i = polynom.polyval(epsilon_i, poly["negative"]) if epsilon_i < 0 else polynom.polyval(
        epsilon_i, poly["positive"])
    solvation_energy_j = polynom.polyval(epsilon_j, poly["negative"]) if epsilon_j < 0 else polynom.polyval(
        epsilon_j, poly["positive"])

    return binding_energy + confinement_energy - solvation_energy_i - solvation_energy_j

def calculate_solvation_energy(epsilon_i, epsilon_j, solvent, temperature=298.0):
    """Calculate free energy of solvation for AIPs with epsilon_i and _j values in given solvent."""

    poly = POLY_COEFFS[solvent][1]

    solvation_energy_i = polynom.polyval(epsilon_i, poly["negative"]) \
        if epsilon_i < 0 else polynom.polyval(epsilon_i, poly["positive"])
    if epsilon_j is not None:
        solvation_energy_j = polynom.polyval(epsilon_j, poly["negative"]) \
            if epsilon_j < 0 else polynom.polyval(epsilon_j, poly["positive"])
        return solvation_energy_i + solvation_energy_j
    else:
        return solvation_energy_i

def calculate_binding_energy(epsilon_i, epsilon_j, solvent, temperature=298.0):
    """Calculate free energy of binding only for AIPs with epsilon_i and _j values in given solvent."""
    t = POLY_COEFFS[solvent][0]
    EvdW = -5.6  # kJmol^-1.
    GAS_CONSTANT = 8.3144598  # Jmol^-1.
    K_vdW = 0.5 * np.exp(-EvdW * np.power(10, 3) / (GAS_CONSTANT * temperature))

    #including repulsion
    interaction_energy = epsilon_i * epsilon_j #if epsilon_i * epsilon_j < 0.0 else 0.0
    K = 0.5 * np.exp(-(interaction_energy + EvdW) * np.power(10, 3) / (GAS_CONSTANT * temperature))

    association_constant = K + K_vdW
    frac_num = -1.0 + np.sqrt(1.0 + 4 * association_constant * t)
    frac_den = 2 * association_constant * t
    binding_energy = 2 * GAS_CONSTANT * temperature * np.log(frac_num / frac_den) / np.power(10, 3)

    c_frac_num = -1.0 + np.sqrt(1.0 + 8 * t)
    c_frac_den = 4 * t
    confinement_energy = - 2 * GAS_CONSTANT * temperature * np.log(c_frac_num / c_frac_den) / np.power(10, 3)

    return binding_energy + confinement_energy

def calculate_association_dual(self, final_df):
    """Calculate free energy of binding in exchage for AIP-solvent interactions
       for AIPs with epsilon_i and epsilon_j values in given solvent."""

    ddG = []; value_L = []; value_R  = []
    for _, row in final_df.iterrows():
        epsilon_i = self.AIP_L[int(row.L_AIP)].value
        epsilon_j = self.AIP_R[int(row.R_AIP)].value
        dg_solv = - calculate_solvation_energy(epsilon_i, self.solvent) - calculate_solvation_energy(epsilon_j, self.solvent)
        ddg = calculate_binding_energy(epsilon_i, epsilon_j, self.solvent)
        if self.AIP_L[int(row.L_AIP)].dual and self.AIP_R[int(row.R_AIP)].dual:
            epsilon_i_alt = self.AIP_L[int(row.L_AIP)].valueMax
            epsilon_j_alt = self.AIP_L[int(row.L_AIP)].valueMax
            ddg_ijalt = calculate_binding_energy(epsilon_i_alt, epsilon_j_alt, self.solvent)
            if ddg_ijalt < ddg:
                ddg = ddg_ijalt
                epsilon_i = epsilon_i_alt
                epsilon_j = epsilon_j_alt
        elif self.AIP_L[int(row.L_AIP)].dual:
            epsilon_i_alt = self.AIP_L[int(row.L_AIP)].valueMax
            ddg_ialt = calculate_binding_energy(epsilon_i_alt, epsilon_j, self.solvent)
            if ddg_ialt < ddg:
                ddg = ddg_ialt
                epsilon_i = epsilon_i_alt
        elif self.AIP_R[int(row.R_AIP)].dual:
            epsilon_j_alt = self.AIP_L[int(row.L_AIP)].valueMax
            ddg_jalt = calculate_binding_energy(epsilon_i, epsilon_j_alt, self.solvent)
            if ddg_jalt < ddg:
                ddg = ddg_jalt
                epsilon_j = epsilon_j_alt
        
        ddG.append(ddg+dg_solv)
        value_L.append(epsilon_i)
        value_R.append(epsilon_j)

    return value_L, value_R, ddG