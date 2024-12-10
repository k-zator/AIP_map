import numpy as np
import pandas as pd
from AIP_interaction_map.asc import calculate_association_energy as ac
from AIP_interaction_map.asc import calculate_association_dual as ac_dual
from AIP_interaction_map.get_atom_aip_lists import get_atom_aip_pairing_lists
from AIP_interaction_map.get_bipartite_pairing import branching_with_bipartite
from AIP_interaction_map.branch_pairing import branching, add_unique, drop_interactions


def getNonPolar(self, SASA_pass=False, L_Atom=None, L_AIP=None, R_Atom=None, R_AIP=None):
    """Determines the non-hydrogen bond interactions of AIPs accoording to a brnaching
       algorithm of choice (default custom branch, otherwise max bipartite algoritm).
       Subject to the max AIP-AIP distance, with atoms being at most the sum of vdW 
       radii + water diameter away. It is possible to make a hierarchy of interactions
       with ones occuring at less than 0.8 A paired first if bypass_at1 flag is used."""
    
    AipPairs = get_atom_aip_pairing_lists(self, SASA_pass, L_Atom, L_AIP, R_Atom, R_AIP)
    if len(AipPairs) == 0:
        return pd.DataFrame()
    final = []
    if SASA_pass == False:
        self.AipPairs = AipPairs

    # the bypass at 1.0 A: if True, add shorest contacts as certain, and only branch
    # pair match the remainder. To be used as a way of speeding up the calculation
    # when the aip network is complex and expensive.
    if self.bypass_at1 == True:
        final_bypass = add_unique(final, AipPairs[AipPairs[:, 5] <= 0.10])
        for row in final_bypass:
            AipPairs = drop_interactions(row, AipPairs[AipPairs[:, 5] > 0.10])
        final_bypass = pd.DataFrame(final_bypass,
        columns=["L", "R", "Atom_Distance", "L_AIP", "R_AIP", "AIP_Distance", "L_frac", "R_frac"])
        for i, row in final_bypass.iterrows():
            x = row.L_frac * row.R_frac
            if x == 1:
                final_bypass.at[i, "L_frac"] = 1.0
            elif x == 0:
                final_bypass.at[i, "L_frac"] = 0.0
            else:
                final_bypass.at[i, "L_frac"] = 0.5
        final_bypass.rename(columns={"L_frac": "Frac"}, inplace=True)
        final_bypass.drop(["R_frac"], axis=1, inplace=True)

    # PAIRING ALGORITHM KEY PART
    if self.branch == True:
        final_df = branching(AipPairs, final)
    else:
        final_df = branching_with_bipartite(AipPairs)
    if self.bypass_at1 == True:
        final_df = pd.concat([final_df, final_bypass],  ignore_index=True)
    #add also energy and atom type information
    type_L = [self.AIP_L[int(L)].type for L in final_df["L_AIP"]]
    type_R = [self.AIP_R[int(R)].type for R in final_df["R_AIP"]]

    if self.state.dual == False:
        value_L = [self.AIP_L[int(L)].value for L in final_df["L_AIP"]]
        value_R = [self.AIP_R[int(R)].value for R in final_df["R_AIP"]]
        ddG_LR = [ac(L, R, self.solvent) for L, R in zip(value_L, value_R)]
    else:
        value_L, value_R, ddG_LR = ac_dual(self, final_df)

    s35_L = [self.Atom_L[L].sasa35 for L in final_df["L"]]
    s35_R = [self.Atom_R[R].sasa35 for R in final_df["R"]]
    s70_L = [self.Atom_L[L].sasa70 for L in final_df["L"]]
    s70_R = [self.Atom_R[R].sasa70 for R in final_df["R"]]
    ddG_scaled = np.array(ddG_LR) * np.array(final_df.Frac)
    final_df = final_df.assign(L_type=type_L, R_type=type_R, 
        L_value=value_L, R_value=value_R, ddG=ddG_LR, ddG_sc=ddG_scaled,
        s35_L=s35_L, s35_R=s35_R, s70_L=s70_L, s70_R=s70_R)
    return final_df.reset_index(drop=True)
