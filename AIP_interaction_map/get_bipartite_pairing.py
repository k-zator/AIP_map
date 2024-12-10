import numpy as np
import pandas as pd
import networkx as nx
from collections import Counter
from AIP_interaction_map.branch_pairing import find_networks



def branching_with_bipartite(AipPairs):

    contacts =  [(AipPairs[:,3][i], AipPairs[:,4][i]) for i in range(0, len(AipPairs))]
    networks = find_networks(contacts)
    network_final = []
    for n in networks:
        nPairs = AipPairs[np.isin(AipPairs[:, 3:8], n).any(axis=1)][:, 3:8]

        for i in nPairs:
            if i[3] == 1.0:
                j = np.array([i[0]+10000, i[1], i[2], 0, 0])
                nPairs = np.vstack((nPairs, j))
            elif i[4] == 1.0:
                j = np.array([i[0], i[1]+10000, i[2], 0, 0])
                nPairs = np.vstack((nPairs, j))

        for_graph = ([(int(i),int(j),{'weight':round(1/k,2)}) for i,j,k in nPairs[:,:3]])

        G = nx.Graph()
        G.add_edges_from(for_graph)
        G_pairs = nx.max_weight_matching(G, weight='weight')

        final = []
        for i in G_pairs:
            if i[0] > 10000:
                pair = np.array([i[0] - 10000, i[1]])
            elif i[1] > 10000:
                pair = np.array([i[0], i[1] - 10000])
            else:
                pair = np.array([i[0], i[1]])
            if pair[0] > max(AipPairs[:,3]):
                pair = np.flip(pair)
            final.append(pair)
        network_final.append(np.array(final))

    network_final = np.concatenate(network_final)
    net_fin_ar = np.array([AipPairs[(AipPairs[:,3] == f[0]) & (AipPairs[:,4] == f[1])][0] for f in network_final])
    final_df = pd.DataFrame(net_fin_ar,
    columns=["L", "R", "Atom_Distance", "L_AIP", "R_AIP", "AIP_Distance", "L_frac", "R_frac"])

    occL = dict(Counter(final_df["L_AIP"]))
    occR = dict(Counter(final_df["R_AIP"]))
    for i, row in final_df.iterrows():
        L = row["L_AIP"]
        if occL[L] == 2:
            final_df.at[i, "L_frac"] = 0.5
        R = row["R_AIP"]
        if occR[R] == 2:
            final_df.at[i, "R_frac"] = 0.5    
            
    for i, row in final_df.iterrows():
        x = row.L_frac * row.R_frac
        if x == 1:
            final_df.at[i, "L_frac"] = 1.0
        elif x == 0:
            final_df.at[i, "L_frac"] = 0.0
        else:
            final_df.at[i, "L_frac"] = 0.5
    final_df.rename(columns={"L_frac": "Frac"}, inplace=True)
    final_df.drop(["R_frac"], axis=1, inplace=True)

    return final_df