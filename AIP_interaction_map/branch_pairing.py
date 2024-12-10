import numpy as np
import pandas as pd
from collections import Counter
from itertools import combinations
from collections import defaultdict
from queue import Queue


def find_networks(contacts):
    """Find self-consistent networks within the system of AIP-AIP pairings"""
    adjacency_list = defaultdict(list)
    visited = set()

    # Build adjacency list from the contact list
    for contact in contacts:
        adjacency_list[contact[0]].append(contact[1])
        adjacency_list[contact[1]].append(contact[0])

    # Perform BFS to find connected components
    networks = []
    for node in adjacency_list:
        if node not in visited:
            network = []
            queue = Queue()
            queue.put(node)
            visited.add(node)
            while not queue.empty():
                current_node = queue.get()
                network.append(current_node)

                for neighbor in adjacency_list[current_node]:
                    if neighbor not in visited:
                        queue.put(neighbor)
                        visited.add(neighbor)
            networks.append(network)
    return networks


def branching(AipPairs, final):
    """Wrapper function around the branching algoritm. Create first set of options for 
       the branching and format the final non-polar df."""
    
    #find independent networks in the AipPairs df
    contacts = [(AipPairs[:, 3][i], AipPairs[:, 4][i]) for i in range(0, len(AipPairs))]
    networks = find_networks(contacts)

    network_final = []
    for n in networks:
        nPairs = AipPairs[np.isin(AipPairs[:, 3:5], n).any(axis=1)]
        occ = dict(Counter(nPairs[:, 3]))
        nonunique_L_list = [l for l, n in occ.items() if n > 1]
        scenario = []
        if sum(nonunique_L_list) > len(nonunique_L_list):
            options = nPairs[nPairs[:, 3] == nonunique_L_list[0]]
            #Reach out to the the BRANCHING ALGORITHM
            branch_out(options, final, nPairs, scenario)
        else:
            #If no branching is present, just add unique AIP-AIP contacts
            a = add_unique(final, nPairs)
            scenario.append(np.array(a))

        #select best scenario for the network, i.e. one with most contacts and smallest sum of distances
        fracx = [(s.T[6]*s.T[7]) for s in scenario]
        no_contacts = np.array([sum([1 if x == 1 else 0.5 for x in frac]) for frac in fracx])
        best_no_contacts = [scenario[index] for index in range(
            len(no_contacts)) if no_contacts[index] == max(no_contacts)]
        scaled_dist = [np.array(s).T[5].sum() for s in best_no_contacts]
        network_final.append(best_no_contacts[scaled_dist.index(min(scaled_dist))])

    final_df = pd.DataFrame(np.concatenate(network_final),
        columns=["L", "R", "Atom_Distance", "L_AIP", "R_AIP", "AIP_Distance", "L_frac", "R_frac"])

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


def branch_out(options, final, AipPairs, scenario, for_L=True):
    """The MAIN branching function.
        It uses an input of the first set of contacts made by an AIP: options,
        and creates branches: different pairing contact lists, each with a varied
        first contact. The subsequent contact are dependent on this contact,
        as contacts to the already interacting AIP are not possible. To make this
        more efficient, the algorithm "follows the branch": looks for the alternative 
        contact that is now made impossible and determines a set of possible contacts
        for this atom. Clearly, if more than one is possible, this creates even more
        branches all running in parallel, so the algorithm is pretty expensive if the
        network is more complex, hence the bypass at 0.8 A option.
        To remove redundancy, only returns contacts lists - scenarios, that have at least
        the largest number of contacts to date and is not yet present."""
    if for_L == True:  # direct to correct columns
        L = 3  # "L_AIP"
        R = 4  # "R_AIP"
    else:
        L = 4  # "R_AIP"
        R = 3  # "L_AIP"

    #main loop looking at addition of choice and branching on alternatives
    choices = find_choices(options)
    for choice in choices:
        aippairs = AipPairs.copy()
        sc1 = final.copy()
        if type(choice) is np.float64:  # the no-choice
            aippairs = aippairs[aippairs[:, L] != choice]
            alternatives = options
        else:
            if len(choice) == 2:  # i.e. two interactions are added together for L_frac == 1
                [sc1.append(ch) for ch in choice]
                for ch in choice:
                    aippairs = drop_interactions(ch, aippairs)
                bl = [np.equal(c, options).all(axis=1) for c in choice]
                alternatives = options[~np.logical_or(bl[0], bl[1])]
            else:
                sc1.append(choice)
                aippairs = drop_interactions(choice, aippairs)
                alternatives = options[~np.equal(choice, options).all(axis=1)]
        
        #the recursive function. The detail depends on whether the choice was None 
        if len(alternatives) > 0 and type(choice) is np.float64:
            _, sc2, aippairs = continue_down_branch(
                alternatives, sc1, aippairs, scenario, emptychoice=True)
        
        elif len(alternatives) > 0:
            _, sc2, aippairs = continue_down_branch(
                alternatives, sc1, aippairs, scenario)
        else:
            sc2 = sc1

            #look for another network that could have been created by the choices made,
            #effectively an independent network but made as a consequence
            occ = dict(Counter(aippairs[:, L]))
            nonunique_L_list = [l for l, n in occ.items() if n > 1]
            if len(nonunique_L_list) == 0:
                sc3 = add_unique(sc2, aippairs)
                if (len(scenario) == 0) or (len(sc3) >= max([len(s) for s in scenario])):
                    # if it is not present in the scenario list already
                    sc3 = np.array(sc3)
                    sc3 = sc3[sc3[:, 0].argsort()]
                    if sum([np.array_equal(sc3, s) for s in scenario]) == 0:
                        scenario.append(sc3)
            else:
                options2 = aippairs[aippairs[:, L] == nonunique_L_list[0]]
                branch_out(options2, sc1, aippairs, scenario)


def find_choices(options, for_L=True):
    """Find all possible combinations of contact choices to be added for an AIP.
       Default choice is L in the option by for_L=True.
       Importantly, also adds a 'no choice' that allows to skip an AIP in a branch."""
    if for_L:
        L = 3  # "L_AIP"
        Lfrac = 6  # "L_frac"
        Rfrac = 7  # "R_frac"
    else:
        L = 4  # "R_AIP"
        Lfrac = 7  # "R_frac"
        Rfrac = 6  # "L_frac"

    choices = []
    to_one = options[options[:, Rfrac] == 1]
    [choices.append(row) for row in to_one.copy()]  # full single interaction
    to_half = options[options[:, Rfrac] == 0.5]
    [choices.append(row) for row in to_half.copy()]  # only the halves

    #treat full AIP as two halves so that it could pair with two distinct AIPs
    if sum(options[:, Lfrac] == 1.0) > 1:
        to_one[:, Rfrac] = [0.5 if (i[Rfrac] != 0) else 0 for i in to_one]
        [choices.append(row) for row in to_one]
        # and added as combinations with others
        options[:, Lfrac] = [0.5 if (i[Rfrac] != 0) else 0 for i in options]
        halves_paired = list(combinations(options, 2))[0]
        [choices.append(np.array(pair)) for pair in halves_paired]

    #adding a "no choice" choice
    choices.append(options[:, L][0])
    return choices


def drop_interactions(choice, aippairs, for_L=True):
    """Depending on interactions made in the choice, all the now impossible interactions are deleted."""
    if for_L:  # direct to correct columns
        L = 3  # "L_AIP"
        Lfrac = 6  # "L_frac"
        R = 4  # "R_AIP"
        Rfrac = 7  # "R_frac"
    else:
        L = 4  # "R_AIP"
        Lfrac = 7  # "R_frac"
        R = 3  # "L_AIP"
        Rfrac = 6  # "L_frac"

    aippairs = aippairs[aippairs[:, L] != choice[L]]
    if choice[Rfrac] == 0.5 and len(aippairs[aippairs[:, R] == choice[R]][:, Rfrac]) > 0 \
            and aippairs[aippairs[:, R] == choice[R]][:, Rfrac][0] == 1.0:
        # this one exception when I add 1.0 to 0.5 contact to what was actually a full AIP
        aippairs[aippairs[:, R] == choice[R], Rfrac] = 0.5
    elif choice[Lfrac] == 0.5 and choice[Rfrac] == 1.0:
        aippairs[aippairs[:, R] == choice[R], Rfrac] = 0.5
    else:
        aippairs = aippairs[aippairs[:, R] != choice[R]]
    return aippairs


def continue_down_branch(alternatives, sc1, aippairs, scenario, emptychoice=False, for_L=True):
    """The ITERATIVE function in the branch_out algorithm, also recursive on itself.
       It receives the alternatives to AIP contact newly made and uses them to investigate
       the downstream of the branch - effects of the choice made by adding the new contact.
       One interation is one step, hence it is made recursive to finish each branch.
       It can recursively return to branch_out as well if the branch is finished
       but other network have been created by the recursion."""
    
    if for_L:  # direct to correct columns
        L = 3  # "L_AIP"
        Lfrac = 6  # "L_frac"
        R = 4  # "R_AIP"
        Rfrac = 7  # "R_frac"
    else:
        L = 4  # "R_AIP"
        Lfrac = 7  # "R_frac"
        R = 3  # "L_AIP"
        Rfrac = 6  # "L_frac"

    for alt in alternatives:  #picks R
        if sum(aippairs[:, R] == alt[R]) > 0:  #given contacts to R exist
            R_contacts = aippairs[aippairs[:, R] == alt[R]]
            choices = find_choices(R_contacts, for_L=False)
            for choice in choices:
                aippairs2 = aippairs.copy()
                sc2 = sc1.copy()
                if type(choice) is np.float64 and emptychoice == False:
                    aippairs2 = aippairs2[aippairs2[:, R] != choice]
                    alts = R_contacts
                    emptychoice = True
                elif type(choice) is np.float64 and emptychoice == True:
                    #ignore furter descent down branch due to second non-interacting choice
                    alts = R_contacts
                    continue
                else:

                    if len(choice) == 2:
                        [sc2.append(c) for c in choice]
                        aippairs2 = aippairs2[aippairs2[:, R] != choice[0, R]]
                        aippairs2 = aippairs2[aippairs2[:, R] != choice[1, R]]
                        if sum(np.isin(aippairs2[:, L], choice[:, L])) > 0:
                            alts = aippairs2[np.isin(aippairs2[:, L], choice[:, L])]
                        else:
                            alts = None
                        for c in choice:
                            if c[Lfrac] == 0.0:
                                pass
                            elif c[Lfrac] == 0.5 \
                                and len(aippairs2[aippairs2[:, L] == c[L], Lfrac]) > 0 \
                                    and aippairs2[aippairs2[:, L] == c[L], Lfrac][0] == 1.0:
                                #this one exception when I add 1.0 to 0.5 contact to what was actually a full AIP
                                aippairs2[aippairs2[:, L] == c[L], Lfrac] = 0.5
                                if alts is not None:
                                    aippairs2[np.array(aippairs2[:, L] == c[L]), Lfrac] = 0.5
                                    alts[np.array(alts[:, L] == c[L]), Lfrac] = 0.5
                            elif c[Rfrac] == 0.5 and c[Lfrac] == 1.0:
                                aippairs2[aippairs2[:, L] == c[L], Lfrac] = 0.5
                                if alts is not None:
                                    alts[np.array(alts[:, L] == c[L]), Lfrac] = 0.5
                            else:
                                aippairs2 = aippairs2[aippairs2[:, L] != c[L]]

                    else:
                        sc2.append(choice)
                        aippairs2 = aippairs2[aippairs2[:, R] != choice[R]]
                        if len(aippairs2[aippairs2[:, L] == choice[L]]) > 0:
                            alts = aippairs2[aippairs2[:, L] == choice[L]]
                        else:
                            alts = None
                        if choice[Lfrac] == 0.0:
                            pass
                        elif choice[Lfrac] == 0.5 \
                                and len(aippairs2[aippairs2[:, L] == choice[L], Lfrac]) > 0 \
                                and aippairs2[aippairs2[:, L] == choice[L], Lfrac][0] == 1.0:
                            #this one exception when I add 1.0 to 0.5 contact to what was actually a full AIP
                            aippairs2[aippairs2[:, L] == choice[L], Lfrac] = 0.5
                            if alts is not None:
                                alts[np.array(alts[:, L] == choice[L]), Lfrac] = 0.5
                        elif choice[Rfrac] == 0.5 and choice[Lfrac] == 1.0:
                            aippairs2[np.array(aippairs2[:, L] == choice[L]), Lfrac] = 0.5
                            if alts is not None:
                                alts[np.array(alts[:, L] == choice[L]), Lfrac] = 0.5
                        else:
                            aippairs2 = aippairs2[aippairs2[:, L] != choice[L]]
                if alts is not None:
                    alts, sc2, aippairs2 = continue_down_branch(
                        alts, sc2, aippairs2, scenario, emptychoice)
                    continue

                #Branch complete - look for another network. If present, hand over to branch_out
                occ = dict(Counter(aippairs2[:, L]))
                nonunique_L_list = [l for l, n in occ.items() if n > 1]
                if len(nonunique_L_list) == 0:
                    sc3 = add_unique(sc2, aippairs2)
                    #if scenario is of same or better number of pairing than others
                    if (len(scenario) == 0) or (len(sc3) >= max([len(s) for s in scenario])):
                        sc3 = np.array(sc3)
                        sc3 = sc3[sc3[:, 0].argsort()]
                        if sum([np.array_equal(sc3, s) for s in scenario]) == 0:
                            scenario.append(sc3)
                else:
                    options2 = aippairs2[aippairs2[:, L]
                                         == nonunique_L_list[0]]
                    branch_out(options2, sc2, aippairs2, scenario)
        else:  #no contacts to alternative R -> look for another tree
            #Branch empty - look for another network. If present, hand over to branch_out
            alts = None
            sc2 = sc1.copy()
            aippairs2 = aippairs.copy()
            occ = dict(Counter(aippairs2[:, L]))
            nonunique_L_list = [l for l, n in occ.items() if n > 1]
            if len(nonunique_L_list) == 0:
                sc3 = add_unique(sc2, aippairs2)
                if (len(scenario) == 0) or (len(sc3) >= max([len(s) for s in scenario])):
                    sc3 = np.array(sc3)
                    sc3 = sc3[sc3[:, 0].argsort()]
                    if sum([np.array_equal(sc3, s) for s in scenario]) == 0:
                        scenario.append(sc3)
            else:
                options2 = aippairs2[aippairs2[:, L] == nonunique_L_list[0]]
                branch_out(options2, sc2, aippairs2, scenario)

    return alts, sc2, aippairs2


def add_unique(scen, aippairs_orig, for_L=True):
    """Adds yet unpaired AIP contacts that are present in AipPairs for L (default)."""
    if for_L:
        L = 3
        R = 4
        Lfrac = 6
        Rfrac = 7
    else:
        L = 4
        R = 3
        Lfrac = 7
        Rfrac = 6

    ar = scen.copy()
    aippairs = aippairs_orig.copy()
    #for L not yet added
    for row in aippairs:
        if row[R] not in np.unique([i[R] for i in ar]) and row[L] not in np.unique([i[L] for i in ar]):
            ar.append(row)
        elif (row[Rfrac] == 1.0) and len(np.array(ar)[np.array(ar)[:, R] == row[R], :]) == 2:
            pass
        elif (row[Lfrac] == 1.0) and len(np.array(ar)[np.array(ar)[:, L] == row[L], :]) == 2:
            pass
        elif (row[Lfrac] == 0.5) and (row[Rfrac] == 1.0) and row[L] not in np.unique([i[L] for i in ar]):
            if np.array(ar)[np.array(ar)[:, R] == row[R]][0][Rfrac] == 1.0:
                ar.append(row)
        elif (row[Rfrac] == 1.0) and (row[Lfrac] == 0.5) and len(np.array(ar)) > 0:
            if len(np.array(ar)[np.array(ar)[:, R] == row[R], :]) > 0:
                if (np.array(ar)[np.array(ar)[:, R] == row[R], :][0][Lfrac] == 0.5):
                    ar.append(row)
        elif (row[Lfrac] == 1.0) and (row[Rfrac] == 0.5) and len(np.array(ar)) > 0:
            if len(np.array(ar)[np.array(ar)[:, R] == row[R], :]) > 0:
                if (np.array(ar)[np.array(ar)[:, R] == row[R], :][0][Lfrac] == 0.5):
                    ar.append(row)
        aippairs = aippairs[aippairs[:, L] != row[L]]
        if (row[Rfrac] == 1.0) and (sum(aippairs[:, R] == row[R]) == 1) and (row[Lfrac] == 0.5):
            aippairs[np.array(aippairs[:, R] == row[R]), Rfrac] = 0.5
        else:
            aippairs = aippairs[aippairs[:, R] != row[R]]

    #for R not yet added
    for row in aippairs:
        if row[L] not in np.unique([i[L] for i in ar]) and row[R] not in np.unique([i[R] for i in ar]):
            ar.append(row)
        elif (row[Lfrac] == 1.0) and len(np.array(ar)[np.array(ar)[:, L] == row[L], :]) == 2:
            pass
        elif (row[Rfrac] == 1.0) and len(np.array(ar)[np.array(ar)[:, R] == row[R], :]) == 2:
            pass
        elif (row[Rfrac] == 0.5) and (row[Lfrac] == 1.0) and row[R] not in np.unique([i[L] for i in ar]):
            if np.array(ar)[np.array(ar)[:, L] == row[L]][0][Lfrac] == 1.0:
                ar.append(row)
        elif (row[Lfrac] == 1.0) and (row[Rfrac] == 0.5) and len(np.array(ar)) > 0:
            if len(np.array(ar)[np.array(ar)[:, L] == row[L], :]) > 0:
                if (np.array(ar)[np.array(ar)[:, L] == row[L], :][0][Lfrac] == 0.5):
                    ar.append(row)
        elif (row[Rfrac] == 1.0) and (row[Lfrac] == 0.5) and len(np.array(ar)) > 0:
            if len(np.array(ar)[np.array(ar)[:, R] == row[R], :]) > 0:
                if (np.array(ar)[np.array(ar)[:, R] == row[R], :][0][Lfrac] == 0.5):
                    ar.append(row)
        aippairs = aippairs[aippairs[:, R] != row[R]]
        if (row[Lfrac] == 1.0) and (sum(aippairs[:, L] == row[L]) == 1) and (row[Rfrac] == 0.5):
            aippairs[np.array(aippairs[:, L] == row[L]), Lfrac] = 0.5
        else:
            aippairs = aippairs[aippairs[:, L] != row[L]]

    return ar
