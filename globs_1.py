import copy
import numpy as np

BOLD = '\033[01m'
RED = "\033[91m"
CYAN = "\033[96m"
YELLOW = "\033[93m"
PURPLE = "\033[35m"
END = "\033[00m"

# STATES = ['s', 'A', 'B', 'T', 'O', 'e']
# STATE_TO_INDEX = {'s':0, 'A':1, 'B':2, 'T':3, 'O':4, 'e':5}

STATES = ['s', '@', 'a', 'A', '&', 'b', 'B', '!', 't', 'T', 'O', 'e']
STATE_TO_INDEX = {l: STATES.index(l) for l in STATES}


NUM_STATES = len(STATES)


# start, negative, positive, polar, hydrophobic, special, end
# EMISSIONS = ['0', '1', '2', '3', '4', '5', '6']
EMISSIONS = ['0', 'G', 'M', 'L', 'N', 'A', 'S', 'F', 'Y', 'T', 'I', 'W', 'Q', 'X', 'U', 'C', 'P', 'D', 'V', 'E', 'H', 'K', 'R', '1']
NUM_EMISSIONS = len(EMISSIONS)

# class Protein:
#
#     def __init__(self, seq, group_seq, structure):
#         self.len = len(seq)
#         self.seq = seq
#         self.group_seq = group_seq
#         self.structure = structure
#         if len(group_seq) != len(structure):
#             print("len group_seq", len(group_seq), "len structure", len(structure))
#             exit(1)
#         if len(seq) != len(structure):
#             print("len seq", len(seq), "len structure", len(structure))
#             exit(1)

def print_emissions(dict):
    print(BOLD+RED+"\t", end='')
    for em in EMISSIONS:
        print(em+"\t", end='')
    print("\n"+END, end="")
    for st in STATES:
        print(BOLD+RED+st+"\t"+END, end="")
        for em in EMISSIONS:
            print(CYAN+str(dict[st][em])+"\t"+END, end='')
        print("\n", end="")
    print("\n", end="")

def print_transitions(dict):
    print(BOLD+RED+"\t", end='')
    for i in STATES:
        print(i+"\t", end='')
    print("\n"+END, end="")
    for i in STATES:
        print(BOLD+RED+i+"\t", end="")
        for j in STATES:
            print(CYAN+str(dict[i][j])+"\t"+END, end='')
        print("\n", end="")
    print("\n", end="")

def init_emissions(proteins, to_print=False):
    """
    Init the emission matrix
    :param proteins:
    :return:
    """
    dict1, dict2, counter = {}, {}, {}
    for em in EMISSIONS:
        counter[em] = 0
    for st in STATES:
        dict1[st] = copy.deepcopy(counter)

    for p in proteins:
        for i in range(p.len):
            dict1[p.structure[i]][p.aa_seq[i]] += 1

    for st in STATES:
        dict2[st] = copy.deepcopy(counter)
        emission_sum = np.sum([a for a in dict1[st].values()]).astype(np.int)
        for i in range(NUM_EMISSIONS):
            if emission_sum == 0:
                dict2[st][str(i)] = 0
                continue
            print(EMISSIONS[i])
            dict2[st][EMISSIONS[i]] = dict1[st][EMISSIONS[i]] / emission_sum

    if to_print:
        print(PURPLE+"### Emissions ###"+END)
        print_emissions(dict2)

    return dict2

def init_emissions_group(proteins, to_print=False):
    """
    Init the emission matrix
    :param proteins:
    :return:
    """
    dict1, dict2, counter = {}, {}, {}
    for em in EMISSIONS:
        counter[em] = 0
    for st in STATES:
        dict1[st] = copy.deepcopy(counter)

    for p in proteins:
        for i in range(p.len):
            dict1[p.structure[i]][p.group_seq[i]] += 1

    for st in STATES:
        dict2[st] = copy.deepcopy(counter)
        emission_sum = np.sum([a for a in dict1[st].values()]).astype(np.int)
        for i in range(NUM_EMISSIONS):
            if emission_sum == 0:
                dict2[st][str(i)] = 0
                continue
            dict2[st][str(i)] = dict1[st][str(i)] / emission_sum

    if to_print:
        print(PURPLE+"### Emissions ###"+END)
        print_emissions(dict2)

    return dict2

def init_transitions(proteins, to_print=False):
    """
    Init the transition matrix
    :param proteins:
    :return:
    """
    dict1, dict2, counter = {}, {}, {}
    for st in STATES:
        counter[st] = 0
    for st in STATES:
        dict1[st] = copy.deepcopy(counter)

    for p in proteins:
        for i in range(p.len-1):
            dict1[p.structure[i]][p.structure[i+1]] += 1

    for i in STATES:
        dict2[i] = copy.deepcopy(counter)
        transition_sum = np.sum([a for a in dict1[i].values()]).astype(np.int)
        for j in STATES:
            if transition_sum == 0:
                dict2[i][j] = 0
                continue
            dict2[i][j] = dict1[i][j] / transition_sum

    if to_print:
        print(PURPLE+"### Transitions ###"+END)
        print_transitions(dict2)

    return dict2

if __name__ == "__main__":
    p1 = Protein("AAAAAAA","0111116","sAAAAAe")
    p2 = Protein("AAAAAAA","0123456","sBBBBBe")
    proteins = [p1, p2]
    emissions = init_emissions(proteins, True)
    transitions = init_transitions(proteins, True)
