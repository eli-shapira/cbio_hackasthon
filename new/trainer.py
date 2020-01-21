from globs import *

def transition_to_matrix(dict):
    transition_matrix = np.zeros((NUM_STATES, NUM_STATES))
    for i in range(NUM_STATES):
        for j in range(NUM_STATES):
            transition_matrix[i, j] = dict[STATES[i]][STATES[j]]
    return np.log(transition_matrix)

def emission_to_matrix(dict):
    emission_matrix = np.zeros((NUM_STATES, NUM_EMISSIONS))
    for i in range(NUM_STATES):
        for j in range(NUM_EMISSIONS):
            emission_matrix[i, j] = dict[STATES[i]][EMISSIONS[j]]
    return np.log(emission_matrix)

def init_emissions(proteins, group=True):
    emissions = np.zeros((NUM_STATES, NUM_EMISSIONS))
    for p in proteins:
        seq = p.group_seq if group else p.aa_seq
        for i in range(p.len):
            emissions[STATE_TO_INDEX[p.structure[i]], int(p.group_seq[i])] += 1

    for j in range(NUM_STATES):
        emission_sum = np.sum(emissions[j,:])
        for k in range(NUM_EMISSIONS):
            if emission_sum == 0:
                emissions[j,k] = 0
                continue
            emissions[j,k] = emissions[j,k] / emission_sum

    return np.log(emissions)

def init_emissions_bi(proteins):
    emissions = np.zeros((NUM_EMISSIONS, NUM_STATES, NUM_EMISSIONS))
    for p in proteins:
        for i in range(1, p.len):
            emissions[int(p.group_seq[i]), STATE_TO_INDEX[p.structure[i]], int(p.group_seq[i-1])] += 1
        emissions[:,0,0] = 1

    for i in range(NUM_EMISSIONS):
        for j in range(NUM_STATES):
            emission_sum = np.sum(emissions[i,j,:])
            for k in range(NUM_EMISSIONS):
                if emission_sum == 0:
                    emissions[i,j,k] = 0
                    continue
                emissions[i,j,k] = emissions[i,j,k] / emission_sum

    return np.log(emissions)

def init_transitions(proteins):
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

    return transition_to_matrix(dict2)
