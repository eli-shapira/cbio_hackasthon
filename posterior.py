from globs import *
import numpy as np
import Protein

try:
    from scipy.misc import logsumexp
except:
    from scipy.special import logsumexp
from numpy import logaddexp


def transition_to_matrix(dict):
    transition_matrix = np.zeros((NUM_STATES, NUM_STATES))
    for i in range(NUM_STATES):
        for j in range(NUM_STATES):
            transition_matrix[i, j] = dict[STATES[i]][STATES[j]]
    return np.log(transition_matrix)


def aa_dict():
    counter = 0
    dict = {}
    for a in EMISSIONS:
        dict[a] = counter
        counter += 1
    print(dict)
    return dict


def emission_to_matrix(dict):
    emission_matrix = np.zeros((NUM_STATES, NUM_EMISSIONS))
    for i in range(NUM_STATES):
        for j in range(NUM_EMISSIONS):
            emission_matrix[i, j] = dict[STATES[i]][EMISSIONS[j]]
    return np.log(emission_matrix)


def calculate_viterbi(seq, transition_matrix, emission_matrix):
    viterbi_matrix = np.full((NUM_STATES, len(seq)), fill_value=np.NINF, dtype=float)
    trace_matrix = np.zeros((NUM_STATES, len(seq)), dtype=int)
    viterbi_matrix[0, 0] = 0
    for j in range(1, len(seq)):
        for i in range(NUM_STATES):
            vec = viterbi_matrix[:, j - 1] + transition_matrix[:, i]
            viterbi_matrix[i, j] = max(vec) + emission_matrix[i, int(seq[j])]
            trace_matrix[i, j] = np.argmax(vec)

    # traceback state sequence
    trace = STATES[-1]
    index = np.argmax(viterbi_matrix[:, -1])
    trace += STATES[int(trace_matrix[index, -1])]
    for j in range(len(seq) - 3):
        index = int(STATE_TO_INDEX[trace[-1]])
        trace += STATES[int(trace_matrix[index, -j - 2])]
    trace += STATES[0]
    trace = trace[::-1]

    return viterbi_matrix, trace


def calculate_forward_group(seq, transition_matrix, emission_matrix):
    forward_matrix = np.full((NUM_STATES, len(seq)), fill_value=np.NINF, dtype=float)
    forward_matrix[0, 0] = 0
    for j in range(1, len(seq)):
        for i in range(NUM_STATES):
            vec = forward_matrix[:, j - 1] + transition_matrix[:, i]
            forward_matrix[i, j] = logsumexp(vec) + emission_matrix[i, int(seq[j])]
    result = logsumexp(forward_matrix[:, -1])
    return forward_matrix, result

def calculate_forward(seq, transition_matrix, emission_matrix):
    aa = aa_dict()
    forward_matrix = np.full((NUM_STATES, len(seq)), fill_value=np.NINF, dtype=float)
    forward_matrix[0, 0] = 0
    for j in range(1, len(seq)):
        for i in range(NUM_STATES):
            vec = forward_matrix[:, j - 1] + transition_matrix[:, i]
            forward_matrix[i, j] = logsumexp(vec) + emission_matrix[i, aa[seq[j]]]
    result = logsumexp(forward_matrix[:, -1])
    return forward_matrix, result

def calculate_backward(seq, transition_matrix, emission_matrix):
    aa = aa_dict()
    backward_matrix = np.full((NUM_STATES, len(seq)), fill_value=np.NINF, dtype=float)
    backward_matrix[:, -1] = 0
    for j in range(len(seq) - 1, 0, -1):
        for i in range(NUM_STATES):
            vec = backward_matrix[:, j] + transition_matrix[i, :] + emission_matrix[:, aa[seq[j]]]
            backward_matrix[i, j - 1] = logsumexp(vec)
    result = backward_matrix[0, 0]
    return backward_matrix, result


def calculate_backward_group(seq, transition_matrix, emission_matrix):
    backward_matrix = np.full((NUM_STATES, len(seq)), fill_value=np.NINF, dtype=float)
    backward_matrix[:, -1] = 0
    for j in range(len(seq) - 1, 0, -1):
        for i in range(NUM_STATES):
            vec = backward_matrix[:, j] + transition_matrix[i, :] + emission_matrix[:, int(seq[j])]
            backward_matrix[i, j - 1] = logsumexp(vec)
    result = backward_matrix[0, 0]
    return backward_matrix, result


def calculate_posterior(seq, transition_matrix, emission_matrix):
    forward_matrix, result1 = calculate_forward(seq, transition_matrix, emission_matrix)
    backward_matrix, result2 = calculate_backward(seq, transition_matrix, emission_matrix)
    return forward_matrix + backward_matrix

def calculate_posterior_group(seq, transition_matrix, emission_matrix):
    forward_matrix, result1 = calculate_forward_group(seq, transition_matrix, emission_matrix)
    backward_matrix, result2 = calculate_backward_group(seq, transition_matrix, emission_matrix)
    return forward_matrix + backward_matrix


def trace_states(seq, posterior_matrix):
    trace = STATES[int(np.argmax(posterior_matrix[:, -1]))]
    for j in range(1, len(seq)):
        trace += STATES[int(np.argmax(posterior_matrix[:, -j - 1]))]
    trace = trace[::-1]
    return trace


if __name__ == "__main__":
    np.set_printoptions(precision=3)
    emissions = init_emissions(proteins, True)
    transitions = init_transitions(proteins, True)
    emissions_matrix = emission_to_matrix(emissions)
    transitions_matrix = transition_to_matrix(transitions)

    # seq = "0112233446"
    # posterior = calculate_posterior(seq, transitions_matrix, emissions_matrix)
    # print(posterior, "\n")

