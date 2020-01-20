from globs import *
import numpy as np
try:
    from scipy.misc import logsumexp
except:
    from  scipy.special import logsumexp
from numpy import logaddexp

def transition_to_matrix(dict):
    transition_matrix = np.zeros((NUM_STATES, NUM_STATES))
    for i in range(NUM_STATES):
        for j in range(NUM_STATES):
            transition_matrix[i,j] = dict[STATES[i]][STATES[j]]
    return np.log(transition_matrix)

def emission_to_matrix(dict):
    emission_matrix = np.zeros((NUM_STATES, NUM_EMISSIONS))
    for i in range(NUM_STATES):
        for j in range(NUM_EMISSIONS):
            emission_matrix[i,j] = dict[STATES[i]][EMISSIONS[j]]
    return np.log(emission_matrix)

def calculate_viterbi(seq, transition_matrix, emission_matrix):
    viterbi_matrix = np.full((NUM_STATES, len(seq)), fill_value=np.NINF, dtype=float)
    trace_matrix = np.zeros((NUM_STATES, len(seq)), dtype=int)
    viterbi_matrix[0, 0] = 0
    for j in range(1, len(seq)):
        for i in range(NUM_STATES):
            vec = viterbi_matrix[:,j-1] + transition_matrix[:,i]
            viterbi_matrix[i, j] = max(vec) + emission_matrix[i, int(seq[j])]
            trace_matrix[i, j] = np.argmax(vec)

    # traceback state sequence
    index = np.argmax(viterbi_matrix[:, -1])
    trace = str(int(trace_matrix[index, -1]))
    for j in range(len(seq)-1):
        index = int(trace[-1])
        trace += str(int(trace_matrix[index, -j-1]))
    trace = trace[::-1]
    print(trace)
    print(trace_matrix)
    
    return viterbi_matrix, trace

def calculate_forward(seq, transition_matrix, emission_matrix):
    forward_matrix = np.full((NUM_STATES, len(seq)), fill_value=np.NINF, dtype=float)
    forward_matrix[0, 0] = 0
    for j in range(1, len(seq)):
        for i in range(NUM_STATES):
            vec = forward_matrix[:, j-1] + transition_matrix[: ,i]
            forward_matrix[i, j] = logsumexp(vec) + emission_matrix[i, int(seq[j])]
    result = logsumexp(forward_matrix[:, -1])
    return forward_matrix, result

def calculate_backward(seq, transition_matrix, emission_matrix):
    backward_matrix = np.full((NUM_STATES, len(seq)), fill_value=np.NINF, dtype=float)
    backward_matrix[:, -1] = 0
    for j in range(len(seq)-1, 0, -1):
        for i in range(NUM_STATES):
            vec = backward_matrix[:, j] + transition_matrix[i, :] + emission_matrix[:, int(seq[j])]
            backward_matrix[i, j-1] = logsumexp(vec)
    result = backward_matrix[0, 0]
    return backward_matrix, result

def calculate_posterior(seq, transition_matrix, emission_matrix):
    forward_matrix, result1 = calculate_forward(seq, transition_matrix, emission_matrix)
    backward_matrix, result2 = calculate_backward(seq, transition_matrix, emission_matrix)
    print(forward_matrix)
    print(backward_matrix)
    return forward_matrix + backward_matrix

def trace_states(seq, posterior_matrix):
    trace = STATES[int(np.argmax(posterior_matrix[:, -1]))]
    for j in range(1, len(seq)):
        trace += STATES[int(np.argmax(posterior_matrix[:, -j-1]))]
    trace = trace[::-1]
    return trace

if __name__ == "__main__":
    np.set_printoptions(precision=3)
    p1 = Protein("AAAAAAA","0111116","sAAAAAe")
    p2 = Protein("AAAAAAA","0222226","sBBBBBe")
    p3 = Protein("AAAAAAA","0333336","sTTTTTe")
    p4 = Protein("AAAAAAA","0444446","sOOOOOe")
    p5 = Protein("AAAAAAA","0111116","sOABTOe")
    p6 = Protein("AAAAAAA","0222226","sOABTOe")
    p7 = Protein("AAAAAAA","0333336","sOABTOe")
    p8 = Protein("AAAAAAA","0444446","sOABTOe")
    p9 = Protein("AAAAAAA","0555556","sOABTOe")
    proteins = [p1, p2, p3, p4, p5, p6, p7, p8, p9]
    emissions = init_emissions(proteins, True)
    transitions = init_transitions(proteins, True)
    emissions_matrix = emission_to_matrix(emissions)
    transitions_matrix = transition_to_matrix(transitions)

    seq = "0112233446"
    #posterior = calculate_posterior(seq, transitions_matrix, emissions_matrix)
    #print(posterior, "\n")

    matrix, trace = calculate_viterbi(seq, transitions_matrix, emissions_matrix)