from globs import *

def calculate_viterbi(seq, transition_matrix, emission_matrix, bi=False):
    viterbi_matrix = np.full((NUM_STATES, len(seq)), fill_value=np.NINF, dtype=float)
    trace_matrix = np.zeros((NUM_STATES, len(seq)), dtype=int)
    viterbi_matrix[0, 0] = 0
    for j in range(1, len(seq)):
        for i in range(NUM_STATES):
            vec = viterbi_matrix[:, j - 1] + transition_matrix[:, i]
            if bi:
                em = emission_matrix[EMISSIONS.index(seq[j]), i, EMISSIONS.index(seq[j-1])]
            else:
                em = emission_matrix[i, EMISSIONS.index(seq[j])]
            viterbi_matrix[i, j] = max(vec) + em
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

    return trace

def calculate_forward(seq, transition_matrix, emission_matrix, bi):
    forward_matrix = np.full((NUM_STATES, len(seq)), fill_value=np.NINF, dtype=float)
    forward_matrix[0, 0] = 0
    for j in range(1, len(seq)):
        for i in range(NUM_STATES):
            vec = forward_matrix[:, j - 1] + transition_matrix[:, i]
            if bi:
                em = emission_matrix[EMISSIONS.index(seq[j]), i, EMISSIONS.index(seq[j-1])]
            else:
                em = emission_matrix[i, EMISSIONS.index(seq[j])]
            forward_matrix[i, j] = logsumexp(vec) + em
    result = logsumexp(forward_matrix[:, -1])
    return forward_matrix, result

def calculate_backward(seq, transition_matrix, emission_matrix, bi):
    backward_matrix = np.full((NUM_STATES, len(seq)), fill_value=np.NINF, dtype=float)
    backward_matrix[:, -1] = 0
    for j in range(len(seq) - 1, 0, -1):
        for i in range(NUM_STATES):
            if bi:
                em = emission_matrix[EMISSIONS.index(seq[j]), :, EMISSIONS.index(seq[j-1])]
            else:
                em = emission_matrix[:, EMISSIONS.index(seq[j])]
            vec = backward_matrix[:, j] + transition_matrix[i, :] + em
            backward_matrix[i, j - 1] = logsumexp(vec)
    result = backward_matrix[0, 0]
    return backward_matrix, result

def calculate_posterior(seq, transition_matrix, emission_matrix, bi=False):
    forward_matrix, result1 = calculate_forward(seq, transition_matrix, emission_matrix, bi)
    backward_matrix, result2 = calculate_backward(seq, transition_matrix, emission_matrix, bi)
    posterior_matrix = forward_matrix + backward_matrix
    trace = trace_states(seq, posterior_matrix)
    return trace

def trace_states(seq, posterior_matrix):
    trace = STATES[int(np.argmax(posterior_matrix[:, -1]))]
    for j in range(1, len(seq)):
        trace += STATES[int(np.argmax(posterior_matrix[:, -j - 1]))]
    trace = trace[::-1]
    return trace

def batch_viterbi(seqs, transition_matrix, emission_matrix, bi=False):
    results = []
    for seq in seqs:
        results.append(calculate_viterbi(seq, transition_matrix, emission_matrix, bi))
    return results

def batch_posterior(seqs, transition_matrix, emission_matrix, bi=False):
    results = []
    for seq in seqs:
        results.append(calculate_posterior(seq, transition_matrix, emission_matrix, bi))
    return results
