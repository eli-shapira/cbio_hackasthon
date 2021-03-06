import numpy as np
from globs import *

print("3")

#structure
def get_transition_prbs(seqs):

    transition_counter = count_transitions(seqs)
    print(transition_counter)
    MLE = calculate_MLE(transition_counter)

    return MLE


def calculate_MLE(transition_counter):

    MLE = np.zeros((NUM_STATES, NUM_STATES))
    for row, i in zip(transition_counter,range(transition_counter.shape[1])):
        for val, j in zip(row,range(row.shape[0])):
            if row.sum() == 0:
                MLE[i, j] = 0
                continue
            MLE[i, j] = val/row.sum()

    return MLE


def count_transitions(seqs):
    transition_counter = np.zeros((NUM_STATES, NUM_STATES))
    for seq in seqs:
        for i in range(len(seq) - 1):
            cur_state_index = STATE_TO_INDEX[seq[i]]
            next_state_index = STATE_TO_INDEX[seq[i + 1]]
            transition_counter[cur_state_index, next_state_index] += 1
    return transition_counter


MLE=get_transition_prbs(["sAAAAABBBBOTTTOOOOe", 'sBABOe'])


def numpy_to_dict(transition):

    transition_prb = {}
    for state_col, i in zip(transition,range(transition.shape[0])):
        inner_dict = {}
        for val,j in zip(state_col,range(transition.shape[1])):
            inner_dict[STATES[j]] = transition[i,j]
        transition_prb[STATES[i]] = inner_dict
    print(transition_prb)
    return transition_prb


def evaluate(true_structure,our_prediciton):

    matches = 0
    total = 0
    for strc_true,strc_our in zip(true_structure, our_prediciton):
        for i in range(len(strc_true)):
            if strc_our[i] == strc_true[i]:
                matches += 1
            total += 1
    return matches/total

print(MLE)
numpy_to_dict(MLE)


def Viterbit(obs, states, s_pro, t_pro, e_pro):
	path = { s:[] for s in states} # init path: path[s] represents the path ends with s
	curr_pro = {}
	for s in states:
		curr_pro[s] = s_pro[s]*e_pro[s][obs[0]]
	for i in range(1, len(obs)):
		last_pro = curr_pro
		curr_pro = {}
		for curr_state in states:
			max_pro, last_sta = max(((last_pro[last_state]*t_pro[last_state][curr_state]*e_pro[curr_state][obs[i]], last_state)
				                       for last_state in states))
			curr_pro[curr_state] = max_pro
			path[curr_state].append(last_sta)

	# find the final largest probability
	max_pro = -1
	max_path = None
	for s in states:
		path[s].append(s)
		if curr_pro[s] > max_pro:
			max_path = path[s]
			max_pro = curr_pro[s]
		# print '%s: %s'%(curr_pro[s], path[s]) # different path and their probability
	return max_path


if __name__ == '__main__':
	obs = ['normal', 'cold', 'dizzy']
	print (Viterbit(obs, states, start_probability, transition_probability, emission_probability))