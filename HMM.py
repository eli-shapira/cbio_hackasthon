import numpy as np



#structure
def get_transition_prbs(seqs):

    states_dict = {'O' : 3 , 'H' : 2, 'T':1,'S':0}
    transition_counter ,starting_counter= count_transitions(seqs, states_dict)

    print(transition_counter)
    MLE = calculate_MLE(transition_counter)
    return MLE

# def calculate_start_prb()

def calculate_MLE(transition_counter):


    MLE = np.zeros((4, 4))

    for row, i in zip(transition_counter,range(transition_counter.shape[1])):
        for val, j in zip(row,range(row.shape[0])):
            MLE[i, j] = val/row.sum()

    return MLE




def count_transitions(seqs, states_dict):
    transition_counter = np.zeros((4, 4))
    starting_transition=np.zeros(4)
    for seq in seqs:
        for i in range(len(seq) - 1):
            if i==1:
                starting_transition[states_dict[seq[i]]]+=1
                continue
            cur_state_index = states_dict[seq[i]]
            next_state_index = states_dict[seq[i + 1]]
            transition_counter[cur_state_index, next_state_index] += 1
    return transition_counter,starting_transition

MLE=get_transition_prbs(["STHSTHSHSHSHSOHTTTOSTSHSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS", 'HHHHOTOOTOTOTOHHOTHOT'])


def numpy_to_dict(transition):

    index_to_state={0:'S',1:'T',2:'H',3:'O'}
    transition_prb={}
    for state_col,i in zip(transition,range(transition.shape[0])):
        inner_dict={}
        for val,j in zip(state_col,range(transition.shape[1])):
            inner_dict[index_to_state[j]]=transition[i,j]
        transition_prb[index_to_state[i]]=inner_dict
    print(transition_prb)


def evaluate(true_structure,our_prediciton):

    matches=0
    total=0
    for strc_true,strc_our in zip(true_structure,our_prediciton):
        for i in range(len(strc_true)):
            if strc_our[i]==strc_true[i]:
                matches+=1
            total+=1
    return matches/total

print(MLE)
numpy_to_dict(MLE)


states = ('Healthy', 'Fever')

observations = ('normal', 'cold', 'dizzy')

start_probability = {'Healthy': 0.6, 'Fever': 0.4}

transition_probability = {
    'Healthy': {'Healthy': 0.7, 'Fever': 0.3},
    'Fever': {'Healthy': 0.4, 'Fever': 0.6},
}

emission_probability = {
    'Healthy': {'normal': 0.5, 'cold': 0.4, 'dizzy': 0.1},
    'Fever': {'normal': 0.1, 'cold': 0.3, 'dizzy': 0.6},
}




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