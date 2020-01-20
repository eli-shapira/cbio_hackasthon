import numpy as np
import copy
from globs import *

class HMM:

    def __init__(self, train_p, test_p):
        """
        :param seq: the sequence
        :param initial_emission: n initial emission matrix: initial_emission.tsv - deﬁnes the emission probabilities
        of the states for each nucleotide. Columns correspond to A,C,G,T in order, left to right.
        The rows describe the motif states in order 1,...,k.
        :param p: transition from a background state (B1,B2) to the next state (M1,Bend).
        :param q: transition from Bstart to B1. q represents a prior knowledge about the probability a sequence
        contains a motif. For q = 1, we get an OOPS model.
        """
        self.train_p = train_p
        self.test_p = test_p
        self.e = self.init_emission(train_p)
        self.t = self.init_transition()
        self.s_prob = self.t['s']
        true = []
        pred = []
        for p in self.test_p:
            aa_list = self.viterbi(p.group_seq, ['A', 'B', 'T', 'O'])
            path = ''.join(aa_list)
            true.append(p.structure[1:-1])
            pred.append(path[:-1])

        err = self.evaluate(true, pred)
        print(err)

    def init_emission(self, train_p):
        """
        Init the emission matrix
        :param train_p:
        :return:
        """
        emission = [0 for i in range(7)]
        A_tr, B_tr, T_tr, O_tr = copy.deepcopy(emission), copy.deepcopy(emission), copy.deepcopy(emission), copy.deepcopy(emission)
        states = {'A': A_tr, 'B': B_tr, 'T': T_tr, 'O': O_tr}
        for p in train_p:
            for i in range(1, p.len - 1):
                # states[p[i]] returns the dictionary of the state
                states[p.structure[i]][int(p.group_seq[i])] += 1

        # e = np.zeros((4, 5), dtype=np.float)
        emission = [0 for i in range(7)]
        A_tr, B_tr, T_tr, O_tr = copy.deepcopy(emission), copy.deepcopy(emission), copy.deepcopy(emission), copy.deepcopy(emission)
        # s_tr = np.array([1, 0, 0, 0, 0, 0], dtype=np.float64)
        # e_tr = np.array([0, 0, 0, 0, 0, 1], dtype=np.float64)
        s_tr = [1, 0, 0, 0, 0, 0]
        e_tr = [0, 0, 0, 0, 0, 1]
        emission_prob = {'s': s_tr, 'A': A_tr, 'B': B_tr, 'T': T_tr, 'O': O_tr, 'e': e_tr}
        counter = 0
        for k, s in states.items():
            sum_of_aa = np.sum([s[i] for i in range(7)]).astype(np.float)
            states[k][0], emission_prob[k][NUM_EMISSIONS - 1] = 0, 0
            for i in range(1, NUM_EMISSIONS - 1):
                emission_prob[k][i] = s[i] / sum_of_aa
                # e[counter, aa - 1] = s[str(aa)] / sum_of_aa
            counter += 1

        return emission_prob

    def init_transition(self):
        transition_counter = self.count_transitions()
        MLE = self.calculate_MLE(transition_counter)
        return self.numpy_to_dict(MLE)


    def calculate_MLE(self, transition_counter):
        MLE = np.zeros((NUM_STATES, NUM_STATES))
        for row, i in zip(transition_counter,range(transition_counter.shape[1])):
            for val, j in zip(row,range(row.shape[0])):
                if row.sum() == 0:
                    MLE[i, j] = 0
                    continue
                MLE[i, j] = val/row.sum()

        return MLE

    def count_transitions(self):
        transition_counter = np.zeros((NUM_STATES, NUM_STATES))
        for prot in self.train_p:
            seq = prot.structure
            for i in range(len(seq) - 1):
                cur_state_index = STATE_TO_INDEX[seq[i]]
                next_state_index = STATE_TO_INDEX[seq[i + 1]]
                transition_counter[cur_state_index, next_state_index] += 1
        return transition_counter

    def numpy_to_dict(self, transition):
        transition_prb = {}
        for state_col, i in zip(transition, range(transition.shape[0])):
            inner_dict = {}
            for val,j in zip(state_col, range(transition.shape[1])):
                inner_dict[STATES[j]] = transition[i, j]
            transition_prb[STATES[i]] = inner_dict
        return transition_prb

    def evaluate(self, true_structure, our_prediciton):
        matches = 0
        total = 0
        for strc_true, strc_our in zip(true_structure, our_prediciton):

            for i in range(len(strc_true)):
                curr_prb = 0
                if len(strc_true) != len(strc_our):
                    print(strc_our)
                    print(strc_true)
                    print(len(strc_true))
                    print(len(strc_our))
                    print("oh no")
                    exit(1)
                if strc_our[i] == strc_true[i]:
                    matches += 1
                total += 1

        return matches/total

    def viterbi(self, obs, states):
        path = {s:[] for s in states}  # init path: path[s] represents the path ends with s
        curr_pro = {}
        for s in states:
            curr_pro[s] = self.s_prob[s] * self.e[s][int(obs[1])]
        for i in range(1, len(obs) - 1):
            last_pro = curr_pro
            curr_pro = {}
            for curr_state in states:
                # for last_state in states:
                    # if last_state == 'e':
                    # print("last state " + last_state)
                    # print(obs[i])
                    # print(last_pro[last_state])
                    # print(self.t[last_state][curr_state])
                    # print(last_pro[last_state]*self.t[last_state][curr_state]*self.e[curr_state][int(obs[i])])
                max_pro, last_sta = max(((last_pro[last_state]*self.t[last_state][curr_state]*self.e[curr_state][int(obs[i])], last_state)
                                         for last_state in states))
                # print(max_pro)
                # print(last_sta)
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

