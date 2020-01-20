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
        self.e = self.init_emission(train_p)
        self.t = self.init_transmission(train_p)
        self.viterbi(test_p)

    class Protein:

        def __init__(self, seq, states):
            self.len = len(seq)
            self.group_seq = seq
            self.structure = states
            if len(seq) != len(states):
                print("len seq")
                print(len(seq))
                print(len(states))
                print(len(states))
                exit(1)

    def init_emission(self, train_p):
        """
        Init the emission matrix
        :param train_p:
        :return:
        """
        A_tr, B_tr, T_tr, O_tr = copy.deepcopy(EMISSIONS), copy.deepcopy(EMISSIONS), copy.deepcopy(EMISSIONS), copy.deepcopy(EMISSIONS)
        states = {'A': A_tr, 'B': B_tr, 'T': T_tr, 'O': O_tr}
        for p in train_p:
            for i in range(p.len):
                # states[p[i]] returns the dictionary of the state
                states[p.structure[i]][p.group_seq[i]] += 1

        # e = np.zeros((4, 5), dtype=np.float)
        A_tr, B_tr, T_tr, O_tr = copy.deepcopy(EMISSIONS), copy.deepcopy(EMISSIONS), copy.deepcopy(EMISSIONS), copy.deepcopy(EMISSIONS)
        emission_prob = {'A': A_tr, 'B': B_tr, 'T': T_tr, 'O': O_tr}
        counter = 0
        for k, s in states:
            sum_of_aa = np.sum([a for a in s.values()]).astype(np.float)
            for i in range(NUM_EMISSIONS):
                emission_prob[k][str(i+1)] = s[str(i+1)] / sum_of_aa
                # e[counter, aa - 1] = s[str(aa)] / sum_of_aa
            counter += 1

        print(emission_prob)
        return emission_prob


if __name__ == "__main__":
    hmm = HMM(1, 2)
