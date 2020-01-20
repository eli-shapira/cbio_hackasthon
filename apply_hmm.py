import numpy as np

class HMM:

    def __init__(self, train_p, test_proteins):
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


    def init_emission(self, train_p):
        """
        Init the emission matrix
        :param train_p:
        :return:
        """
        aa_groups = {key: 0 for key in range(1, 6)}
        O, H, S, T = aa_groups, aa_groups, aa_groups, aa_groups
        states = {'O': O, 'H': O, 'S': O, 'T': O}
        for p in train_p:
            for i in range(p.len):
                # states[p[i]] returns the dictionary of the state
                states[p.structure[i]][p.group_seq] += 1

        e = np.zeros((4, 5), dtype=np.float)
        for state in states:
            sum_of_aa = np.sum([c for c in t.values for t in state.values()]).astype(np.float)
            for aa in range(1, 6):
                e[state, aa] = state[aa] / sum_of_aa

        return e

