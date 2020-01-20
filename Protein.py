class Protein:

    def __init__(self, name, aa_seq, group_seq, keywords, structure):
        self.name = name
        self.aa_seq = aa_seq[:-1]
        self.len = len(self.aa_seq)
        self.group_seq = group_seq[:-1]
        self.group_seq = '0' + self.group_seq + '1'
        self.keywords = keywords
        self.structure = self.replace_states(structure)[:-1]
        self.to_drop = False
        other_count = self.structure.count('O')
        if other_count / self.len > 0.85:
            self.to_drop = False
        self.structure = 's' + self.structure + 'e'

    def replace_states(self, s):
        return s.replace('H', 'A').replace('S', 'B').replace('o', 'O')
