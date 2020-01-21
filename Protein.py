def add_to_struct3(addition):
    if addition[0] == 'O':
        return addition
    elif addition[0] == 'A':
        return '@' + 'a' + addition[2:]
    elif addition[0] == 'B':
        return '&' + 'b' + addition[2:]
    elif addition[0] == 'T':
        return '!' + 't' + addition[2:]

def build_structure3(struct):
    struct3 = ''
    check = ''
    last_struct_seq = struct[0]
    for i in range(1, len(struct)):
        if last_struct_seq[-1] == struct[i]:
            last_struct_seq += struct[i]
        else:
            check += last_struct_seq
            struct3 += add_to_struct3(last_struct_seq)
            last_struct_seq = struct[i]
    check += last_struct_seq
    struct3 += add_to_struct3(last_struct_seq)
    return struct3

def revert_structure3(struct3):
    struct3 = struct3.replace('@', 'A')
    struct3 = struct3.replace('a', 'A')
    struct3 = struct3.replace('&', 'B')
    struct3 = struct3.replace('b', 'B')
    struct3 = struct3.replace('!', 'T')
    struct3 = struct3.replace('t', 'T')
    return struct3

class Protein:

    def __init__(self, name, aa_seq, group_seq, keywords, structure):
        self.name = name
        self.aa_seq = aa_seq[:-1]
        self.aa_seq = '0' + self.aa_seq + '1'
        self.len = len(self.aa_seq)
        self.group_seq = group_seq[:-1]
        self.group_seq = '0' + self.group_seq + '6'
        self.keywords = keywords.split(';')
        self.keywords[-1] = self.keywords[-1].split('.\n')[0]
        self.structure = self.replace_states(structure)[:-1]
        self.to_drop = False
        other_count = self.structure.count('O')
        if other_count / self.len > 0.85:
            self.to_drop = True
        self.structure = 's' + self.structure + 'e'
        self.to_3_states()

    def replace_states(self, s):
        return s.replace('H', 'A').replace('S', 'B').replace('o', 'O')

    def to_3_states(self):
        # print(self.structure)
        self.structure = 's' + build_structure3(self.structure[1:-1]) + 'e'

    def to_1_states(self):
        # print(self.structure)
        self.structure = 's' + revert_structure3(self.structure[1:-1]) + 'e'
