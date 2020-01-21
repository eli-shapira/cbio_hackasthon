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

        OMIT = {"Reference proteome", "3D-structure", "Direct protein sequencing"}
        self.keywords = keywords.split(';')
        self.keywords[-1] = self.keywords[-1].split('.\n')[0]
        for i in range(len(self.keywords)):
            self.keywords[i] = self.keywords[i].strip()
        self.keywords = [self.keywords[i] for i in range(len(self.keywords)) if self.keywords[i] not in OMIT]
        
        self.structure = self.replace_states(structure)[:-1]
        self.to_drop = False
        other_count = self.structure.count('O')
        if other_count / self.len > 0.8:
            self.to_drop = True
        self.structure = 's' + self.structure + 'e'
        self.is_3_states = False

    def replace_states(self, s):
        return s.replace('H', 'A').replace('S', 'B').replace('o', 'O')

    def evaluate_prediction(self, pred):
        self.prediction = pred
        matches = 0
        total = 0
        for strc_true, strc_our in zip(self.structure, self.prediction):
            for i in range(len(strc_true)):
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
        self.score = matches/total
        print("evaluated",self.name[:-1],":",self.score)
        return matches/total

    def to_3_states(self):
        if not self.is_3_states:
            self.structure = 's' + build_structure3(self.structure[1:-1]) + 'e'
        self.is_3_states = True

    def to_1_states(self):
        self.is_3_states = False
        self.structure = 's' + revert_structure3(self.structure[1:-1]) + 'e'
