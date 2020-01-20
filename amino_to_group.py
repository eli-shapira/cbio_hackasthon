import sys

groups = {
    'special': ['C', 'G', 'P', 'U'],
    'negative': ['D', 'E'],
    'positive': ['H', 'K', 'R'],
    'polar': ['N', 'S', 'T', 'Q'],
    'hydro': ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
}

groups_number = {
    '1': ['D', 'E'],
    '2': ['H', 'K', 'R'],
    '3': ['N', 'S', 'T', 'Q'],
    '4': ['A', 'F', 'I', 'L', 'M', 'V', 'W', 'Y'],
    '5': ['C', 'G', 'P', 'U']
}

amino_to_group = {
    'A': '4',
    'B': '',
    'C': '5',
    'D': '1',
    'E': '1',
    'F': '4',
    'G': '5',
    'H': '2',
    'I': '4',
    'J': '',
    'K': '2',
    'L': '4',
    'M': '4',
    'N': '3',
    'O': '',
    'P': '5',
    'Q': '3',
    'R': '2',
    'S': '3',
    'T': '3',
    'U': '5',
    'V': '4',
    'W': '4',
    'X': '',
    'Y': '4',
    'Z': ''
}


def amino_seq_to_group(amino_seq):
    group_seq = ''
    for s in amino_seq:
        group_seq += amino_to_group[s]
    return group_seq


def transform_file(path):
    with open(path, 'r') as org_file, open(path + '_group', 'w') as dest_file:
        line = org_file.readline()
        while line:
            new_line = line
            if not line.startswith('+') and line.upper() == line:
                new_line = amino_seq_to_group(line.strip()) + '\n'
                assert len(new_line) == len(line)
            dest_file.write(new_line)
            line = org_file.readline()


if __name__ == "__main__":
    assert sys.argv.__len__() > 1
    transform_file(sys.argv[1])
