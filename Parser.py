from Protein import Protein
import re
import HMM, apply_hmm
import sys
import random
from sklearn.model_selection import train_test_split


#patters
p_name = r"protein:)(\w+)(\s*)"
states_pattern = r"[oHST]+\s+"


aa_set = {'A', 'C' , 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y'}
group_set = {'o', 'H', 'S', 'T'}


def parse_file(path):
    """
    cleans the vm file lines from garbage (didnt clean the spaces)
    :param path:
    :return:
    """

    #states
    aa_state = False
    group_state = False
    states_state = False
    keywords_state = False

    proteins = []
    name = "No one"
    aa_seq = ""
    group_seq = ""
    structure = ""
    keywords = ""

    file = open(path)
    for line in file:
        if aa_state:
            if line[0] in aa_set:
                aa_seq += line
            else:
                aa_state = False
                group_state = True
                group_seq += line

        elif group_state:
            char = line[0]
            if char.isdigit():
                group_seq += line
            else:
                group_state = False
                states_state = True
                structure += line
        elif states_state:
            if re.fullmatch(states_pattern, line):  # maybe we can make more efficient with adding a special symbol like '$' at the beginning of each states line
                structure += line
            else:
                states_state = False
                keywords_state = True
                keywords += line

        elif keywords_state:
            if line.startswith('++'):
                keywords_state = False
                protein = Protein(name, aa_seq, group_seq, keywords, structure)
                if not protein.to_drop:
                    proteins.append(protein)
                name = None
                aa_seq = ""
                group_seq = ""
                structure = ""
                keywords = ""
            else:
                keywords += line

        elif line.startswith('protein:'):
            name = line[8:]
            aa_state = True

        else:
            raise ValueError("Wrong string pattern!")

    return proteins


if __name__ == '__main__':
    p_list = parse_file('prot_data_yeast.txt')
    # for p in p_list:
    #     print(p.keywords)

    ind = int(len(p_list) / 2)
    # indeces = [i for i in range(len(p_list))]
    # train_i = random.sample(indeces, 300)
    # test_i = i


    print(len(p_list))
    test_p = p_list[2300:]
    train_p = p_list[:2300]
    hmm = apply_hmm.HMM(train_p, test_p)