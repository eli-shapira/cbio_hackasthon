import Protein
import re
import HMM, apply_hmm
import sys
import random
from sklearn.model_selection import train_test_split
from posterior import*
from globs import*


#patters
p_name = r"protein:)(\w+)(\s*)"
states_pattern = r"[oHST]+\s+"


aa_set = {'A', 'C' , 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', 'X'}
group_set = {'o', 'H', 'S', 'T'}


def replace_symbols(line):
    line = line.replace('H', 'A')
    line = line.replace('S', 'B')
    line = line.replace('o', 'O')
    return line


def parse_file(path):
    """
    cleans the vm file lines from garbage (didnt clean the spaces)
    :param path:
    :return:
    """

    #states
    aa_state = False
    group_state= False
    states_state = False
    keywords_state = False

    proteins = []
    name = "No one"
    aa_seq= ""
    group_seq= ""
    structure = ""
    keywords= ""

    #print(name)

    with open(path) as file:
        for line in file.readlines():
            #print(line)
            line = line.strip()
            if aa_state:
                if all(aa in aa_set for aa in line):
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
                    structure += replace_symbols(line)
            elif states_state:
                if re.fullmatch(states_pattern,line):  # maybe we can make more efficient with adding a special symbol like '$' at the beginning of each states line
                    structure += replace_symbols(line)
                else:
                    states_state = False
                    keywords_state = True
                    keywords += line

            elif keywords_state:
                if line.startswith('++'):
                    keywords_state = False
                    # print("I'm Here")
                    # print(name)
                    # print(aa_seq)
                    # print(group_seq)
                    # print(structure)
                    assert group_seq.isdigit(), 'assert 1 ' + group_seq
                    assert all(s in 'ABTO' for s in structure), 'assert 2'
                    protein = Protein.Protein(name, aa_seq, group_seq, keywords, structure)
                    if not protein.to_drop:
                        proteins.append(protein)
                    name = None
                    aa_seq= ""
                    group_seq= ""
                    structure = ""
                    keywords= ""
                else:
                    keywords += line

            elif line.startswith('protein:'):
                name = line[8:]
                print(name)
                aa_state = True

            else:
                raise ValueError("Wrong string pattern!")

    return proteins


def evaluate(true_structure, our_prediciton):
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


if __name__ == '__main__':
    proteins = parse_file('prot_data_human')
    print(len(proteins))
    # print([p.name for p in proteins])
    # for p in p_list:
    #     print(p.keywords)

    p_train = proteins[:4000]
    p_test = proteins[4000:]

    emissions = init_emissions_group(p_train, True)
    transitions = init_transitions(p_train, True)
    emissions_matrix = emission_to_matrix(emissions)
    transitions_matrix = transition_to_matrix(transitions)

    # seq = "0112233446"
    # posterior = calculate_posterior(seq, transitions_matrix, emissions_matrix)
    # print(posterior, "\n")

    true = []
    pred = []
    for p in p_test:
        matrix = calculate_posterior_group(p.group_seq, transitions_matrix, emissions_matrix)
        trace = trace_states(p.group_seq, matrix)
        # matrix, trace = calculate_viterbi(p.group_seq, transitions_matrix, emissions_matrix)
        trace = Protein.revert_structure3(trace)
        pred.append(trace)
        true.append(p.structure)
        try:
            print(trace[:60])
            print(p.structure[:60])
            print("======================================")
        except IndexError:
            continue

    err = evaluate(true, pred)
    print(err)

    # test_p = p_list[-2:]
    # train_p = p_list[:-2]
    # hmm = apply_hmm.HMM(train_p, test_p)
