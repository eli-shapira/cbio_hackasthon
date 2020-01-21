from globs import *
from parser import *
from trainer import *
from algos import *
from Protein import revert_structure3


def main():

    PATH = 'data/prot_data_human'
    proteins = parse_file(PATH)

    if 'a' in STATES:
        for p in proteins:
            p.to_3_states()

    p_train = proteins[:10000]
    p_test = proteins[501:800]

    transitions = init_transitions(p_train)
    emissions_bi = init_emissions_bi(p_train)
    #emissions = init_emissions(p_train)
    
    if 'a' in STATES:
        for p in p_test:
            p.to_1_states()
    true = [p.structure for p in p_test]
    seqs = [p.group_seq for p in p_test]
    print("testing")
    pred_bi = batch_posterior(seqs, transitions, emissions_bi, True)
    #pred = batch_posterior(seqs, transitions, emissions)
    #pred_viterbi = batch_viterbi(seqs, transitions, emissions)
    pred = pred_bi

    for i in range(len(seqs)):
        if 'a' in STATES:
            pred[i] = revert_structure3(pred[i])
        try:
            print(pred[i][:60])
            print(true[i][:60])
            print("======================================")
        except IndexError:
            continue

    err = evaluate(true, pred)
    print(err)





if __name__ == "__main__":
    main()
