from globs import *
from parser import *
from trainer import *
from algos import *


def main():

    PATH = 'data/prot_data_human'
    proteins = parse_file(PATH)

    p_train = proteins[:10000]
    p_test = proteins[501:800]

    transitions = init_transitions(p_train)
    emissions_bi = init_emissions_bi(p_train)
    #emissions = init_emissions(p_train)
    
    true = [p.structure for p in p_test]
    seqs = [p.group_seq for p in p_test]
    print("testing")
    pred_bi = batch_posterior(seqs, transitions, emissions_bi, True)
    #pred = batch_posterior(seqs, transitions, emissions)
    #pred_viterbi = batch_viterbi(seqs, transitions, emissions)
    pred = pred_bi

    for i in range(len(seqs)):
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
