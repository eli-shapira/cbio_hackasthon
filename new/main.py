from globs import *
from parser import *
from trainer import *
from algos import *

def filter_by_keyword(proteins, kw):
    return [p for p in proteins if kw in p.keywords]
    
def print_sorted_dict(dct):    
    for kw in sorted(dct, key=dct.get)[::-1]:
        print(RED,dct[kw],END,kw)

def keyword_histogram(kw_lists):
    dct = {}
    for lst in kw_lists:
        for i in range(len(lst)):
            if lst[i] in dct.keys():
                dct[lst[i]] += 1
            else:
                dct[lst[i]] = 1
    print_sorted_dict(dct)

    return dct

def score_histogram(proteins, all_keywords):
    counter = {}
    for p in proteins:
        for kw in p.keywords:
            if kw in counter.keys():
                counter[kw] += p.score
            else:
                counter[kw] = p.score
    for kw in counter:
        counter[kw] = counter[kw] / all_keywords[kw]
    for kw in sorted(counter, key=counter.get)[::-1]:
        print(CYAN,all_keywords[kw],RED,counter[kw],END,kw)

    return counter


def main():

    PATH = 'data/prot_data_human'
    proteins = parse_file(PATH)
    print("len proteins:", len(proteins))

    p_train = proteins[:3000]
    p_test = proteins[100:200]


    transitions = init_transitions(p_train)
    #emissions_bi = init_emissions_bi(p_train)
    emissions = init_emissions(p_train)
    
    print("testing")
    kw_lists = [p.keywords for p in p_test]
    all_keywords = keyword_histogram(kw_lists)
    
    #pred_bi = batch_posterior(seqs, transitions, emissions_bi, True)
    #pred_viterbi = batch_viterbi(seqs, transitions, emissions)
    for p in p_test:
        pred = calculate_posterior(p.group_seq, transitions, emissions)
        p.evaluate_prediction(pred)

    score_histogram(p_test, all_keywords)



if __name__ == "__main__":
    main()
