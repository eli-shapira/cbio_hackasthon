from globs import *
from parser import *
from trainer import *
from algos import *
from Protein import *
from random import shuffle

def print_50(str1, str2):
    i,j = 0,0
    while i < len(str1)-1:
        j = min(i+50, len(str1))
        print(str1[i:j])
        print(str2[i:j], '\n')
        i = j

def filter_by_keyword(proteins, kw):
    print("filtering proteins by keyword:",RED,kw,END)
    filtered = [p for p in proteins if kw in p.keywords]
    print("filtered",len(filtered),"from",len(proteins))
    return [p for p in proteins if kw in p.keywords]

def filter_by_name(proteins, name):
    for p in proteins:
        if p.name == name:
            return p

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

def average_score(proteins):
    count = 0
    for p in proteins:
        count += p.score
    avg = count / len(proteins)
    print("average score:",RED,avg,END)
    return count / len(proteins)

def main():

    PATH = 'data/prot_data_human'
    proteins = parse_file(PATH)
    shuffle(proteins)
    print("len proteins:", len(proteins))

    #fil = filter_by_keyword(proteins, "Zinc-finger")
    fil = filter_by_keyword(proteins, "Transcription regulation")
    fil = filter_by_keyword(fil, "Transcription")
    #fil = proteins

    idx = len(fil)//10
    p_train = fil[:idx*9]
    p_test = fil[idx*9+1::]
    print("training:",(idx*9), "test",(len(fil)-idx*9))

    if 'a' in STATES:
        for p in proteins:
            p.to_3_states()

    transitions = init_transitions(p_train)
    #emissions_bi = init_emissions_bi(p_train)
    emissions = init_emissions(p_train)

    np.savetxt('fam_em.out', emissions, delimiter=',')
    np.savetxt('fam_trans.out', transitions, delimiter=',')

    if 'a' in STATES:
        for p in p_test:
            p.to_1_states()
    
    print("testing")
    kw_lists = [p.keywords for p in p_test]
    all_keywords = keyword_histogram(kw_lists)
    
    #pred_bi = batch_posterior(seqs, transitions, emissions_bi, True)
    #pred_viterbi = batch_viterbi(seqs, transitions, emissions)
    
    for p in p_test:
        pred = calculate_posterior(p.group_seq, transitions, emissions)
        p.evaluate_prediction(revert_structure3(pred))
        if p.score > 0.8:
            print_50(p.structure, p.prediction)

    score_histogram(p_test, all_keywords)
    average_score(p_test)


if __name__ == "__main__":
    main()
