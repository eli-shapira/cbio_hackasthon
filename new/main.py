from globs import *
from parser import *
from trainer import *
from algos import *
import matplotlib.pyplot as plt
import random
# from viz import calculate_error
from Protein import revert_structure3

def calculate_error(y1, y2, y3, name):

    # if len(y)>300:
    #     y = random.choices(y, k=300)
        # y=y[:300]
    min_len = min([len(y2), len(y3), len(y1)])
    x=np.arange(min_len)
    plt.ylabel("posterior error")
    plt.xlabel("protein  sequences")
    plt.ylim(top=1, bottom=0)
    plt.title("posterior accuracies on sequences - " + name)
    # ax = plt.figure()
    # ax.set_xticklabels([])
    plt.plot(x,y1[:min_len], label='short')
    plt.plot(x,y2[:min_len], label='mid')
    plt.plot(x,y3[:min_len], label='long')
    plt.savefig(name + "_accChart.png")

    plt.show()

def filter_by_keyword(proteins, kw):
    print("filtering proteins by keyword:",RED,kw,END)
    filtered = [p for p in proteins if kw in p.keywords]
    print("filtered",len(filtered),"from",len(proteins))
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
    # PATH = 'data/prot_data_yeast'
    # proteins += parse_file(PATH)
    print("len proteins:", len(proteins))

    #fil = filter_by_keyword(proteins, "Zinc-finger")
    #fil = filter_by_keyword(proteins, "Transcription regulation")
    #fil = filter_by_keyword(fil, "Transcription")
    fil = proteins

    idx = len(fil)//10
    p_train = fil[:idx*9]
    p_test = fil[idx*9+1::]
    print("training:",(idx*9), "test",str(len(fil)-idx*9))

    if 'a' in STATES:
        for p in proteins:
            p.to_3_states()


    transitions = init_transitions(p_train)
    #emissions_bi = init_emissions_bi(p_train)
    emissions = init_emissions(p_train)
    np.set_printoptions(precision=2, linewidth=250)
    print(np.exp(transitions))
    print(np.exp(emissions))
    # return
    if 'a' in STATES:
        for p in p_test:
            p.to_1_states()

    p_test = sorted(p_test, key=lambda x: x.len)
    
    true = [p.structure for p in p_test]
    seqs = [p.group_seq for p in p_test]
    print("testing")
    kw_lists = [p.keywords for p in p_test]
    all_keywords = keyword_histogram(kw_lists)
    
    #pred_bi = batch_posterior(seqs, transitions, emissions_bi, True)
    #pred_viterbi = batch_viterbi(seqs, transitions, emissions)
    acc = 0
    sum_len = 0
    for p in p_test:
        print(len(p.group_seq))
        pred = calculate_posterior(p.group_seq, transitions, emissions)
        # pred = calculate_viterbi(p.group_seq, transitions, emissions)
        # pred = revert_structure3(pred)
        # print(pred[:100])
        eva = p.evaluate_prediction(pred)
        acc += eva * p.len
        sum_len += p.len

    score_histogram(p_test, all_keywords)
    average_score(p_test)
    print(acc / sum_len)
    errors = [p.score for p in p_test]
    break_point = int(len(errors) / 3)
    err1 = sorted(errors[:break_point])
    err2 = sorted(errors[break_point:-break_point])
    err3 = sorted(errors[-break_point:])
    # min_len = min(err1, err2, )
    # calculate_error(errors, 'random')
    # calculate_error(err1, 'short')
    # calculate_error(err2, 'mid')
    # calculate_error(err3, 'long')
    calculate_error(err1, err2, err3, 'compare by sequence length')
    # calculate_error(errors[-300:], 'human and yeats', 'long, human and yeats', '')


if __name__ == "__main__":
    main()
