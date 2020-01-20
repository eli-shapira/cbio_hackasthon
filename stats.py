import argparse
import matplotlib.pyplot as plt
from collections import Counter
import copy
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('to_stats', help='File path with the data to do stats on')
    parser.add_argument('ds', help='the dataset to indicate this group (e.g. yeast)')
    parser.add_argument('action', help='the action to perform to indicate this group', choices=['lengths', 'aa_dist', 'aag_dist'])
    return parser.parse_args()


def count_length(struct_seq, lengths):
    prev_state = 'o'
    last_len = 0
    for i in range(len(struct_seq)):
        cur_state = struct_seq[i]
        if prev_state != cur_state:
            if prev_state != 'o':
                if last_len not in lengths[prev_state]:
                    lengths[prev_state][last_len] = 0
                lengths[prev_state][last_len] += 1
            last_len = 0
        else:
            last_len += 1
        prev_state = cur_state
    return lengths


def count_length_file(path, ds):
    lengths = {'T': {}, 'H': {}, 'S': {}}
    with open(path, 'r') as fd:
        line = fd.readline().strip()
        while line:
            if all(k in ['o', 'T', 'H', 'S'] for k in Counter(line).keys()):
                lengths = count_length(line, lengths)
            line = fd.readline().strip()
    print(lengths)
    title = ' lengths histogram %s' % ds
    plt.figure('TURN' + title)
    plt.bar(list(lengths['T']), lengths['T'].values(), color='g')
    plt.savefig('stats/TURN' + title + '.png')
    plt.figure('ALPHA HELIX' + title)
    plt.bar(list(lengths['H']), lengths['H'].values(), color='b')
    plt.savefig('stats/ALPHA HELIX' + title + '.png')
    plt.figure('BETA STRAND' + title)
    plt.bar(list(lengths['S']), lengths['S'].values(), color='r')
    plt.savefig('stats/BETA STRAND' + title + '.png')


def count_aag_dist(aag_struct_seq, counts):
    for t in aag_struct_seq:
        counts[t[1]][t[0]] += 1
    return counts


def count_aag_dist_file(path, ds):
    aag = {'1': 0, '2': 0, '3': 0, '4': 0, '5': 0}
    counts = {'T': copy.deepcopy(aag), 'H': copy.deepcopy(aag), 'S': copy.deepcopy(aag), 'o': copy.deepcopy(aag)}
    with open(path, 'r') as fd:
        line = fd.readline().strip()
        while line:
            if line.startswith('protein'):
                while not line.startswith('+'):
                    if all(k in ['o', 'T', 'H', 'S'] for k in Counter(line).keys()):
                        tags = line
                    elif line.isdigit():
                        aags = line
                    line = fd.readline().strip()
                counts = count_aag_dist(zip(aags, tags), counts)
            line = fd.readline().strip()
    df = pd.DataFrame(counts)
    print(df)
    df_struct_sum = df.sum(axis=0)
    print(df_struct_sum)
    df_struct_norm = df / df_struct_sum # normalized by count of structures
    print(df_struct_norm)
    print(df_struct_norm.sum(axis=0))
    df_aag_sum = df.sum(axis=1)
    print(df_aag_sum)
    df_aag_norm = (df.T / df_aag_sum).T # normalized by count of AA groups
    print(df_aag_norm)
    print(df_aag_norm.sum(axis=1))
    
    title = "count vs structure (log) %s" %ds
    plot = df.T.plot(logy=True, title=title)
    fig = plot.get_figure()
    fig.savefig('stats/' + title + '.png')

    title = "count vs AA group (log) %s" %ds
    plot = df.plot(logy=True, title=title)
    fig = plot.get_figure()
    fig.savefig('stats/' + title + '.png')

    title = "count vs AA group normalized by sum per structure %s" %ds
    plot = df_struct_norm.plot(title=title)
    fig = plot.get_figure()
    fig.savefig('stats/' + title + '.png')

    title = "count vs structure normalized by sum per structure %s" %ds
    plot = df_struct_norm.T.plot(title=title)
    fig = plot.get_figure()
    fig.savefig('stats/' + title + '.png')

    title = "count vs AA group normalized by sum per AA group %s" %ds
    plot = df_aag_norm.plot(logy=True, title=title)
    fig = plot.get_figure()
    fig.savefig('stats/' + title + '.png')

    title = "count vs structure normalized by sum per AA group %s" %ds
    plot = df_aag_norm.T.plot(logy=True, title=title)
    fig = plot.get_figure()
    fig.savefig('stats/' + title + '.png')

    plt.show()
    
    


if __name__ == "__main__":
    args = parse_args()
    if args.action == 'lengths':
        count_length_file(args.to_stats, args.ds)
    if args.action == 'aag_dist':
        count_aag_dist_file(args.to_stats, args.ds)
