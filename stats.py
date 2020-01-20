import argparse
import matplotlib.pyplot as plt
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('to_stats', help='File path with the data to do stats on')
    parser.add_argument('ds', help='the dataset to indicate this group (e.g. yeast)')
    parser.add_argument('action', help='the action to perform to indicate this group', choices=['lengths'])
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
    plt.savefig('TURN' + title + '.png')
    plt.figure('ALPHA HELIX' + title)
    plt.bar(list(lengths['H']), lengths['H'].values(), color='b')
    plt.savefig('ALPHA HELIX' + title + '.png')
    plt.figure('BETA STRAND' + title)
    plt.bar(list(lengths['S']), lengths['S'].values(), color='r')
    plt.savefig('BETA STRAND' + title + '.png')
    


if __name__ == "__main__":
    args = parse_args()
    if args.action == 'lengths':
        count_length_file(args.to_stats, args.ds)
