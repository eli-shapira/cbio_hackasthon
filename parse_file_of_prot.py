from itertools import groupby
import pprint
import argparse
from amino_to_group import amino_seq_to_group


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('to_parse', help='File path with the data as download')
    parser.add_argument('ds', help='the dataset to indicate this group (e.g. yeast)')
    return parser.parse_args()


def parse_one(prot_data):
    seq = ''
    structures = []
    str_dict = {}
    key_words = ''

    for i in range(len(prot_data)):
        if prot_data[i].startswith('SQ'):
            i += 1
            while not prot_data[i].startswith('//'):
                seq += prot_data[i].strip()
                i += 1
            seq = seq.replace(' ', '')
        elif prot_data[i].startswith('FT'):
            tag = '-'
            line = prot_data[i][2:].strip()
            if line.startswith('TURN'):
                line = line.replace('TURN', '').strip()
                tag = 'T'
            elif line.startswith('STRAND'):
                line = line.replace('STRAND', '').strip()
                tag = 'S'
            elif line.startswith('HELIX'):
                line = line.replace('HELIX', '').strip()
                tag = 'H'
            else:
                continue
            ind = line.split('..')
            start = int(ind[0])
            end = int(ind[1])
            for j in range(start, end+1):
                structures.append((tag, j))
                assert j not in str_dict
                str_dict[j] = tag
        elif prot_data[i].startswith('KW'):
            key_words += prot_data[i][2:].strip()
    tags = ''
    for k in range(1, len(seq) + 1):
        if k in str_dict:
            tags += str_dict[k]
        else:
            tags += 'o'
    return '%s\n%s\n%s\n%s\n++\n' % (seq, amino_seq_to_group(seq), tags, key_words)


def parse_file(path, ds):
    with open(path, 'r') as prot_file:
        lines = prot_file.readlines()
    prots = (x[1] for x in groupby(lines, lambda line: line.startswith("protein:")))
    with open('prot_data_' + ds, 'w') as tag_seq_file:
        for prot_lines in prots:
            prot_name = next(prot_lines)
            prot_data = [s.strip() for s in next(prots)]
            tag_seq_file.write(prot_name)
            tag_seq_file.write(parse_one(prot_data))


if __name__ == "__main__":
    args = parse_args()
    parse_file(args.to_parse, args.ds)
