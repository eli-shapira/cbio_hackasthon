from Protein import Protein
import re
import sys
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

    file = open(path)
    for line in file:
        #print(line)
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
            if re.fullmatch(states_pattern,line):  # maybe we can make more efficient with adding a special symbol like '$' at the beginning of each states line
                structure += line
            else:
                states_state = False
                keywords_state = True
                keywords += line

        elif keywords_state:
            if line.startswith('++'):
                keywords_state = False
                print("I'm Here")
                print(name)
                print(aa_seq)
                print(group_seq)
                print(structure)
                proteins.append(Protein(name, aa_seq, group_seq, keywords,structure))
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

if __name__ == '__main__':
    path = sys.argv[1]
    p_list = parse_file(path)
    print(len(p_list))
    for p in p_list:
        p.printP()

    # str_list = [ "sadasd", "dasd", "123213124", "safsaf", "safasa"]
    # for str1 in str_list:
    #     if str1[0].isdigit():
    #         print(str1)