import copy
import numpy as np
try:
    from scipy.misc import logsumexp
except:
    from scipy.special import logsumexp
from numpy import logaddexp

BOLD = '\033[01m'
RED = "\033[91m"
CYAN = "\033[96m"
YELLOW = "\033[93m"
PURPLE = "\033[35m"
END = "\033[00m"

SIMPLE_STATES = ['s', 'A', 'B', 'T', 'O', 'e']
MULTI_STATES = ['s', '@', 'a', 'A', '&', 'b', 'B', '!', 't', 'T', 'O', 'e']
STATES = MULTI_STATES
STATE_TO_INDEX = {l: STATES.index(l) for l in STATES}
NUM_STATES = len(STATES)

# start, negative, positive, polar, hydrophobic, special, end
EMISSIONS_GROUP = ['0', '1', '2', '3', '4', '5', '6']
EMISSIONS_RAW = ['0', 'G', 'M', 'L', 'N', 'A', 'S', 'F', 'Y', 'T', 'I', 'W', 'Q', 'C', 'P', 'D', 'V', 'E', 'H', 'K', 'R', '1']
EMISSIONS = EMISSIONS_GROUP
NUM_EMISSIONS = len(EMISSIONS)

def print_emissions(dict):
    print(BOLD+RED+"\t", end='')
    for em in EMISSIONS:
        print(em+"\t", end='')
    print("\n"+END, end="")
    for st in STATES:
        print(BOLD+RED+st+"\t"+END, end="")
        for em in EMISSIONS:
            print(CYAN+str(dict[st][em])+"\t"+END, end='')
        print("\n", end="")
    print("\n", end="")

def print_transitions(dict):
    print(BOLD+RED+"\t", end='')
    for i in STATES:
        print(i+"\t", end='')
    print("\n"+END, end="")
    for i in STATES:
        print(BOLD+RED+i+"\t", end="")
        for j in STATES:
            print(CYAN+str(dict[i][j])+"\t"+END, end='')
        print("\n", end="")
    print("\n", end="")
