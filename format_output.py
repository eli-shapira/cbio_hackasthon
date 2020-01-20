BOLD = '\033[01m'
RED = "\033[91m" 
CYAN = "\033[96m"
YELLOW = "\033[93m"
PURPLE = "\033[35m"
END = "\033[00m"


MAX_CHARS = 50

from termcolor import cprint

color_dict = { 'A':'green', 'B':'blue', 'O':'yellow', 'T':'magenta'}

def print_line(a, b, i, L):

    for j in range(i, L):
        if a[j] != b[j]:
            cprint(a[j],color_dict[a[j]], 'on_red',end ="")
        else:
            cprint(a[j],color_dict[a[j]], end ="")

    print()
    for j in range(i, L):
        if a[j] != b[j]:
            cprint(b[j],color_dict[b[j]], 'on_red' ,end ="") #color_dict[b[j]]
        else:
            cprint(b[j],color_dict[b[j]], end ="") #color_dict[b[j]]

    print('\n')

def print_seq(a,b):

    if len(a) != len(b):
        raise Exception("The length of the 2 sequences should be identical!")

    i=0
    while(i<len(a)-MAX_CHARS):
        # print(a[i:i+MAX_CHARS])
        # print(b[i:i+MAX_CHARS] + "\n")
        print_line(a,b,i,i+MAX_CHARS),

        i += MAX_CHARS

    print_line(a,b,i,len(a))




str1 = "OOOOOOOOOOTTTTTATTTTBBBBOOABBBAAABAOOOOOOATATOOOOOOOOOOTTTTTATTTTBBBBOOABBBAAABAOOOOOOATAT"
str2 = "OOOOOOOOOOTTTTTTTTTTBBBBBBBBBBAAAAAOOOOOOATATOOOOOOOOOOTTTTTATTTTBBBBOOABBBAAABAOAOOOOATAT"



print_seq(str1, str2)
