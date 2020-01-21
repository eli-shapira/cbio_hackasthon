import matplotlib
import numpy as np
import matplotlib.pyplot as plt
# sphinx_gallery_thumbnail_number = 2



#Can show a similiary matrix , emission , transition.
def show_heat_map(accuracy_2d,x_axis,y_axis,name,color='r'):
    """
    Recieves 2d np array , and axis names and prints the matrix with the heat method.
    :param accuracy_2d:
    :param x_axis:
    :param y_axis:
    :param color: the color of the text
    :return:
    """
    fig, ax = plt.subplots()
    im = ax.imshow(accuracy_2d,cmap='Blues_r')
    # We want to show all ticks...
    ax.set_xticks(np.arange(len(x_axis)))
    ax.set_yticks(np.arange(len(y_axis)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(x_axis)
    ax.set_yticklabels(y_axis)
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor"     )
    # Loop over data dimensions and create text annotations.
    for i in range(len(y_axis)):
        for j in range(len(x_axis)):
            text = ax.text(j, i, accuracy_2d[i, j],
                           ha="center", va="center", color=color)
    ax.set_title(name)
    fig.tight_layout()
    plt.savefig(name)
    plt.show()





def create_sizes(protiens):
    """
    Recieves protiens and count the distrubtions of the states , then shows a pie chart of them
    :param protiens:
    :return:
    """
    counter_a=0
    counter_b=0
    counter_t=0
    counter_o=0

    for seq in protiens:
        for char in seq:
            if char=='A':
                counter_a+=1
            if char=='B':
                counter_b+=1
            if char=='T':
                counter_t+=1
            if char=='O':
                counter_o+=1
    return [counter_a,counter_b,counter_t,counter_o]

def create_pie_chart(protiens,protien_name):
    sizes=create_sizes(protiens)
    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
    labels = 'Alpha', 'Beta', 'Turn', 'Others'
   # only "explode" the 2nd slice (i.e. 'Hogs')
    colors = ['red', 'blue', 'green', 'grey']

    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
            shadow=True, startangle=90,colors=colors)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig("Second_dist_pie"+protien_name+".png")
    plt.show()
create_pie_chart(["TTTTTTOOAOBABABABABA"],"dAVID")
from globs import *

def evaluate( true_structure, our_prediciton , name):
    """
    Evaluates the accuracy btw the true and our pred , also shows a heat confusion matrix which indicate what's
    the most popular error.
    :param true_structure:
    :param our_prediciton:
    :return:
    """
    matches = 0
    total = 0
    confusion_matrix = np.zeros((4, 4))
    for strc_true, strc_our in zip(true_structure, our_prediciton):

        for i in range(len(strc_true)):
            char_our = strc_our[i]
            char_true = strc_true[i]
            if strc_our[i] == strc_true[i]:
                matches += 1
            else:
                confusion_matrix[STATE_TO_INDEX[char_our]-1, STATE_TO_INDEX[char_true]-1] += 1
            total += 1
    show_heat_map(confusion_matrix, ["Alpha","Beta","Turn","Other"],["Alpha","Beta","Turn","Other"],"confusion matrix for : "+name)
    return matches / total

from posterior import *


def evluate_species(organisem_1_ems,organisem_1_trans,organisem_2,orgs_names):


    true = []
    pred = []
    for p in organisem_2:
        matrix = calculate_posterior_group(p.group_seq, organisem_1_trans, organisem_1_ems)
        trace = trace_states(p.group_seq, matrix)
        # matrix, trace = calculate_viterbi(p.group_seq, transitions_matrix, emissions_matrix)
        pred.append(trace)
        true.append(p.structure)
    err = evaluate(true, pred,orgs_names)
    print(err)
    return err

def evaluate_all_species(thetas, species_protiens):
    acc_matrix=np.zeros((len(species_protiens),len(species_protiens)))

    i=0
    names=["Human", "Yeast", "Dolphin", "Fruit fly"]
    for theta,name_1  in zip(thetas,names):
        j=0
        for specie_protiens,name_2 in zip(species_protiens,names):
            amount_of_p=len(species_protiens)
            acc_matrix[i,j]=evluate_species(theta[0],theta[1],specie_protiens[int(amount_of_p*0.9):],name_1+" "+ name_2)

            j+=1
        i+=1
    acc_matrix=np.round(acc_matrix,2)
    print(acc_matrix)
    show_heat_map(acc_matrix,["Human","Yeast"],["Human","Yeast"],"Posterior acc")
    return acc_matrix
def get_theta(species_protiens):


    theta_list=[]
    aa=["Human","Yeast","Dolphin","Fruit fly"]
    for spiece_protiens,name in zip(species_protiens,aa):
        emissions = init_emissions_group(spiece_protiens, True)
        transitions = init_transitions(spiece_protiens, True)
        emissions_matrix = emission_to_matrix(emissions)
        transitions_matrix = transition_to_matrix(transitions)
        states_strings=["start","Alpha","Beta","Turn","Other","End"]
        show_heat_map(np.round(np.exp(transitions_matrix),4),states_strings,states_strings,"Transition Matrix for "+name)
        # show_heat_map(emissions_matrix
        curr_theta=(emissions_matrix,transitions_matrix)
        theta_list.append(curr_theta)
    return theta_list

def postier_accuracy_between_species(species_protiens):
    """
    Receives protiens of different orgs , t
    :param species_protiens:
    :return:
    """
    all_thetas=get_theta(species_protiens)

    evluate_species(all_thetas[0][0],all_thetas[0][1],species_protiens[1],"Human")



def calculate_error(y ,train,protien,protien_name):
    if len(y)>300:
        y=y[:300]
    x=np.arange(len(y))
    plt.ylabel("posterior error")
    plt.xlabel("protein  sequences")
    plt.title("posterior error of 300 sequences trained on the  " +train+" organism \n  and tested on the " + protien + " organism")

    plt.plot(x,y)
    plt.savefig("accChart.png")

    plt.show()



