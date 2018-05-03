#import numpy as np

def process_input(infile):
    readfile = open(infile, 'r')

    readdata = []
    for line in readfile:
        line = line.rstrip('\n')
        line = line.split('&')
        readdata.append(line)

    categories = readdata[0][1:]
    proc_values = [float(line[1]) for line in readdata[1:]]

    #proc_dict={}
    sys_dict={}
    proc_names=[]
    sys_types=[]
    for line in readdata[1:]:
        if ':' not in line[0]:
            proc_names.append(line[0])
            #for cat_idx, cat_name in enumerate(categories):
            #    proc_dict.update({(line[0],cat_name):line[cat_idx+1]})
        else:
            sys=line[0].split(':')
            #print sys
            if sys[1] not in sys_types:
                sys_types.append(sys[1])
            for cat_idx, cat_name in enumerate(categories):
                #print cat_idx, cat_name
                #sys_dict[(sys[0],cat_name)] = {sys[1]:line[cat_idx+1]}
                if (sys[0],cat_name) not in sys_dict:
                    sys_dict[(sys[0],cat_name)]={}
                sys_dict[(sys[0],cat_name)].update({sys[1]:line[cat_idx+1]})

    return(categories, proc_names, proc_values, sys_types, sys_dict)
