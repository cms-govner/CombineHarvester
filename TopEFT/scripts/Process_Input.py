#import numpy as np

def process_input(infile):
    readfile = open(infile, 'r')

    readdata = []
    for line in readfile:
        line = line.lstrip('_')
        line = line.rstrip('\n')
        line = line.split('&')
        readdata.append(line)

    #Lists of names
    categories = readdata[0][1:] #e.g. 2los_ee_2j_1b
    proc_names=[] #e.g. ttH
    sys_types=[] #e.g. Lumi

    proc_dict={} #process,category:nominal rate
    sys_dict={} #process,category:{systematic type:rate}

    for line in readdata[1:]:
        #Fill nominal value dict
        if ':' not in line[0]:
            proc_names.append(line[0])
            for cat_idx, cat_name in enumerate(categories):
                proc_dict.update({(line[0],cat_name):float(line[cat_idx+1])})
        #Fill systematic value dict
        else:
            sys=line[0].split(':')
            if sys[1] not in sys_types:
                sys_types.append(sys[1])
            for cat_idx, cat_name in enumerate(categories):
                if (sys[0],cat_name) not in sys_dict:
                    sys_dict[(sys[0],cat_name)]={}
                sys_dict[(sys[0],cat_name)].update({sys[1]:float(line[cat_idx+1])})

    return(categories, proc_names, proc_dict, sys_types, sys_dict)
