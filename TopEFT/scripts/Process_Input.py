#import numpy as np
import math

def process_input(infile):
    readfile = open(infile, 'r')

    readdata = []
    for line in readfile:
        line = line.lstrip('_')
        line = line.rstrip('\n')
        #line = line.replace('ttH','ttX')
        line = line.split('&')
        readdata.append(line)

    #Lists of names
    categories = ["C_"+cat for cat in readdata[0][1:]] #e.g. 2los_ee_2j_1b
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

    #Find the most sensitive categories to the process signal strength by S/sqrt(B)
    #Summing all signals together to make life easier
    #Only use the top 20 categories
    #In the event of no background, ignore them
    sig_names = proc_names[-7:] #Signal Processes
    bkgd_names = proc_names[:-7] #Background Processes
    SorB_arr = []
    for cat in categories:
        bkgd = 0.
        sig = 0.
        for proc in bkgd_names:
            bkgd += proc_dict[proc,cat]
        for proc in sig_names:
            sig += proc_dict[proc,cat]
        SorB = sig/math.sqrt(bkgd) if bkgd!=0. else 0. # Some categories have no background.
        SorB_arr.append((cat,SorB))
    SorB_arr.sort(key=lambda tup: tup[1], reverse=True)
    SorB_arr = SorB_arr[:20]
    categories_best = [tuple[0] for tuple in SorB_arr]


    #return(categories, proc_names, proc_dict, sys_types, sys_dict)
#    return(categories[:50], proc_names[:19], proc_dict, sys_types, sys_dict)
#    return(categories[:50], ['ttW','ttZ','ttH','tZq','ttbar_dilepton'], proc_dict, sys_types, sys_dict)
    return(categories_best, ['ttW','ttZ','ttH','tZq']+bkgd_names, proc_dict, sys_types, sys_dict)
