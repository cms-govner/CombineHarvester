#import numpy as np
import math
import ROOT
ROOT.gSystem.Load('$CMSSW_BASE/src/CombineHarvester/TopEFT/interface/TH1EFT_h.so')

def process_input(infile):
    readfile = open(infile, 'r')

    readdata = []
    for line in readfile:
        line = line.lstrip('_')
        line = line.rstrip('\n')
        line = line.split('&')
        readdata.append(line)

    #Lists of names
    categories = ["C_"+cat for cat in readdata[0][1:]] #e.g. 2los_ee_2j_1b
    proc_names=[] #e.g. ttH
    sys_types=[] #e.g. Lumi

    proc_dict={} #process,category:nominal rate
    sys_dict={} #process,category:{systematic type:ratio to nominal rate}

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
    #In the event of no background, ignore the category
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
    return(categories_best, ['ttH','tZq']+bkgd_names, proc_dict, sys_types, sys_dict)

def process_root_input(infile,operator):
    readfile = ROOT.TFile.Open(infile)

    #List of processes to look for
    #proc_known = ['singlet_tWchan','singletbar_tWchan','singlet_tchan','singletbar_tchan','singletop_schan','WZ','ZZ','WW','WWW','WWZ','WZZ','ZZZ','DYlowM','DY','WJets','ttJets','ttW','ttZ','ttH','tZq']
    sig_base = ['ttH','tllq','ttll','ttlnu']
    proc_known = ['charge_flips','fakes','WZ','ZZ','WW','WWW','WWZ','WZZ','ZZZ']
    data_known = ['data_doubleEle','data_muonEle','data_doubleMu','data_singleEle','data_singleMu']

    #Lists of names to pass on
    categories=[] #e.g. 2los_ee_2j_1b
    data_names=[] #e.g. doubleEle
    proc_names=[] #e.g. ttH
    sys_types=[] #e.g. Lumi

    #Dicts for rates
    data_dict={} #process,category:rate
    proc_dict={} #process,category:nominal rate
    sys_dict={} #process,category:{systematic type:ratio to nominal rate}


    #Main parsing loop
    for key in readfile.GetListOfKeys():
        hist = readfile.Get(key.GetName())

        #Get categorical information from histogram name
        histname = hist.GetName().split('.')
        category,systematic,process = '','',''
        if(len(histname)==3): [category,systematic,process] = histname
        if(len(histname)==2): [category,process] = histname

        #Debug
        #if process == 'ttll_cpt':
        #    for bin in range(1,hist.GetNbinsX()):
        #        category_njet = 'C_{0}_{1}j'.format(category,bin)
        #        bin_yield = round(hist.GetBinContent(1+bin,ROOT.WCPoint()),4)
        #        print process,category_njet,bin_yield
        #if process == 'ttZ':
        #    for bin in range(1,hist.GetNbinsX()):
        #        category_njet = 'C_{0}_{1}j'.format(category,bin)
        #        bin_yield = round(hist.GetBinContent(1+bin,ROOT.WCPoint()),4)
        #        print process,category_njet,bin_yield

        #Logic for data
        if process in data_known:
            if process not in data_names: data_names.append(process.replace('data_',''))
            for bin in range(1,hist.GetNbinsX()):
                category_njet = 'C_{0}_{1}j'.format(category,bin)
                if category_njet not in categories: categories.append(category_njet)
                bin_yield = round(hist.GetBinContent(1+bin,ROOT.WCPoint()),4)
                data_dict.update({(process,category_njet):bin_yield})

        #Process yields only for the specificied operator
        sig_known = [sig + "_" + operator for sig in sig_base]
        if process not in proc_known+sig_known: continue
        if process in sig_known:
            #print "Process:",process,"Hist:",hist.GetName()
            process = process.rsplit("_",1)[0]
        #process = process.replace('tllq','tZq')
        #process = process.replace('ttll','ttZ')
        #process = process.replace('ttlnu','ttW')

        #Logic for the systematic histograms
        if systematic != '':
            #print "Systematic Hist:",hist.GetName()
            if systematic not in sys_types: sys_types.append(systematic)

            for bin in range(1,hist.GetNbinsX()): #Don't look at 0jet bins
                #Check category exists. This is probably redundant as if it doesn't, it should give an error when calculating the ratio
                category_njet = 'C_{0}_{1}j'.format(category,bin)
                if category_njet not in categories: categories.append(category_njet)

                bin_yield = round(hist.GetBinContent(1+bin,ROOT.WCPoint()),4)
                bin_ratio = 1.0
                #Calculate ratio to nominal
                #If the nominal yield is zero, set the ratio to 1.0
                #If the systematic yield is 0, set the ratio to 0.0001 (Combine doesn't like 0)
                if proc_dict[(process,category_njet)] != 0:
                    bin_ratio = bin_yield/proc_dict[(process,category_njet)]
                    bin_ratio = max(bin_ratio,0.0001)
                #Create sys_dict key if it doesn't exit; can't edit a dict object that doesn't exist yet
                if not sys_dict.has_key((process,category_njet)):
                    sys_dict[(process,category_njet)] = {}
                sys_dict[(process,category_njet)][systematic] = bin_ratio

        #Logic for nominal rates
        else:
            #print "Nominal Hist:",hist.GetName()
            if process not in proc_names: proc_names.append(process)
            for bin in range(1,hist.GetNbinsX()):
                category_njet = 'C_{0}_{1}j'.format(category,bin)
                if category_njet not in categories: categories.append(category_njet)
                bin_yield = round(hist.GetBinContent(1+bin,ROOT.WCPoint()),4)
                proc_dict.update({(process,category_njet):bin_yield})

    # Only analyze categories with at least a few background events to prevent negative yields
    categories_nonzero = []
    for cat in categories:
        bkgd = 0.
        for proc in proc_names[:-4]:
            bkgd += proc_dict[proc,cat]
        if bkgd >= 0.0001: categories_nonzero.append(cat)
        #else:
            #print "Background less than 1!"

    return(categories_nonzero, data_names, data_dict, proc_names, proc_dict, sys_types, sys_dict)

        
