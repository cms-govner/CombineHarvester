import math
import ROOT
ROOT.gSystem.Load('$CMSSW_BASE/src/CombineHarvester/TopEFT/interface/TH1EFT_h.so')

def process_root_input(infile,operator):
    print "Setting up..."
    readfile = ROOT.TFile.Open(infile)

    #List of processes to look for
    #proc_known = ['singlet_tWchan','singletbar_tWchan','singlet_tchan','singletbar_tchan','singletop_schan','WZ','ZZ','WW','WWW','WWZ','WZZ','ZZZ','DYlowM','DY','WJets','ttJets','ttW','ttZ','ttH','tZq']
    sgnl_known = ['ttH','tllq','ttll','ttlnu']
    sgnl_histnames = [sgnl + '_' + operator for sgnl in sgnl_known]
    bkgd_known = ['charge_flips','fakes','WZ','ZZ','WW','WWW','WWZ','WZZ','ZZZ']
    data_known = ['data_doubleEle','data_muonEle','data_doubleMu','data_singleEle','data_singleMu']

    #Lists of names to pass on
    categories=[] #e.g. 2los_ee_2j_1b
    data_names=[] #e.g. doubleEle
    sgnl_names=[] #e.g. ttH
    bkgd_names=[] #e.g. ttH
    sys_types=[] #e.g. Lumi

    #Dicts for rates
    data_dict={} #process,category:rate
    nom_dict={} #process,category:nominal rate
    sys_dict={} #process,category:{systematic type:ratio to nominal rate}

    #For Plot
    plot_values=[[],[],[]]
    plot_sum = dict.fromkeys(['C_2lss_p_ee_1b_4j', 'C_2lss_p_ee_1b_5j', 'C_2lss_p_ee_1b_6j', 'C_2lss_p_ee_1b_7j', 'C_3l_ppp_1b_2j', 'C_3l_ppp_1b_3j', 'C_3l_ppp_1b_4j', 'C_3l_ppp_1b_5j', 'C_3l_ppp_1b_6j', 'C_3l_ppp_1b_7j'],0)

    print "Looping through histograms..."
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

        #Logic for data NOT VERIFIED
        if process in data_known:
            if process not in data_names: data_names.append(process.replace('data_',''))
            for bin in range(1,hist.GetNbinsX()):
                category_njet = 'C_{0}_{1}j'.format(category,bin)
                if category_njet not in categories: categories.append(category_njet)
                bin_yield = round(hist.GetBinContent(1+bin,ROOT.WCPoint()),4)
                data_dict.update({(process,category_njet):bin_yield})

        #Obtain signal yields only for the specificied operator
        if process not in bkgd_known+sgnl_histnames: continue
        if process in sgnl_histnames:
            process = process.rsplit("_",1)[0]
            if process not in sgnl_names: sgnl_names.append(process)
            #process = process.replace('tllq','tZq')
            #process = process.replace('ttll','ttZ')
            #process = process.replace('ttlnu','ttW')
        if process in bkgd_known:
            if process not in bkgd_names: bkgd_names.append(process)

        #Logic for the nominal histograms
        if systematic == '':
            #print "Nominal Hist:",hist.GetName()
            debug_count=0
            for bin in range(1,hist.GetNbinsX()):
                category_njet = 'C_{0}_{1}j'.format(category,bin)
                if category_njet not in categories: categories.append(category_njet)
                bin_yield = round(hist.GetBinContent(1+bin,ROOT.WCPoint()),4)
                nom_dict.update({(process,category_njet):bin_yield})
                debug_count += bin_yield
                if category_njet in ['C_2lss_p_ee_1b_4j', 'C_2lss_p_ee_1b_5j', 'C_2lss_p_ee_1b_6j', 'C_2lss_p_ee_1b_7j', 'C_3l_ppp_1b_2j', 'C_3l_ppp_1b_3j', 'C_3l_ppp_1b_4j', 'C_3l_ppp_1b_5j', 'C_3l_ppp_1b_6j', 'C_3l_ppp_1b_7j']:
                    plot_sum[category_njet] += bin_yield
                    if process in sgnl_known:
                        WCeq1 = round(hist.GetBinContent(1+bin,ROOT.WCPoint("EFTrwgt_"+operator+"_1",100.)),4)
                        plot_values[0].append(process)
                        plot_values[1].append(category_njet)
                        plot_values[2].append(WCeq1-bin_yield)

            #print process,category_njet,debug_count
            #if process in ['WW','WZ','ZZ']: print process,category_njet,debug_count
            #if process in ['WWW','WWZ','WZZ','ZZZ']: print category_njet,debug_count
            #if process in ['fakes']: print category_njet,debug_count
            #if process in ['charge_flips']: print category_njet,debug_count

        #Logic for systematic histograms
        else:
            #print "Systematic Hist:",hist.GetName()
            if systematic not in sys_types: sys_types.append(systematic)

            for bin in range(1,hist.GetNbinsX()): #Don't look at 0jet bins
                #Check category exists. This is probably redundant as if it doesn't, it should give an error when calculating the ratio
                category_njet = 'C_{0}_{1}j'.format(category,bin)
                if category_njet not in categories: categories.append(category_njet)

                #Calculate ratio to nominal
                #If the nominal yield is zero, set the ratio to 1.0
                #If the systematic yield is 0, set the ratio to 0.0001 (Combine doesn't like 0)
                bin_ratio = 1.0
                if nom_dict[(process,category_njet)] != 0:
                    bin_yield = round(hist.GetBinContent(1+bin,ROOT.WCPoint()),4)
                    bin_ratio = bin_yield/nom_dict[(process,category_njet)]
                    bin_ratio = max(bin_ratio,0.0001)

                #Create sys_dict key if it doesn't exit; can't edit a dict object that doesn't exist yet
                if not sys_dict.has_key((process,category_njet)):
                    sys_dict[(process,category_njet)] = {}
                sys_dict[(process,category_njet)][systematic] = bin_ratio

    # Only analyze categories with at least a few background events to prevent negative yields
    print "Getting final categories..."
    categories_nonzero = []
    for cat in categories:
        bkgd = 0.
        for proc in bkgd_names:
            bkgd += nom_dict[proc,cat]
        if bkgd >= 0.0001:
            #if cat != "C_3l_ppp_1b_8j":
            categories_nonzero.append(cat)
        else:
            print "Skipping",cat,"for low bkgd yield."
            #print "Background less than 1!"
    #print categories_nonzero

    # Make a plot of relative signal yield shifts per category
    if(0):
        c1 = ROOT.TCanvas('c1','Hist Canvas')
        pad1 = ROOT.TPad()
        hist = ROOT.TH1D('histogram','Yield Shift',40,0,40)
        for idx,entry in enumerate(plot_values[2]): hist.SetBinContent(1+idx,entry/plot_sum[plot_values[1][idx]])
        #for idx,entry in enumerate(plot_values[2]): hist.SetBinContent(1+idx,entry)
        hist.SetLineWidth(2)
        hist.Draw()
        hist.SetStats(0)
        c1.Update()
        c1.Print("YieldPlot_"+operator+".pdf",'.pdf')
        print "Made plot 'YieldPlot_"+operator+".pdf"
        print plot_values

    return(categories_nonzero, data_names, data_dict, sgnl_names, bkgd_names, nom_dict, sys_types, sys_dict)
