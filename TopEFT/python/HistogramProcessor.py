import logging
import math
import ROOT
ROOT.gSystem.Load('$CMSSW_BASE/src/CombineHarvester/TopEFT/interface/TH1EFT_h.so')

class HistogramProcessor(object):
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.sgnl_known = ['ttH','tllq','ttll','ttlnu']
        self.sgnl_histnames = [sgnl + '_' + '16D' for sgnl in self.sgnl_known]
        self.bkgd_known = ['charge_flips','fakes','WZ','ZZ','WW','WWW','WWZ','WZZ','ZZZ','singlet_tWchan','singletbar_tWchan','convs']
        self.data_known = ['data_doubleEle','data_muonEle','data_doubleMu','data_singleEle','data_singleMu']

        # Initialize reweight point for fake data
        # Default to all operators = 1
        WCPoint_string = 'fakedata'
        self.operators_fakedata = [
            'ctW','ctp','cpQM','ctZ','ctG','cbW','cpQ3','cptb',
            'cpt','cQl3i','cQlMi','cQei','ctli','ctei','ctlSi','ctlTi'
        ]
        for op in self.operators_fakedata:
            WCPoint_string += '_{op}_1'.format(op=op)
            #if op == 'ctZ': WCPoint_string += '_{op}_1'.format(op=op)
            #if op == 'ctW': WCPoint_string += '_{op}_0'.format(op=op)
        self.rwgt_pt = ROOT.WCPoint(WCPoint_string,1.0)
        self.sm_pt = ROOT.WCPoint()

    def process(self,infile,fake_data):
        self.logger.info("Setting up...")
        readfile = ROOT.TFile.Open(infile)

        #Lists of names to pass on
        categories=[] #e.g. 2los_ee_2j_1b
        data_names=[] #e.g. doubleEle
        sgnl_names=[] #e.g. ttH
        bkgd_names=[] #e.g. ttH
        sys_types=[] #e.g. Lumi

        #Dicts for rates
        data_dict={} #process,category:rate
        fakedata_dict={} #process,category:rate
        nom_dict={} #process,category:nominal rate
        sys_dict={} #process,category:{systematic type:ratio to nominal rate}

        debug_processes = []    # Ex: 'ttll_cpt','ttZ'
        debug_categories = []   # Ex: 'C_2lss_m_emu_2b_5j'

        if fake_data:
            self.logger.info("Using operators %s for fake data.",self.operators_fakedata)
        else:
            self.logger.info("Not using fake data. #CURRENTLY SM PREDICTION, not real data#")

        self.logger.info("Looping through histograms...")

        #Main parsing loop
        for key in readfile.GetListOfKeys():
            hist = readfile.Get(key.GetName())

            #Get categorical information from histogram name
            histname = hist.GetName().split('.')
            category,systematic,process = '','',''
            if(len(histname)==3): [category,systematic,process] = histname
            if(len(histname)==2): [category,process] = histname

            if process in debug_processes:
                for bin in range(1,hist.GetNbinsX()):
                    category_njet = 'C_{0}_{1}j'.format(category,bin)
                    #bin_yield = round(hist.GetBinContent(1+bin,ROOT.WCPoint()),4)
                    bin_yield = round(hist.GetBinContent(1+bin,self.sm_pt),4)
                    self.logger.debug("%s %s %s",process,category_njet,str(bin_yield))

            #Logic for data NOT VERIFIED
            if process in self.data_known:
                if process not in data_names: data_names.append(process.replace('data_',''))
                for bin in range(1,hist.GetNbinsX()):
                    category_njet = 'C_{0}_{1}j'.format(category,bin)
                    if category_njet not in categories: categories.append(category_njet)
                    #bin_yield = round(hist.GetBinContent(1+bin,ROOT.WCPoint()),4)
                    bin_yield = round(hist.GetBinContent(1+bin,self.sm_pt),4)
                    data_dict.update({(process,category_njet):bin_yield})

            # Logic for MC yields below
            if process not in self.bkgd_known+self.sgnl_histnames: continue

            # Obtain signal yields, being sure to use correct histograms
            if process in self.sgnl_histnames:
                process = process.rsplit("_",1)[0]
                if process not in sgnl_names: sgnl_names.append(process)
                #process = process.replace('tllq','tZq')
                #process = process.replace('ttll','ttZ')
                #process = process.replace('ttlnu','ttW')
            if process in self.bkgd_known:
                if process not in bkgd_names: bkgd_names.append(process)

            # Logic for the nominal histograms
            if systematic == '':
                self.logger.debug("Nominal Hist: %s",hist.GetName())
                for bin in range(1,hist.GetNbinsX()):
                    # MC Nominal
                    category_njet = 'C_{0}_{1}j'.format(category,bin)
                    if category_njet not in categories: categories.append(category_njet)
                    bin_yield = round(hist.GetBinContent(1+bin,self.sm_pt),4)
                    nom_dict.update({(process,category_njet):bin_yield})
                    # Fake data
                    if fake_data:
                        fakedata_bin_yield = round(hist.GetBinContent(1+bin,self.rwgt_pt),4)
                        fakedata_dict.update({(process,category_njet):fakedata_bin_yield})

                # Get MCStats uncertainty for the nominal histograms
                for bin in range(1,hist.GetNbinsX()): #Don't look at 0jet bins
                    #Check category exists. This is probably redundant as if it doesn't, it should give an error when calculating the ratio
                    category_njet = 'C_{0}_{1}j'.format(category,bin)
                    if category_njet not in categories: categories.append(category_njet)

                    #Calculate ratio to nominal
                    #If the nominal yield is zero, set the ratio to 1.0
                    #If the systematic yield is 0, set the ratio to 0.0001 (Combine doesn't like 0)
                    bin_ratio = 1.0
                    if nom_dict[(process,category_njet)] != 0:
                        if process in self.sgnl_known: # Determined by TH1EFT fit
                            bin_yield = round(hist.GetBinFit(1+bin).evalPointError(self.sm_pt),4)
                        if process in self.bkgd_known: # Determined by sqrt(yield)
                            bin_yield = round(math.sqrt(max(0,hist.GetBinContent(1+bin,self.sm_pt))),4)
                        bin_ratio = 1+bin_yield/nom_dict[(process,category_njet)]
                        bin_ratio = max(bin_ratio,0.0001)

                    #Create sys_dict key if it doesn't exit; can't edit a dict object that doesn't exist yet
                    if not sys_dict.has_key((process,category_njet)):
                        sys_dict[(process,category_njet)] = {}
                    sys_dict[(process,category_njet)]['MCSTATS'] = bin_ratio

            # Logic for systematic histograms
            else:
                self.logger.debug("Systematic Hist: %s",hist.GetName())
                if systematic not in sys_types: sys_types.append(systematic)

                for bin in range(1,hist.GetNbinsX()): #Don't look at 0jet bins
                    # Check category exists. This is probably redundant as if it doesn't, it should 
                    #   give an error when calculating the ratio
                    category_njet = 'C_{0}_{1}j'.format(category,bin)
                    if category_njet not in categories: categories.append(category_njet)

                    #Calculate ratio to nominal
                    #If the nominal yield is zero, set the ratio to 1.0
                    #If the systematic yield is 0, set the ratio to 0.0001 (Combine doesn't like 0)
                    bin_ratio = 1.0
                    if nom_dict[(process,category_njet)] != 0:
                        bin_yield = round(hist.GetBinContent(1+bin,self.sm_pt),4)
                        bin_ratio = bin_yield/nom_dict[(process,category_njet)]
                        bin_ratio = max(bin_ratio,0.0001)

                        #Special case for fake rate uncertainty; average effect over all njets bins
                        if systematic in ['FRUP','FRDOWN']:
                            bin_yield = hist.Integral(1,hist.GetNbinsX())
                            bin_ratio = round(bin_yield/readfile.Get(category+'.'+process).Integral(1,hist.GetNbinsX()),4)
                            bin_ratio = max(bin_ratio,0.0001)

                    #Create sys_dict key if it doesn't exit; can't edit a dict object that doesn't exist yet
                    if not sys_dict.has_key((process,category_njet)):
                        sys_dict[(process,category_njet)] = {}
                    sys_dict[(process,category_njet)][systematic] = bin_ratio

                #Append systematic type to list, but don't split into UP/DOWN
                sys_type = systematic
                if systematic.endswith('UP'): sys_type = systematic[:-2]
                elif systematic.endswith('DOWN'): sys_type = systematic[:-4]
                if sys_type not in sys_types: sys_types.append(sys_type)

        # Only analyze categories with at least a few background events to prevent negative yields
        self.logger.info("Getting final categories...")
        categories_nonzero = []
        for cat in categories:
            bkgd = 0.
            for proc in bkgd_names:
                if cat in debug_categories: self.logger.debug("%s %s",proc,str(nom_dict[proc,cat]))
                bkgd += nom_dict[proc,cat]
            sgnl = 0.
            for proc in sgnl_names:
                sgnl += nom_dict[proc,cat]
                if cat in debug_categories: self.logger.debug("%s %s",proc,str(nom_dict[proc,cat]))
            if sgnl >= 0.01:
                if bkgd >= 0.01:
                    categories_nonzero.append(cat)
                else:
                    self.logger.info("Skipping %s for low background yield.",cat)
            else:
                self.logger.info("Skipping %s for low signal yield.",cat)
            if cat in debug_categories: self.logger.debug("%s %s %s",str(bkgd),str(sgnl),str(bkgd+sgnl))

        #categories_nonzero = ['C_2lss_m_emu_2b_5j'] # Override for debug
        self.logger.info("Categories: %s",str(categories_nonzero))
        self.logger.info("Signals: %s",str(sgnl_names))
        self.logger.info("Backgrounds: %s",str(bkgd_names))
        #print data_dict
        #print nom_dict

        return(categories_nonzero, data_names, data_dict, fakedata_dict, sgnl_names, bkgd_names, nom_dict, sys_types, sys_dict)

    ##############################
    # Getter/Setter methods
    ##############################

    def setOperators(self,lst):
        self.operators_fakedata = lst

    def setReweightPoint(self,d):
        self.rwgt_pt.setSMPoint()   # Reset all WC strengths to 0
        for wc,v in d.iteritems():
            if not wc in self.operators_fakedata:
                continue
            self.rwgt_pt.setStrength(wc,v)

    def setSignalProcesses(self,lst):
        self.sgnl_known = lst

    def setBackgroundProcesses(self,lst):
        self.bkgd_known = lst

    def setDatasets(self,lst):
        self.data_known = lst

    def getOperators(self):
        return self.operators_fakedata

    def getReweightPoint(self):
        return self.rwgt_pt

    def getSignalProcesses(self):
        return self.sgnl_known

    def getBackgroundProcesses(self):
        return self.bkgd_known

    def getDatasets(self):
        return self.data_known
