import sys
import logging
import CombineHarvester.CombineTools.ch as ch
from HistogramProcessor import HistogramProcessor

class DatacardMaker(object):
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.hp = HistogramProcessor()

        # Initialize CombineHarvester instance
        self.cb = ch.CombineHarvester()
        self.cb.SetVerbosity(0)

        self.debug = 0
        self.eras = ['2017']    #Data eras
        self.chan = ['']        #Indistiguishable process subcategories i.e. for MC only

        self.outf = "EFT_MultiDim_Datacard.txt"

    def make(self,infile,fake_data=False):
        self.logger.info("Parsing input file!")
        (categories, data_names, data_dict, fakedata_dict, sgnl_names, bkgd_names, nom_dict, sys_types, sys_dict) = self.hp.process(infile,fake_data)
        self.logger.info("Done parsing input file")
        self.logger.info("Now creating Datacard!")

        cats = list(enumerate(categories))  #Process bins. Must be list of tuple like this

        #Fill observation
        obs_rates={}
        #Fake data
        if fake_data:
            self.logger.info("Using fake data")
            for cat in categories:
                cat_asimov = 0
                for proc in sgnl_names+bkgd_names:
                    cat_asimov += fakedata_dict[proc,cat]
                obs_rates[cat]=cat_asimov
        #Asimov data
        else:
            for cat in categories:
                cat_asimov = 0
                for proc in sgnl_names+bkgd_names:
                    cat_asimov += nom_dict[proc,cat]
                obs_rates[cat]=cat_asimov

        #Actual data (hists not currently filled, but placeholders are present)
        #obs_rates={}
        #for cat in categories:
        #    obs_rates[cat]=0
        #for (proc, cat) in data_dict.keys():
        #   if cat in categories: obs_rates[cat] += data_dict[proc,cat]

        #Initialize structure of each observation (data) and process (MC signal and background)
        #Can be adjusted very heavily to specialize each era/process/category
        # ['*'] is Higgs mass hypothesis, which we ignore
        self.cb.AddObservations( ['*'], ['Top_EFT'], self.eras, self.chan, cats )

        self.cb.AddProcesses( ['*'], ['Top_EFT'], self.eras, self.chan, sgnl_names, cats, True) 
        self.cb.AddProcesses( ['*'], ['Top_EFT'], self.eras, self.chan, bkgd_names, cats, False)

        #Fill the data rates
        self.logger.info("Adding observation rates...")
        self.cb.ForEachObs(lambda x: x.set_rate(obs_rates[x.bin()])) #Not split by process/channel, obviously

        #Function for making sure bin yield > 0
        def checkRate(x):
            if nom_dict[x.process(),x.bin()]>0:
                x.set_rate(nom_dict[x.process(),x.bin()])

        #Fill the nominal MC rates
        self.logger.info("Adding MC rates...")
        #self.cb.ForEachProc(lambda x: x.set_rate(nom_dict[x.process(),x.bin()]))
        self.cb.ForEachProc(lambda x: checkRate(x))

        #Round systematics (Only for debug purposes when viewing datacard! Keep full accuracy otherwise!)
        #for outerkey,outervalue in sys_dict.items():
        #    for innerkey,innervalue in outervalue.items():
        #        sys_dict[outerkey][innerkey]=round(innervalue,4)

        #Fill systematic rates
        for proc in sgnl_names+bkgd_names:
            self.logger.info("Adding systematics for %s...",proc)

            #PDF/Q2 rate uncertainty values
            PDFrate = 1.0
            Q2rate = [1.0,1.0]
            if proc=='ttH': # Done
                PDFrate = 1.036
                Q2rate  = [0.908,1.058]
            if proc=='ttlnu': # Done
                PDFrate = 1.02
                Q2rate  = [0.88,1.13]
            if proc=='ttll': # Done
                PDFrate = 1.03
                Q2rate  = [0.88,1.10]
            if proc=='tllq': # V+jets?? Uncertain!
                PDFrate = 1.04
                Q2rate  = [0.99,1.01]
            if proc=='tHq': # V+jets?? Uncertain!
                PDFrate = 1.037
                Q2rate  = [0.92,1.06]
            if proc in ['singlet_tWchan','singletbar_tWchan']: # Done
                PDFrate = 1.03
                Q2rate  = [0.98,1.03]
            if proc in ['Diboson']: # Done
                PDFrate = 1.02
                Q2rate  = [0.98,1.02]
            if proc in ['Triboson']: # See ATLAS paper https://arxiv.org/pdf/1610.05088.pdf
                PDFrate = 1.042
                Q2rate  = [0.974,1.026]
            if proc in ['charge_flips','fakes']: # Data-driven, so none (These won't go in the datacard)
                PDFrate = 1.0
                Q2rate  = [1.0,1.0]
            if proc in ['ttGJets']: # Unknown; conservative here
                PDFrate = 1.5
                Q2rate = [0.90,1.10]

            for cat in categories:
                if '2lss' in cat:
                    if '4j' in cat:
                        PSISRDOWN = 0.95
                        PSISRUP = 1.05
                    if '5j' in cat:
                        PSISRDOWN = 0.975
                        PSISRUP = 1.025
                    if '6j' in cat:
                        PSISRDOWN = 1.025
                        PSISRUP = 0.975
                    if '7j' in cat:
                        PSISRDOWN = 1.05
                        PSISRUP = 0.95
                if '3l' in cat:
                    if '2j' in cat:
                        PSISRDOWN = 0.95
                        PSISRUP = 1.05
                    if '3j' in cat:
                        PSISRDOWN = 0.975
                        PSISRUP = 1.025
                    if '4j' in cat:
                        PSISRDOWN = 1.025
                        PSISRUP = 0.975
                    if '5j' in cat:
                        PSISRDOWN = 1.05
                        PSISRUP = 0.95
                if '4l' in cat:
                    if '1j' in cat:
                        PSISRDOWN = 0.95
                        PSISRUP = 1.05
                    if '2j' in cat:
                        PSISRDOWN = 0.975
                        PSISRUP = 1.025
                    if '3j' in cat:
                        PSISRDOWN = 1.025
                        PSISRUP = 0.975
                    if '4j' in cat:
                        PSISRDOWN = 1.05
                        PSISRUP = 0.95
                if nom_dict[proc,cat]:
                    #MCStats uncertainty (fully correlated, taken from nominal bin errors)
                    #self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'MCStats','lnN',ch.SystMap()( sys_dict[(proc,cat)]['MCSTATS']))
                    if proc=='fakes': self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'FR_stats','lnN',ch.SystMap()( sys_dict[(proc,cat)]['MCSTATS']))
                    #Lumi uncertainty (fully correlated, flat rate, identical for all categories)
                    self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'lumi_13TeV_2017','lnN',ch.SystMap()( 1.023 ))
                    #Charge Flip rate uncertainty (fully correlated, flat rate, identical for all categories, Charge Flip process only)
                    if proc=='charge_flips': self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'ChargeFlips','lnN',ch.SystMap()( 1.30 ))
                    #PDF rate uncertainty (correlated within process, flat rate, identical for all categories within process)
                    if proc in ['ttH']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'pdf_ggttH','lnN',ch.SystMap()( 1/PDFrate ))
                    if proc in ['ttll','ttGJets']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'pdf_gg','lnN',ch.SystMap()( 1/PDFrate ))
                    if proc in ['ttlnu','tllq','Diboson','Triboson']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'pdf_qq','lnN',ch.SystMap()( 1/PDFrate ))
                    if proc in ['tHq']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'pdf_qgtHq','lnN',ch.SystMap()( 1/PDFrate ))
                    #Q2 rate uncertainty (correlated within process, flat rate, identical for all categories within process)
                    if proc in ['ttH']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'QCDscale_ttH','lnN',ch.SystMap()( Q2rate ))
                    if proc in ['ttll','ttlnu']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'QCDscale_ttbar','lnN',ch.SystMap()( Q2rate ))
                    if proc in ['tllq']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'QCDscale_V','lnN',ch.SystMap()( Q2rate ))
                    if proc in ['Diboson']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'QCDscale_VV','lnN',ch.SystMap()( Q2rate ))
                    if proc in ['Triboson']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'QCDscale_VVV','lnN',ch.SystMap()( Q2rate ))
                    if proc in ['ttGJets']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'QCDscale_ttG','lnN',ch.SystMap()( Q2rate ))
                    #PSISR for anatest14 ONLY
                    if proc not in ['fakes','charge_flips']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'PSISR','lnN',ch.SystMap()( [PSISRDOWN,PSISRUP] ))
                    #Standard uncertainties with usual UP/DOWN variations
                    #Includes FR, JES, CERR1, CERR2, HF, HFSTATS1, HFSTATS2, LF, LFSTATS1, LFSTATS2, MUR, MUF, LEPID, TRG, PU, PSISR
                    for sys in sys_types:
                        # Use CMS-standard names for uncertainties
                        sys_name = sys
                        if sys == 'JES': sys_name = 'CMS_scale_j'
                        if sys == 'TRG': sys_name = 'CMS_eff_em'
                        if sys in ['HF','HFSTATS1','HFSTATS2','LF','LFSTATS1','LFSTATS2']:
                            sys_name = sys.lower()

                        if sys+'UP' not in sys_dict[(proc,cat)].keys(): continue
                        self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,sys_name,'lnN',ch.SystMap()( [sys_dict[(proc,cat)][sys+'DOWN'], sys_dict[(proc,cat)][sys+'UP']] ))

        #Printout of signal and background yields (debug)
        if self.debug:
            background = {}
            signal = {}
            print ""
            self.logger.debug("\nCategory, signal yield, background yield:")
            for cat in categories:
                background[cat] = 0
                signal[cat] = 0
                for proc in bkgd_names:
                    background[cat] += nom_dict[proc,cat]
                for proc in sgnl_names:
                    signal[cat] += nom_dict[proc,cat]
                print cat,signal[cat],background[cat]

        #self.cb.PrintAll() #Print the datacard
        self.logger.info("Writing datacard '%s'...",self.outf)
        self.cb.WriteDatacard(self.outf)

    ##############################
    # Getter/Setter methods (mostly just pass-throughs)
    ##############################

    def setOperators(self,lst):
        self.hp.setOperators(lst)

    def setReweightPoint(self,d):
        self.hp.setReweightPoint(d)

    def setSignalProcesses(self,lst):
        self.hp.setSignalProcesses(lst)

    def setBackgroundProcesses(self,lst):
        self.hp.setBackgroundProcesses(lst)

    def setDatasets(self,lst):
        self.hp.setDatasets(lst)

    def getOperators(self):
        return self.hp.getOperators()

    def getReweightPoint(self):
        return self.hp.getReweightPoint()

    def getSignalProcesses(self):
        return self.hp.getSignalProcesses()

    def getBackgroundProcesses(self):
        return self.hp.getBackgroundProcesses()

    def getDatasets(self):
        return self.hp.getDatasets()

if __name__ == "__main__":
    log_file = 'DatacardMaker.log'

    FORMAT1 = '%(message)s'
    FORMAT2 = '[%(levelname)s] %(message)s'
    FORMAT3 = '[%(levelname)s][%(name)s] %(message)s'

    frmt1 = logging.Formatter(FORMAT1)
    frmt2 = logging.Formatter(FORMAT2)
    frmt3 = logging.Formatter(FORMAT3)

    logging.basicConfig(
        level=logging.DEBUG,
        format=FORMAT2,
        filename=log_file,
        filemode='w'
    )

    # Configure logging to also output to stdout
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(frmt2)
    logging.getLogger('').addHandler(console)

    # Check for coefficient argument
    fake_data = False
    if len(sys.argv) == 2:
        if sys.argv[1] in ['True', 'true', '1']:
            fake_data = True
        elif sys.argv[1] in ['False', 'false', '0']:
            fake_data = False
        else:
            logging.error("Value of argument 1 unrecognized!")
            sys.exit()
    if len(sys.argv) > 2:
        logging.error("Only one argument allowed!")
        sys.exit()

    # Run datacard maker
    dm = DatacardMaker()
    dm.make('../hist_files/anatest14.root',fake_data)

    logging.info("Logger shutting down!")
    logging.shutdown()



