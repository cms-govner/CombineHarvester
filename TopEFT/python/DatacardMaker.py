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

        #Fill the nominal MC rates
        self.logger.info("Adding MC rates...")
        self.cb.ForEachProc(lambda x: x.set_rate(nom_dict[x.process(),x.bin()]))

        #Fill systematic rates
        for proc in sgnl_names+bkgd_names:
            self.logger.info("Adding systematics for %s...",proc)

            #PDF/Q2 rate uncertainty values
            PDFrate = 1.0
            Q2rate = [1.0,1.0]
            if proc=='ttH':
                PDFrate = 1.036
                Q2rate  = [0.908,1.058]
            if proc=='ttlnu':
                PDFrate = 1.02
                Q2rate  = [0.88,1.13]
            if proc=='ttll':
                PDFrate = 1.03
                Q2rate  = [0.88,1.10]
            if proc in ['singlet_tWchan','singletbar_tWchan']:
                PDFrate = 1.03
                Q2rate  = [0.98,1.03]
            if proc=='tllq': # V+jets??
                PDFrate = 1.04
                Q2rate  = [0.99,1.01]
            if proc in ['WZ','ZZ','WW']:
                PDFrate = 1.02
                Q2rate  = [0.98,1.02]
            if proc in ['charge_flips','fakes']:
                PDFrate = 1.0
                Q2rate  = [1.0,1.0]

            for cat in categories:
                #MCStats uncertainty (fully correlated, taken from nominal bin errors)
                self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'MCStats','lnN',ch.SystMap()( sys_dict[(proc,cat)]['MCSTATS']))
                #Lumi uncertainty (fully correlated, flat rate, identical for all categories)
                self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'Lumi','lnN',ch.SystMap()( 1.025 ))
                #Charge Flip rate uncertainty (fully correlated, flat rate, identical for all categories, Charge Flip process only)
                if proc=='charge_flips': self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'ChargeFlips','lnN',ch.SystMap()( 1.30 ))
                #PDF rate uncertainty (correlated within process, flat rate, identical for all categories within process)
                if proc=='ttH': self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'PDF_ggttH','lnN',ch.SystMap()( PDFrate ))
                if proc in ['ttll']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'PDF_gg','lnN',ch.SystMap()( PDFrate ))
                if proc in ['ttlnu','tllq','WZ','ZZ','WW']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'PDF_qq','lnN',ch.SystMap()( PDFrate ))
                if proc in ['singlet_tWchan','singletbar_tWchan']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'PDF_qg','lnN',ch.SystMap()( PDFrate ))
                #Q2 rate uncertainty (correlated within process, flat rate, identical for all categories within process)
                if proc=='ttH': self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'Q2_ttH','lnN',ch.SystMap()( Q2rate ))
                if proc in ['ttll','ttlnu']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'Q2_tt','lnN',ch.SystMap()( Q2rate ))
                if proc in ['singlet_tWchan','singletbar_tWchan']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'Q2_t','lnN',ch.SystMap()( Q2rate ))
                if proc=='tllq': self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'Q2_V','lnN',ch.SystMap()( Q2rate ))
                if proc in ['WZ','ZZ','WW']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'Q2_VV','lnN',ch.SystMap()( Q2rate ))
                #Standard uncertainties with usual UP/DOWN variations
                #Includes FR, JES, CERR1, CERR2, HF, HFSTATS1, HFSTATS2, LF, LFSTATS1, LFSTATS2, MUR, MUF, LEPID
                for sys in sys_types:
                    if sys+'UP' not in sys_dict[(proc,cat)].keys(): continue
                    self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,sys,'lnN',ch.SystMap()( [sys_dict[(proc,cat)][sys+'UP'], sys_dict[(proc,cat)][sys+'DOWN']]))

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

    ################################################################################################

    def setOperators(self,lst):
        self.hp.operators_known = lst

    def setReweightPoint(self,d):
        self.hp.setReweightPoint(d)

    def setSignalProcesses(self,lst):
        self.hp.sgnl_known = lst

    def setBackgroundProcesses(self,lst):
        self.bkgd_known = lst


if __name__ == "__main__":
    log_file = 'out.log'

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
    dm.make('../hist_files/anatest10.root',fake_data)

    logging.info("Logger shutting down!")
    logging.shutdown()



