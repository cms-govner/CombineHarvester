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
        self.eras = ['2016']    #Data eras
        self.chan = ['']        #Indistiguishable process subcategories i.e. for MC only

        self.outf = "EFT_MultiDim_Datacard.txt"

    def make(self,infile,fake_data=True):
        self.logger.info("Parsing input file!")
        (categories, data_names, data_dict, fakedata_dict, sgnl_names, bkgd_names, nom_dict, sys_types, sys_dict) = self.hp.process(infile,fake_data)
        self.logger.info("Done parsing input file")
        self.logger.info("Now creating Datacard!")

        cats = list(enumerate(categories))  #Process bins. Must be list of tuple like this

        #Fill observation
        obs_rates={}
        #Fake data
        if fake_data:
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
            for cat in categories:
                #Lumi uncertainty (fully correlated, identical for all categories)
                self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'Lumi','lnN',ch.SystMap()( 1.025 ))
                #MCStats uncertainty (fully correlated)
                self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'MCStats','lnN',ch.SystMap()( sys_dict[(proc,cat)]['MCSTATS']))
                #PDF uncertainty (fully correlated, identical for all signal, only signal)
                if proc in sgnl_names: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'PDF','lnN',ch.SystMap()( 1.015 ))
                #Charge Flip rate uncertainty (fully correlated, identical for all categories, Charge Flip only)
                if proc=='charge_flips': self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'ChargeFlips','lnN',ch.SystMap()( 1.30 ))
                #Fake rate uncertainty (fully correlated, asymmetric, Fakes only)
                if proc=='fakes': self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'Fakes','lnN',ch.SystMap()( [sys_dict[(proc,cat)]['FRUP'], sys_dict[(proc,cat)]['FRDOWN']] ))
                #JES uncertainty (fully correlated, asymmetric)
                self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'JES','lnN',ch.SystMap()( [sys_dict[(proc,cat)]['JESUP'], sys_dict[(proc,cat)]['JESDOWN']]))
                if (proc,cat) in sys_dict.keys(): # probably unnecessary safeguard
                    #Other systematics -- To be added
                    for sys_type in sys_types:
                        #Fully correlated systematics listed by rate
                        if sys_type in []:
                            self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,sys_type,'lnN',ch.SystMap()( sys_dict[(proc,cat)][sys_type] ))
                        #Systematics listed by their fluctuations
                        elif sys_type in ['pdfUP']:
                            unc = (nom_dict[proc,cat]+sys_dict[proc,cat][sys_type])/nom_dict[proc,cat] if sys_dict[proc,cat][sys_type]>0. else 1.0000
                            self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,sys_type.rstrip('UP'),'lnN',ch.SystMap()(unc))
                        #Shape systematics (not used in counting experiment)
                        #elif sys_type in ['MCStatUP','MCStatDOWN','Q2UP','Q2DOWN']:
                        #Superfluous systematics
                        #Fully uncorrelated systematics
                        #else:
                            #self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,proc+cat+':'+sys_type,'lnN',ch.SystMap()( float(sys_dict[(proc,cat)][sys_type]) ))

        #Printout of signal and background yields (debug)
        if self.debug:
            background = {}
            signal = {}
            print ""
            self.logger.debug("Category, signal yield, background yield:")
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
    fake_data = True
    if len(sys.argv) == 2:
        fake_data = sys.argv[1]
    dm = DatacardMaker()
    dm.make('../data/anatest9.root',fake_data)

    logging.info("Logger shutting down!")
    logging.shutdown()



