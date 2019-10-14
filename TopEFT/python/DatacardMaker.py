import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
import argparse
import CombineHarvester.CombineTools.ch as ch
from HistogramProcessor import HistogramProcessor

class DatacardMaker(object):
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.hp = HistogramProcessor()

        # Initialize CombineHarvester instance
        self.cb = ch.CombineHarvester()
        self.cb.SetVerbosity(0)

        self.debug = 1
        self.yieldTeX = 0
        
        self.eras = ['2017']    #Data eras
        self.chan = ['']        #Indistiguishable process subcategories i.e. for MC only

        self.outf = "EFT_MultiDim_Datacard.txt"
        
        self.minyield = 0.0000 # Minimum nominal process+category yield required

    def make(self,infile,fake_data=False,central=False):
        self.logger.info("Parsing input file!")
        (categories, data_names, data_dict, fakedata_dict, sgnl_names, bkgd_names, nom_dict, sys_types, sys_dict) = self.hp.process(infile,fake_data,central)
        self.logger.info("Done parsing input file")
        self.logger.info("Now creating Datacard!")

        cats = list(enumerate(categories))  #Process bins. Must be list of tuple like this

        #Fill observation
        obs_rates={}
        #Fake data
        if fake_data:
            self.logger.info("Using fake data")
            for cat in categories:
                cat_fakeyield = 0
                for proc in sgnl_names+bkgd_names:
                    if nom_dict[proc,cat] > self.minyield: # Only use this bin if the process has meaningful yield
                        cat_fakeyield += fakedata_dict[proc,cat]
                obs_rates[cat]=cat_fakeyield
        #Actual data
        else:
            obs_rates={}
            for cat in categories:
                obs_rates[cat]=0
            for (proc, cat) in data_dict.keys():
                if cat in categories: obs_rates[cat] += data_dict[proc,cat]

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
            if nom_dict[x.process(),x.bin()]>self.minyield:
                x.set_rate(nom_dict[x.process(),x.bin()])

        #Fill the nominal MC rates
        self.logger.info("Adding MC rates...")
        self.cb.ForEachProc(lambda x: checkRate(x))

        #Round systematics (Only for debug purposes when viewing datacard! Keep full accuracy otherwise!)
        #for outerkey,outervalue in sys_dict.items():
        #    for innerkey,innervalue in outervalue.items():
        #        sys_dict[outerkey][innerkey]=round(innervalue,2)

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
            if proc in ['convs']: # Unknown; conservative here
                PDFrate = 1.05
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

                if nom_dict[proc,cat]>self.minyield: # Only add systematic if yield is nonzero.
                    #MCStats uncertainty, currently unused (fully correlated, taken from nominal bin errors)
                    #self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'MCStats','lnN',ch.SystMap()( sys_dict[(proc,cat)]['MCSTATS']))
                    #FR_stats uncertainty (fully uncorrelated, taken from MC stats error)
                    if proc=='fakes': self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'FR_stats'+cat.strip('C'),'lnN',ch.SystMap()( sys_dict[(proc,cat)]['MCSTATS']))
                    #Lumi uncertainty (fully correlated, flat rate, identical for all categories)
                    self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'lumi_13TeV_2017','lnN',ch.SystMap()( 1.023 ))
                    #Charge Flip rate uncertainty (fully correlated, flat rate, identical for all categories, Charge Flip process only)
                    if proc=='charge_flips': self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'ChargeFlips','lnN',ch.SystMap()( 1.30 ))
                    #PDF rate uncertainty (correlated within parent process, flat rate, identical for all categories within process)
                    if proc in ['ttH']:   self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'pdf_ggttH','lnN',ch.SystMap()( 1/PDFrate ))
                    if proc in ['ttll']:  self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'pdf_gg','lnN',ch.SystMap()( 1/PDFrate ))
                    if proc in ['ttlnu']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'pdf_qq','lnN',ch.SystMap()( 1/PDFrate ))
                    if proc in ['tllq']:  self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'pdf_qq','lnN',ch.SystMap()( 1/PDFrate ))
                    if proc in ['tHq']:   self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'pdf_qgtHq','lnN',ch.SystMap()( 1/PDFrate ))
                    if proc in ['Diboson','Triboson']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'pdf_qq','lnN',ch.SystMap()( 1/PDFrate ))
                    if proc in ['convs']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'pdf_gg','lnN',ch.SystMap()( 1/PDFrate ))
                    #Q2 rate uncertainty (correlated within parent process, flat rate, identical for all categories within process)
                    if proc in ['ttH']:   self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'QCDscale_ttH','lnN',ch.SystMap()( Q2rate ))
                    if proc in ['ttll']:  self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'QCDscale_ttbar','lnN',ch.SystMap()( Q2rate ))
                    if proc in ['ttlnu']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'QCDscale_ttbar','lnN',ch.SystMap()( Q2rate ))
                    if proc in ['tllq']:  self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'QCDscale_V','lnN',ch.SystMap()( Q2rate ))
                    if proc in ['tHq']:   self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'QCDscale_tHq','lnN',ch.SystMap()( Q2rate ))
                    if proc in ['Diboson']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'QCDscale_VV','lnN',ch.SystMap()( Q2rate ))
                    if proc in ['Triboson']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'QCDscale_VVV','lnN',ch.SystMap()( Q2rate ))
                    if proc in ['convs']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'QCDscale_ttG','lnN',ch.SystMap()( Q2rate ))
                    #PSISR for anatest14&15. After that, it's included as up/down systematic. The hist file should overwrite this line anyway.
                    if proc not in ['fakes','charge_flips']: self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,'PSISR','lnN',ch.SystMap()( [PSISRDOWN,PSISRUP] ))
                    #Standard uncertainties with usual UP/DOWN variations
                    #Includes FR_shape, JES, CERR1, CERR2, HF, HFSTATS1, HFSTATS2, LF, LFSTATS1, LFSTATS2, MUR, MUF, LEPID, TRG, PU, PSISR
                    for sys in sys_types:
                        # Use CMS-standard names for uncertainties
                        sys_name = sys
                        if sys == 'JES': sys_name = 'CMS_scale_j'
                        if sys == 'TRG': sys_name = 'CMS_eff_em'
                        if sys in ['HF','HFSTATS1','HFSTATS2','LF','LFSTATS1','LFSTATS2']:
                            sys_name = sys.lower()

                        if sys+'UP' not in sys_dict[(proc,cat)].keys(): continue
                        #if sys in ['CERR1','CERR2']: # DEBUG for CERR. Might not want this!
                        #    sysup = sys_dict[(proc,cat)][sys+'UP']
                        #    sysdown = sys_dict[(proc,cat)][sys+'DOWN']
                        #    sysavg = 1+(abs(sysup-1)+abs(sysdown-1))/2
                        #    self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,sys_name,'lnN',ch.SystMap()( sysavg ))
                        #else:
                        
                        # Largest of errors
                        #if sys in ['JES'] and ((sys_dict[(proc,cat)][sys+'DOWN']-1)*(sys_dict[(proc,cat)][sys+'UP']-1))>0:
                        #if sys in ['JES']:
                        #    uperr = sys_dict[(proc,cat)][sys+'UP']-1
                        #    downerr = sys_dict[(proc,cat)][sys+'DOWN']-1
                        #    symerr = max(abs(uperr),abs(downerr))
                        #    if uperr<0:
                        #        sym_JES = 1-symerr
                        #    else:
                        #        sym_JES = 1+symerr
                        #    #print("Symmetrizing errors: {}, {} for rate {}".format(downerr,uperr,nom_dict[(proc,cat)]))
                        #    self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,sys_name,'lnN',ch.SystMap()( sym_JES ))
                        # Partial Symm, Average.
                        #if sys in ['JES'] and ((sys_dict[(proc,cat)][sys+'DOWN']-1)*(sys_dict[(proc,cat)][sys+'UP']-1))>0:
                        #    sym_JES = (sys_dict[(proc,cat)][sys+'DOWN']+sys_dict[(proc,cat)][sys+'UP'])/2
                        #    self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,sys_name,'lnN',ch.SystMap()( sym_JES ))
                        # Full Symm, Average.
                        if sys in ['JES']:
                            uperr = sys_dict[(proc,cat)][sys+'UP']-1
                            downerr = sys_dict[(proc,cat)][sys+'DOWN']-1
                            symerr = (abs(uperr)+abs(downerr))/2
                            if uperr<0:
                                sym_JES = 1-symerr
                            else:
                                sym_JES = 1+symerr
                            self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,sys_name,'lnN',ch.SystMap()( sym_JES ))
                        else:
                            self.cb.cp().process([proc]).bin([cat]).AddSyst(self.cb,sys_name,'lnN',ch.SystMap()( [sys_dict[(proc,cat)][sys+'DOWN'], sys_dict[(proc,cat)][sys+'UP']] ))

        # Make nuisance groups for easy testing in combine
        self.cb.SetGroup('TheoryNuisances',['^pdf.*','^QCDscale.*'])
        self.cb.SetGroup('bTagNuisances',['^hf.*','^lf.*','^CERR.*'])
        self.cb.SetGroup('DDBkgdNuisances',['FR_shape','ChargeFlips','^FR_stats.*'])
        self.cb.SetGroup('SystematicNuisances',['^.*'])
        
        # Printout of signal and background yields (debug)
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
                print cat,signal[cat]+background[cat],signal[cat],background[cat]
        
        # Printout of yield table in TeX form
        if self.yieldTeX:
            print ""
            print "TeX for yield table:"
            # First sum up the njet categories
            categories_combj = [cat.rsplit('_',1)[0].replace('_','\_') for cat in categories]
            categories_combj = list(set(categories_combj))
            categories_combj.sort()
            combj_yields = {}
            combj_data = {}
            for cat in categories:
                for proc in bkgd_names+sgnl_names:
                    combj_yields[(proc,cat.rsplit('_',1)[0].replace('_','\_'))]=0
                    combj_data[cat.rsplit('_',1)[0].replace('_','\_')]=0
            for cat in categories:
                combj_data[cat.rsplit('_',1)[0].replace('_','\_')] += obs_rates[cat]
                for proc in bkgd_names+sgnl_names:
                    combj_yields[(proc,cat.rsplit('_',1)[0].replace('_','\_'))] += nom_dict[proc,cat]

            # Print category line
            print '& '+' & '.join(categories_combj)+' \\\\'
            print '\hline'
            # Print background yields
            for proc in bkgd_names:
                line = proc
                for cat in categories_combj:
                    line += ' & '+'%.2f'%combj_yields[proc,cat]
                print line+' \\\\'
            print '\hline'
            # Print sum of backgrounds
            line = "Sum Background:"
            for cat in categories_combj:
                sum_bkgd = 0
                for proc in bkgd_names:
                    sum_bkgd += combj_yields[proc,cat]
                line += ' & '+'%.2f'%sum_bkgd
            print line+' \\\\'
            del line
            print '\hline'
            # Print signal yields
            for proc in sgnl_names:
                line = proc
                for cat in categories_combj:
                    line += ' & '+'%.2f'%combj_yields[proc,cat]
                print line+' \\\\'
            print '\hline'
            # Print sum of signals
            line = "Sum Signal:"
            for cat in categories_combj:
                sum_sgnl = 0
                for proc in sgnl_names:
                    sum_sgnl += combj_yields[proc,cat]
                line += ' & '+'%.2f'%sum_sgnl
            print line+' \\\\'
            del line
            print '\hline\hline'
            # Print sum of signal and background
            line = "Sum S+B:"
            for cat in categories_combj:
                sum = 0
                for proc in sgnl_names+bkgd_names:
                    sum += combj_yields[proc,cat]
                line += ' & '+'%.2f'%sum
            print line+' \\\\'
            del line
            # Print data
            line = "Data:"
            for cat in categories_combj:
                line += ' & '+'%.0f'%combj_data[cat]
            print line+' \\\\'
            del line
            print ""

        #self.cb.PrintAll() # Print the entire datacard to terminal
        self.logger.info("Writing datacard '%s'...",self.outf)
        self.cb.WriteDatacard(self.outf)

        # Debug, plotting systematics
        if(0):
            SystogramPoints = ([],[])
            clusterPoints = ([],[],[])
            
            sys_fulltypes = []
            for ip,proc in enumerate(sgnl_names+bkgd_names):
                for cat in categories:
                    for sys in sys_dict[proc,cat]:
                        if sys != 'MCSTATS':
                            if sys not in sys_fulltypes: sys_fulltypes.append(sys)
                            if nom_dict[proc,cat] < 0: print("Found a nom zero!",nom_dict[proc,cat])
                            if sys_dict[proc,cat][sys] < 0: print("Found a sys zero!")
                            SystogramPoints[0].append(nom_dict[proc,cat])
                            SystogramPoints[1].append(sys_dict[proc,cat][sys])
                            clusterPoints[0].append(ip)
                            clusterPoints[1].append(sys_fulltypes.index(sys))
                            clusterPoints[2].append(sys_dict[proc,cat][sys])
                        
            plt.hist2d(x=SystogramPoints[0],y=SystogramPoints[1],bins=[70,80],range=[[0,35],[0,40]],cmin=1)
            plt.colorbar()
            plt.savefig('SystogramRvS.png')
            plt.close()

            plt.hist2d(x=SystogramPoints[0],y=SystogramPoints[1],bins=[100,80],range=[[0.0000001,1],[0,40]],cmin=1)
            plt.colorbar()
            plt.savefig('SystogramRvS_zoom.png')
            plt.close()

            plt.hist2d(x=np.log10(SystogramPoints[0]),y=SystogramPoints[1],bins=[100,80],range=[[-7,-0.1],[0,40]],cmin=1)
            plt.colorbar()
            plt.xlim(left=-7,right=-0.1)
            plt.savefig('SystogramRvS_zoomlog.png')
            plt.close()

            newclusterPoints = ([],[])        
            newclusterPoints[0].extend([clusterPoints[0][i] for i,kappa in enumerate(clusterPoints[2]) if kappa>2])
            newclusterPoints[1].extend([clusterPoints[1][i] for i,kappa in enumerate(clusterPoints[2]) if kappa>2])
            ax=plt.subplot(111)
            h = ax.hist2d(x=newclusterPoints[0],y=newclusterPoints[1],bins=[len(sgnl_names+bkgd_names),len(sys_fulltypes)],range=[[0,len(sgnl_names+bkgd_names)],[0,len(sys_fulltypes)]],cmin=1)
            plt.xticks([x+0.5 for x in range(len(sgnl_names+bkgd_names))],sgnl_names+bkgd_names,rotation=45)
            plt.yticks([y+0.5 for y in range(len(sys_fulltypes))],sys_fulltypes)
            plt.grid(1)
            #ax.set_xticks(range(len(sgnl_names+bkgd_names)))
            #ax.set_xticklabels(sgnl_names+bkgd_names)
            #ax.set_yticks(range(len(sys_fulltypes)))
            #ax.set_yticklabels(sys_fulltypes)
            plt.colorbar(h[3],ax=ax)
            plt.savefig('ClusterPvS.png')

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

    # Check for fake data argument
    parser = argparse.ArgumentParser()
    parser.add_argument('--fakedata', help="Flag to use fake data instead of real data", default=False, action='store_true')
    parser.add_argument('--central', help="Flag to use central signal samples instead of private", default=False, action='store_true')
    args = parser.parse_args()

    # Run datacard maker
    dm = DatacardMaker()
    #dm.make('../hist_files/anatest23_v3.root',fake_data)
    dm.make('../hist_files/anatest24_MergeLepFl.root',args.fakedata,args.central)
    #dm.make('../hist_files/TOP-19-001_unblinded_v1.root',fake_data)

    logging.info("Logger shutting down!")
    logging.shutdown()



