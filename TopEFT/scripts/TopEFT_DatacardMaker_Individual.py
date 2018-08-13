import sys
import CombineHarvester.CombineTools.ch as ch
import CombineHarvester.TopEFT.Process_Input_Individual as pi

#Check for coefficient argument
if len(sys.argv) == 1:
    raise Exception("Must have WC as an argument.")
if len(sys.argv) > 2:
    raise Exception("Only one argument (the WC) is accepted.")
operators_known = ['ctZ','ctW','ctp','ctl1','ctG','cQe1','cpt','cptb','cpQM','cpQ3']
if sys.argv[1] not in operators_known:
    raise Exception("Specified WC not in known list of WC's.")

debug = 0
    

#Parse input file
print "Parsing input file for operator",sys.argv[1],"..."
(categories, data_names, data_dict, sgnl_names, bkgd_names, nom_dict, sys_types, sys_dict) = pi.process_root_input('../data/anatest7.root',sys.argv[1])
print "Done parsing input file."
print "Now creating Datacard..."

#Initialize CombineHarvester instance
cb = ch.CombineHarvester()
cb.SetVerbosity(0)

#Declare lists of analysis divisions
eras = ['2016'] #Data eras

#Manually set signal and background processes... chance for automation here
#sig_procs = proc_names[-4:] #Signal Processes
#bkgd_procs = proc_names[:-4] #Background Processes

chan = [''] # Indistiguishable process subcategories i.e. for MC only

cats = list(enumerate(categories)) #Process bins. Must be list of tuple like this

#Asimov data
obs_rates={}
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
cb.AddObservations( ['*'], ['Top_EFT'], eras, chan, cats )

cb.AddProcesses( ['*'], ['Top_EFT'], eras, chan, sgnl_names, cats, True) 
cb.AddProcesses( ['*'], ['Top_EFT'], eras, chan, bkgd_names, cats, False)

#Fill the data rates
print "Adding observation rates..."
cb.ForEachObs(lambda x: x.set_rate(obs_rates[x.bin()])) #Not split by process/channel, obviously

#Fill the nominal MC rates
print "Adding MC rates..."
cb.ForEachProc(lambda x: x.set_rate(nom_dict[x.process(),x.bin()]))

#Fill systematic rates
for proc in sgnl_names+bkgd_names:
    print "Adding systematics for",proc,"..."
    for cat in categories:
        #Lumi uncertainty (fully correlated, identical for all categories)
        cb.cp().process([proc]).bin([cat]).AddSyst(cb,'Lumi','lnN',ch.SystMap()( 1.025 ))
        #JES uncertainty (fully correlated, asymmetric)
        cb.cp().process([proc]).bin([cat]).AddSyst(cb,'JES','lnN',ch.SystMap()( [sys_dict[(proc,cat)]['JESUP'], sys_dict[(proc,cat)]['JESDOWN']]))
        if (proc,cat) in sys_dict.keys(): # probably unnecessary safeguard
            #Other systematics -- To be added
            for sys_type in sys_types:
                #Fully correlated systematics listed by rate
                if sys_type in []:
                    cb.cp().process([proc]).bin([cat]).AddSyst(cb,sys_type,'lnN',ch.SystMap()( sys_dict[(proc,cat)][sys_type] ))
                #Systematics listed by their fluctuations
                elif sys_type in ['pdfUP']:
                    unc = (nom_dict[proc,cat]+sys_dict[proc,cat][sys_type])/nom_dict[proc,cat] if sys_dict[proc,cat][sys_type]>0. else 1.0000
                    cb.cp().process([proc]).bin([cat]).AddSyst(cb,sys_type.rstrip('UP'),'lnN',ch.SystMap()(unc))
                #Shape systematics (not used in counting experiment)
                #elif sys_type in ['MCStatUP','MCStatDOWN','Q2UP','Q2DOWN']:
                #Superfluous systematics
                #Fully uncorrelated systematics
                #else:
                    #cb.cp().process([proc]).bin([cat]).AddSyst(cb,proc+cat+':'+sys_type,'lnN',ch.SystMap()( float(sys_dict[(proc,cat)][sys_type]) ))

#Printout of signal and background yields (debug)
if debug:
    background = {}
    signal = {}
    print "\nCategory, signal yield, background yield:"
    for cat in categories:
        background[cat] = 0
        signal[cat] = 0
        for proc in bkgd_names:
            background[cat] += nom_dict[proc,cat]
        for proc in sgnl_names:
            signal[cat] += nom_dict[proc,cat]
        print cat,signal[cat],background[cat]

#cb.PrintAll() #Print the datacard
print "Writing datacard '{}'...".format("Datacard_root_"+sys.argv[1]+".txt")
cb.WriteDatacard('Datacard_root_'+sys.argv[1]+'.txt')

