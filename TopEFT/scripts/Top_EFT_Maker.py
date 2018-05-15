import CombineHarvester.CombineTools.ch as ch
from Process_Input import process_input

#Parse input file
print "Parsing input file..."
(categories, proc_names, proc_dict, sys_types, sys_dict) = process_input('../data/test_input_for_datacard_2.txt')

#Initialize CombineHarvester instance
cb = ch.CombineHarvester()
cb.SetVerbosity(0)

#Declare lists of analysis divisions
eras = ['2016'] #Data eras

#sig_procs = proc_names[-7:] #Signal Processes
#bkgd_procs = proc_names[:-7] #Background Processes
#sig_procs = proc_names[-1:] #Signal Processes
#bkgd_procs = proc_names[:-1] #Background Processes
sig_procs = proc_names[:4] #Signal Processes
bkgd_procs = proc_names[4:] #Background Processes

chan = [''] # Indistiguishable process subcategories i.e. for MC only

cats = list(enumerate(categories)) #Process bins. Must be list of tuple

#Placeholder for cross-sections/rates until we have an input
obs_rates={}
for cat in categories:
    cat_asimov = 0
    for proc in proc_names:
        cat_asimov += proc_dict[proc,cat]
    obs_rates[cat]=cat_asimov

#Initialize structure of each observation (data) and process (MC signal and background)
#Can be adjusted very heavily to specialize each era/process/category
cb.AddObservations( ['*'], ['Top_EFT'], eras, chan, cats )

cb.AddProcesses( ['*'], ['Top_EFT'], eras, chan, sig_procs, cats, True) 
cb.AddProcesses( ['*'], ['Top_EFT'], eras, chan, bkgd_procs, cats, False)

#Fill the nominal rates
print "Adding observation rates..."
cb.ForEachObs(lambda x: x.set_rate(obs_rates[x.bin()])) #Not split by process/channel, obviously
print "Adding MC rates..."
cb.ForEachProc(lambda x: x.set_rate(proc_dict[x.process(),x.bin()]))

#Fill systematic rates
for proc in proc_names:
    print "Adding systematics for",proc,"..."
    for cat in categories:
        #print "    in category",cat,"..."
        if (proc,cat) in sys_dict.keys(): # probably unnecessary safeguard
            for sys_type in sys_types:
                #Fully correlated systematics listed by rate FIXME
                if sys_type in []:
                    cb.cp().process([proc]).bin([cat]).AddSyst(cb,sys_type,'lnN',ch.SystMap()( sys_dict[(proc,cat)][sys_type] ))
                #Systematics listed by their fluctuations
                elif sys_type in ['LumiUP','pdfUP']:
                    unc = (proc_dict[proc,cat]+sys_dict[proc,cat][sys_type])/proc_dict[proc,cat] if sys_dict[proc,cat][sys_type]>0. else 1.0000
                    cb.cp().process([proc]).bin([cat]).AddSyst(cb,sys_type.rstrip('UP'),'lnN',ch.SystMap()(unc))
                #Shape systematics (not used in counting experiment)
                #elif sys_type in ['MCStatUP','MCStatDOWN','Q2UP','Q2DOWN']:
                #Superfluous systematics
                #elif sys_type in ['LumiDOWN','pdfDOWN']:
                #Fully uncorrelated systematics
                #else:
                    #cb.cp().process([proc]).bin([cat]).AddSyst(cb,proc+cat+':'+sys_type,'lnN',ch.SystMap()( float(sys_dict[(proc,cat)][sys_type]) ))


#cb.PrintAll()
print "Writing datacard '{}'...".format("Datacard_test.txt")
#cb.WriteDatacard('EFT_datacard.txt')
cb.WriteDatacard('Datacard_test.txt')

