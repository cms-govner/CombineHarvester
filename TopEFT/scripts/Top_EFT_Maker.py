import CombineHarvester.CombineTools.ch as ch
from Process_Input import process_input

#Parse input file
print "Parsing input file..."
(categories, proc_names, proc_values, sys_dict, sys_types) = process_input('../data/test_input_for_datacard_2.txt')

#Initialize CombineHarvester instance
cb = ch.CombineHarvester()
cb.SetVerbosity(0)

#Declare lists of analysis divisions
eras = ['2016'] #Data eras

bkgd_procs = proc_names[:-7]
sig_procs = proc_names[-7:] #Signal Processes

chan = [''] # Indistiguishable process subcategories i.e. for MC only

cats = list(enumerate(categories)) #Process bins. Must be list of tuple

#Placeholder for cross-sections/rates until we have an input
obs_rates = dict(zip(categories,range(0,len(categories)))) #Placeholder
mc_rate = dict(zip(categories,proc_values)) #Not a placeholder

#Initialize structure of each observation (data) and process (MC signal and background)
#Can be adjusted very heavily to specialize each era/process/category
cb.AddObservations( ['*'], ['Top_EFT'], eras, chan, cats )

cb.AddProcesses( ['*'], ['Top_EFT'], eras, chan, sig_procs, cats, True) 
cb.AddProcesses( ['*'], ['Top_EFT'], eras, chan, bkgd_procs, cats, False)

#Fill the rates
print "Adding observation rates..."
cb.ForEachObs(lambda x: x.set_rate(obs_rates[x.bin()])) #Not split by process/channel, obviously
#Modify to be a for loop over list of processes; add list items to input of lookup
print "Adding MC rates..."
cb.ForEachProc(lambda x: x.set_rate(mc_rate[x.bin()]))

#Add systematics
#These are split by category
for proc in proc_names:
    print "Adding systematics for",proc,":"
    for cat in categories:
        print "    in category",cat,"..."
        if (proc,cat) in sys_dict.keys(): # probably unnecessary
            for sys_type in sys_dict[(proc,cat)].keys():
                #Fully correlated systematics
                if sys_type in ['LumiUP','LumiDOWN','pdfUP','pdfDOWN','Q2UP','Q2DOWN','MCStatUP','MCStatDOWN']:
                    cb.cp().process([proc]).bin([cat]).AddSyst(cb,sys_type,'lnN',ch.SystMap()( float(sys_dict[(proc,cat)][sys_type]) ))
                #Fully uncorrelated systematics (MCStatUP/DOWN)
                else:
                    cb.cp().process([proc]).bin([cat]).AddSyst(cb,proc+cat+':'+sys_type,'lnN',ch.SystMap()( float(sys_dict[(proc,cat)][sys_type]) ))


#cb.PrintAll()
cb.WriteDatacard('EFT_datacard.txt')

