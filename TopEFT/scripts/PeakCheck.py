import ROOT
ROOT.gSystem.Load('$CMSSW_BASE/src/CombineHarvester/TopEFT/interface/TH1EFT_h.so')

infile = '$CMSSW_BASE/src/CombineHarvester/TopEFT/hist_files/anatest23_v3_MergeLepFl.root'
readfile = ROOT.TFile.Open(infile)

sgnl_known = ['ttH','tllq','ttll','ttlnu','tHq']
sgnl_histnames = [sgnl + '_' + '16D' for sgnl in sgnl_known] # For private samples
#bkgd_known = ['charge_flips','fakes','Diboson','Triboson','convs']
bkgd_known = ['charge_flips','fakes','WZ','WWW','ttGJets']

def name_bin(category,bin):
    # For standard histogram files
    if "2lss" in category:
        return 'C_{0}_{1}{2}j'.format(category, 'ge' if bin==4 else '', bin+3)
    if "3l" in category:
        return 'C_{0}_{1}{2}j'.format(category, 'ge' if bin==4 else '', bin+1)
    if "4l" in category:
        return 'C_{0}_{1}{2}j'.format(category, 'ge' if bin==4 else '', bin)

# Lists of names to pass on
categories=[] #e.g. C_2lss_ee_2j_1b
sys_types=[] #e.g. Lumi

# Dicts for rates
peak_dict = {}
nonpeak_dict = {}

# Main parsing loop
for key in readfile.GetListOfKeys():
    hist = readfile.Get(key.GetName())

    # Get categorical information from histogram name
    histname = hist.GetName().split('.')
    category,systematic,process = '','',''
    if(len(histname)==3): [category,systematic,process] = histname
    if(len(histname)==2): [category,process] = histname
    
    # Rename processes to be fed into DatacardMaker
    # For renaming central sample processes:
    #process = process.replace('tZq','tllq')
    #process = process.replace('ttZ','ttll')
    #process = process.replace('ttW','ttlnu')
    # For accurate naming of backgrounds:
    #process = process.replace('ttGJets','convs')
    #process = process.replace('WZ','Diboson')
    #process = process.replace('WWW','Triboson')
    # For accurate naming of systematics:
    systematic = systematic.replace('FR','FR_FF') # Don't touch this without changing below!

    # For standard histogram files
    maxbin = 4
    # For anatest19
    #if '2lss' in category: maxbin=2
    #if '3l' in category: maxbin=2
    #if '4l' in category: maxbin=4

    if process in ['fakes','charge_flips','data']: continue
    if systematic == '':
        # Find peak bin
        bin_yields = []
        for bin in range(1,maxbin+1):
            category_njet = name_bin(category,bin)
            if category_njet not in categories: categories.append(category_njet)
            bin_yield = round(hist.GetBinContent(bin,ROOT.WCPoint()),8)
            jesdown_yield = round(readfile.Get(category+'.JESDOWN.'+process).GetBinContent(bin,ROOT.WCPoint()),8)
            jesup_yield = round(readfile.Get(category+'.JESUP.'+process).GetBinContent(bin,ROOT.WCPoint()),8)
            bin_yields.append(bin_yield)
            nonpeak_dict[(process,category_njet)] = (bin_yield,jesdown_yield,jesup_yield)
        peakbin = bin_yields.index(max(bin_yields))
        peakcat = name_bin(category,peakbin+1)
        if not (peakbin==0 and process in bkgd_known): # Ignore bottom bin peaks for backgrounds
            peak_dict[(process,peakcat)] = nonpeak_dict.pop((process,peakcat))
        else:
            print "Skipping process {} bin {} for bottom bin peak".format(process,category)

ss=0.
ss_2l=0.
ss_3l=0.
ss_4l=0.
os=0.
os_2l=0.
os_3l=0.
os_4l=0.
for key in peak_dict:
    if peak_dict[key][0] <=0: continue
    ratiodown = peak_dict[key][1]/peak_dict[key][0]
    ratioup = peak_dict[key][2]/peak_dict[key][0]
    if (ratiodown>1 and ratioup>1) or (ratiodown<1 and ratioup<1):
        ss = ss+1
        if '2l' in key[1]: ss_2l = ss_2l+1
        if '3l' in key[1]: ss_3l = ss_3l+1
        if '4l' in key[1]: ss_4l = ss_4l+1
    else:
        os = os+1
        if '2l' in key[1]: os_2l = os_2l+1
        if '3l' in key[1]: os_3l = os_3l+1
        if '4l' in key[1]: os_4l = os_4l+1
print "NPeaks same direction: ",ss,"  Occurance = {}%".format(100*ss/(ss+os))
print "NPeaks 2l same direction: ",ss_2l,"  Occurance = {}%".format(100*ss_2l/(ss_2l+os_2l))
print "NPeaks 3l same direction: ",ss_3l,"  Occurance = {}%".format(100*ss_3l/(ss_3l+os_3l))
print "NPeaks 4l same direction: ",ss_4l,"  Occurance = {}%".format(100*ss_4l/(ss_4l+os_4l))
print "NPeaks opposite direction:",os,"  Occurance = {}%".format(100*os/(ss+os))
print "NPeaks 2l opposite direction: ",os_2l,"  Occurance = {}%".format(100*os_2l/(ss_2l+os_2l))
print "NPeaks 3l opposite direction: ",os_3l,"  Occurance = {}%".format(100*os_3l/(ss_3l+os_3l))
print "NPeaks 4l opposite direction: ",os_4l,"  Occurance = {}%".format(100*os_4l/(ss_4l+os_4l))
ss=0.
ss_2l=0.
ss_3l=0.
ss_4l=0.
os=0.
os_2l=0.
os_3l=0.
os_4l=0.
for key in nonpeak_dict:
    if nonpeak_dict[key][0] <=0: continue
    ratiodown = nonpeak_dict[key][1]/nonpeak_dict[key][0]
    ratioup = nonpeak_dict[key][2]/nonpeak_dict[key][0]
    if (ratiodown>1 and ratioup>1) or (ratiodown<1 and ratioup<1):
        ss = ss+1
        if '2l' in key[1]: ss_2l = ss_2l+1
        if '3l' in key[1]: ss_3l = ss_3l+1
        if '4l' in key[1]: ss_4l = ss_4l+1
    else:
        os = os+1
        if '2l' in key[1]: os_2l = os_2l+1
        if '3l' in key[1]: os_3l = os_3l+1
        if '4l' in key[1]: os_4l = os_4l+1
print "NOffpeaks same direction: ",ss,"  Occurance = {}%".format(100*ss/(ss+os))
print "NOffpeaks 2l same direction: ",ss_2l,"  Occurance = {}%".format(100*ss_2l/(ss_2l+os_2l))
print "NOffpeaks 3l same direction: ",ss_3l,"  Occurance = {}%".format(100*ss_3l/(ss_3l+os_3l))
print "NOffpeaks 4l same direction: ",ss_4l,"  Occurance = {}%".format(100*ss_4l/(ss_4l+os_4l))
print "NOffpeaks opposite direction:",os,"  Occurance = {}%".format(100*os/(ss+os))
print "NOffpeaks 2l opposite direction: ",os_2l,"  Occurance = {}%".format(100*os_2l/(ss_2l+os_2l))
print "NOffpeaks 3l opposite direction: ",os_3l,"  Occurance = {}%".format(100*os_3l/(ss_3l+os_3l))
print "NOffpeaks 4l opposite direction: ",os_4l,"  Occurance = {}%".format(100*os_4l/(ss_4l+os_4l))





