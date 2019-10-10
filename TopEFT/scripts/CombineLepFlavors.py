import ROOT
ROOT.gSystem.Load('$CMSSW_BASE/src/CombineHarvester/TopEFT/interface/TH1EFT_h.so')

class CombineLepFlavors(object):
    def __init__(self):
        self.Filename = "anatest24.root"

    def execute(self):
        File = ROOT.TFile.Open('$CMSSW_BASE/src/CombineHarvester/TopEFT/hist_files/{}'.format(self.Filename))
        outFile = ROOT.TFile.Open('$CMSSW_BASE/src/CombineHarvester/TopEFT/hist_files/{}_MergeLepFl.root'.format(self.Filename.strip('.root')),'RECREATE')
        outFile.cd()

        for key in File.GetListOfKeys():
            hist = File.Get(key.GetName())
            histname = key.GetName().split('.')
            category,systematic,process = '','',''
            if(len(histname)==3): [category,systematic,process] = histname
            if(len(histname)==2): [category,process] = histname
            if '2l' in category:
                if 'ee' in category:
                    emu = File.Get(key.GetName().replace('ee','emu'))
                    mumu = File.Get(key.GetName().replace('ee','emu'))
                    hist.Add(emu)
                    hist.Add(mumu)
                    hist.SetName(key.GetName().replace('_ee',''))
                    outFile.cd()
                    hist.Write()
            else:
                outFile.cd()
                hist.Write()
                
        outFile.Close()

if __name__ == "__main__":

    CLF = CombineLepFlavors()
    CLF.execute()
