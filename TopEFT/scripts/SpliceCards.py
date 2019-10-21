import ROOT
ROOT.gSystem.Load('$CMSSW_BASE/src/CombineHarvester/TopEFT/interface/TH1EFT_h.so')

class CardSplicer(object):
    def __init__(self):
        self.Filename2 = "anatest16_MergeLepFl.root"
        self.Filename1 = "TOP-19-001_unblinded_v1_MergeLepFl.root"

    def splice(self):
        File1 = ROOT.TFile.Open('$CMSSW_BASE/src/CombineHarvester/TopEFT/hist_files/{}'.format(self.Filename1))
        File2 = ROOT.TFile.Open('$CMSSW_BASE/src/CombineHarvester/TopEFT/hist_files/{}'.format(self.Filename2))
        outFile = ROOT.TFile.Open('$CMSSW_BASE/src/CombineHarvester/TopEFT/hist_files/spliced.root','RECREATE')
        outFile.cd()

        for key in File1.GetListOfKeys():
            hist = File1.Get(key.GetName())
            histname = hist.GetName()
            if 'tHq' in histname:
                outFile.cd()
                hist.Write()

        for key in File2.GetListOfKeys():
            hist = File2.Get(key.GetName())
            histname = hist.GetName()
            if not 'tHq' in histname:
                outFile.cd()
                hist.Write()

if __name__ == "__main__":

    CS = CardSplicer()
    CS.splice()
