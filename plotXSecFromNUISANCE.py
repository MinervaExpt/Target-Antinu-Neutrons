import ROOT
import array
import sys

def isMultiNeutron(mytree):

    #need 1 muon and no mesons and photons > 10 MeV in final state

    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg  = mytree.pdg

    nNeutrons = 0

    for p in range(0,nfsp):
        if(pdg[p]==2112 and Efsp[p] > 0.010): nNeutrons += 1

    return nNeutrons > 1

def getEAvail(mytree):
  nfsp = mytree.nfsp
  Efsp = mytree.E
  pdg  = mytree.pdg

  Eavail = 0
  for whichPart in range(0, nfsp):
    myPDG = pdg[whichPart]
    energy = Efsp[whichPart]
    if abs(myPDG) == 211: Eavail += energy - 139.57/1000
    if myPDG == 2212:     Eavail += energy - 938.28/1000
    if myPDG == 111:      Eavail += energy
    if myPDG == 22:       Eavail += energy

  return max(0, Eavail)

def getMuonMomentum(mytree):
    muon_mom = ROOT.TVector3()

    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg  = mytree.pdg
    px = mytree.px
    py = mytree.py
    pz = mytree.pz
    
    for p in range(0,nfsp):
        if(pdg[p]==-13):
            muon_mom.SetX(px[p])
            muon_mom.SetY(py[p])
            muon_mom.SetZ(pz[p])
            break
    return muon_mom

def isInclusiveRHC(mytree):
    nu_pdg = mytree.PDGnu
    lep_pdg = mytree.PDGLep #TODO: I don't normally check the true lepton PDG

    if(nu_pdg == -14 and lep_pdg == -13):
        return True
    
    return False

def isGoodMuon(mytree):
    #<20 degree muon
    #2 to 20 muon mom

    muonmom = getMuonMomentum(mytree)

    goodMuonMom   = muonmom.Mag()>2 and muonmom.Mag()<20
    goodMuonAngle = muonmom.Theta()*180/3.1415 < 20

    if goodMuonMom and goodMuonAngle: return True
    
    return False

mytree = ROOT.TChain("FlatTree_VARS")
for filename in sys.argv[1:]:
  mytree.Add(filename)

ptbins = [0, 0.05, 0.11, 0.17, 0.24, 0.31, 0.4, 0.5, 0.62, 0.77, 0.95, 1.18, 1.5]

mypt = ROOT.TH1D("pt","pt",len(ptbins)-1,array.array("d",ptbins))

for e in mytree:
    if not (isInclusiveRHC(e) and isGoodMuon(e) and getEAvail(e) < 0.1 and isMultiNeutron(e)): continue
    coslep = e.CosLep
    elep= e.ELep
    fScaleFactor = e.fScaleFactor
  
    P = ROOT.TMath.Sqrt(elep*elep-0.105*0.105)
    Pl = coslep*P
    Pt = ROOT.TMath.Sqrt(1-coslep*coslep)*P

    mypt.Fill(Pt,fScaleFactor)

mypt.Scale(1./mytree.GetNtrees())
mypt.GetXaxis().SetTitle("Muon p_{t} (GeV)")
mypt.GetYaxis().SetTitle("d#sigma/dp_{t} cm^2/GeV/nucleon")

outFileName = sys.argv[1]
outFileName = outFileName[outFileName.find("tune") + 5:outFileName.rfind("tune")] + ".root"
myoutput = ROOT.TFile(outFileName,"CREATE")
mypt.Write()

#Remember to bin-width normalize before plotting!
