import ROOT
import array
import sys

def isFSSignal(mytree):

    #need 1 muon and no mesons and photons > 10 MeV in final state

    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg  = mytree.pdg

    nNeutrons = 0
    nProtons = 0
    nMuons = 0
    nPhotons = 0
    nMesons = 0
    nHeavy = 0

    for p in range(0,nfsp):
        if(pdg[p]==2112 and (Efsp[p]-(939.57/1000.0)) > 0.010): 
            nNeutrons += 1
        elif(pdg[p]==2212 and (Efsp[p]-(938.272/1000.0)) > 0.120):
            nProtons += 1
        elif(abs(pdg[p])==13):
            nMuons += 1
        elif(abs(pdg[p])==211 or abs(pdg[p])==321 or abs(pdg[p])==323 or pdg[p]==111 or pdg[p]==130 or pdg[p]==310 or pdg[p]==311 or pdg[p]==313):
            nMesons += 1
        elif(pdg[p]==3112 or pdg[p]==3122 or pdg[p]==3212 or pdg[p]==3222 or pdg[p]==4112 or pdg[p]==4122 or pdg[p]==4222 or pdg[p]==411 or pdg[p]==421):
            nHeavy += 1
        elif(pdg[p]==22 and Efsp[p] > 0.010):
            nPhotons += 1


    return (nNeutrons > 0 and nMuons==1 and nMesons==0 and nProtons==0 and nPhotons==0 and nHeavy==0)

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

    goodMuonMom   = muonmom.Mag()>1.5 and muonmom.Mag()<20
    goodMuonAngle = muonmom.Theta()*180/3.1415 <= 17

    if goodMuonMom and goodMuonAngle: return True
    
    return False

mytree = ROOT.TChain("FlatTree_VARS")
for filename in sys.argv[2:]:
  mytree.Add(filename)

ptbins = [0, 0.125, 0.25, 0.38, 0.515, 0.64, 0.78, 0.94, 1.15, 1.5]

mypt = ROOT.TH1D("pt","pt",len(ptbins)-1,array.array("d",ptbins))

for e in mytree:
    #print "E"
    if not (isInclusiveRHC(e) and isGoodMuon(e) and isFSSignal(e)): continue
    #print "Got an event! YAYYYYYY!"
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

outFileName = sys.argv[1]+".root"
#outFileName = outFileName[outFileName.find("tune") + 5:outFileName.rfind("tune")] + ".root"
myoutput = ROOT.TFile(outFileName,"CREATE")
mypt.Write()

#Remember to bin-width normalize before plotting!
