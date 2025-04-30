import ROOT as rt
import math

def deltaR(gen_eta, gen_phi, reco_eta, reco_phi):
    delta_eta = gen_eta - reco_eta
    delta_phi = gen_phi - reco_phi
    # Adjust delta_phi to account for periodicity (wrap into [-pi, pi])
    delta_phi = (delta_phi + math.pi) % (2 * math.pi) - math.pi
    dr = math.sqrt(delta_eta**2 + delta_phi**2)
    return dr


rt.gROOT.SetBatch(True)

cc = rt.TChain("deepntuplizer/tree")
cc.Add("output_pfc_0.root")

maxEvt=-1



hNdaughter = rt.TH1F("hNdaughter","",10,0,10)
hallb_compare = rt.TH2F("hallb_compare","",3,0,3,3,0,3)
hallc_compare = rt.TH2F("hallc_compare","",3,0,3,3,0,3)

hNumBoth = rt.TH1F("hNumBoth","",10,0,10)
hNumMine = rt.TH1F("hNumMine","",10,0,10)
hMinDrStd = rt.TH1F("hNumDrStd","",100,0,1)

hdrB_Both = rt.TH1F("hdrB_Both","",100,0,3)
hdrB_Mine = rt.TH1F("hdrB_Mine","",100,0,3)
hdrB_Std = rt.TH1F("hdrB_Std","",100,0,3)

hDenMatched = rt.TH1F("hDenMatched","",10,0,10)
hNumMatched = rt.TH1F("hNumMatched","",10,0,10)

h2DfractionNum = rt.TH2F("h2DfractionNum","",9,1,10,100,0,1)
h2DdrNum = rt.TH2F("h2DdrNum","",10,0,10,100,0,3)
h2DptNum = rt.TH2F("h2DptNum","",9,1,10,40,0,40)


nboth=0.0
nonly_mine=0.0
nonly_int=0.0
nsuccess=0.
nfail=0.
ievt=0.
for evt in cc:
   if ievt%100==0: print(ievt,"/",cc.GetEntries())
   ievt+=1
   if ievt==maxEvt: break
   lastIdx=-1
   nDaughterB=-1
 #  basedNdaughters={1:[],2:[],3:[],4:[],5:[],6:[]}
 #  idx_list=[]
   for iB in range(evt.nGenBmeson):
     if lastIdx!=evt.genBmeson_daughter_mesonIdx[iB]:
        hNdaughter.Fill(nDaughterB)
#        basedNdaughters[nDaughterB] = basedNdaughters[nDaughterB] + idx_list
        lastIdx=evt.genBmeson_daughter_mesonIdx[iB]
        nDaughterB=1
  #      idx_list=[iB]
     else:
        nDaughterB+=1
   #     idx_list.append(iB)
   
   
   if evt.nGenBmeson==0: print("0 evt")
   hallb_compare.Fill(evt.isMatchedB+evt.isMatchedBB+evt.isMatchedBAndC,evt.isB+evt.isBB+evt.isGBB+evt.isLeptonicB+evt.isLeptonicB_C)
   hallc_compare.Fill(evt.isMatchedC+evt.isMatchedBAndC,evt.isC+evt.isCC+evt.isGCC)
   
   minDR_B=100.
   for iB in range(evt.nGenBmeson):
     if minDR_B< deltaR(evt.jet_eta,evt.jet_phi,evt.genBmeson_eta[iB],evt.genBmeson_phi[iB]):
        continue
     minDR_B = deltaR(evt.jet_eta,evt.jet_phi,evt.genBmeson_eta[iB],evt.genBmeson_phi[iB])

   if evt.isMatchedB+evt.isMatchedBB+evt.isMatchedBAndC+evt.isMatchedC>1:
      print("issue matchedB:",evt.isMatchedB,"BB:",evt.isMatchedBB,"B+C:",evt.isMatchedBAndC,"C:",evt.isMatchedC)
   hDenMatched.Fill(evt.nMatchedDaughtersB)
   if evt.isB+evt.isBB+evt.isGBB+evt.isLeptonicB+evt.isLeptonicB_C>0:
      hNumMatched.Fill(evt.nMatchedDaughtersB)
   h2DfractionNum.Fill(evt.nMatchedDaughtersB,evt.nMatchedDaughtersB*1.0/evt.nConstituent)
   h2DdrNum.Fill(evt.nMatchedDaughtersB,minDR_B)
   h2DptNum.Fill(evt.nMatchedDaughtersB,evt.jet_pt)
   if evt.isB+evt.isBB+evt.isGBB+evt.isLeptonicB+evt.isLeptonicB_C>0 and evt.isMatchedB+evt.isMatchedBB+evt.isMatchedBAndC==0:
      nonly_int+=1.0
      hdrB_Std.Fill(minDR_B)
   if evt.isB+evt.isBB+evt.isGBB+evt.isLeptonicB+evt.isLeptonicB_C==0 and evt.isMatchedB+evt.isMatchedBB+evt.isMatchedBAndC>0:
      nonly_mine+=1.0
      hNumMine.Fill(evt.nMatchedDaughtersB)
      hdrB_Mine.Fill(minDR_B)
   if evt.isB+evt.isBB+evt.isGBB+evt.isLeptonicB+evt.isLeptonicB_C>0 and evt.isMatchedB+evt.isMatchedBB+evt.isMatchedBAndC>0:
      nboth+=1.0
      hNumBoth.Fill(evt.nMatchedDaughtersB)
      hdrB_Both.Fill(minDR_B)

print("both",nboth,"only mine",nonly_mine,"only int",nonly_int)

rt.gStyle.SetOptStat(0)

c0 = rt.TCanvas("c0","",800,600)
hNdaughter.Draw()
hNdaughter.SetTitle("Gen daughters;N;")
c0.SaveAs("comlabels_Nd.png")

c1 = rt.TCanvas("c1","",800,600)
hallb_compare.Scale(1./hallb_compare.Integral())
hallb_compare.Draw("COLZ TEXT");
hallb_compare.SetTitle("daughter-based vs standard match; track-based; standard")
c1.SetLogz()
c1.SaveAs("comlabels_btag.png")

c2 = rt.TCanvas("c2","",800,600)
hallc_compare.Scale(1./hallc_compare.Integral())
hallc_compare.Draw("COLZ TEXT");
hallc_compare.SetTitle("daughter-based vs standard match; track-based; standard")
c2.SetLogz()
c2.SaveAs("comlabels_ctag.png")


c3 = rt.TCanvas("c3","",800,600)
hNumMine.Scale(1./hNumMine.Integral())
hNumBoth.Scale(1./hNumBoth.Integral())
hNumMine.Draw("HIST");
hNumBoth.SetLineColor(2)
hNumBoth.Draw("HIST sames")
hNumMine.SetTitle("# matched daughter; #matched; ")
c3.SaveAs("comlabels_nmatch.png")


c4b = rt.TCanvas("c4b","",800,600)
hdrB_Both.Draw();
hdrB_Mine.SetLineColor(2)
hdrB_Mine.Draw("sames")
hdrB_Std.SetLineColor(3)
hdrB_Std.Draw("sames")
hdrB_Both.SetTitle("DR (jet,B); DR; ")
c4b.SaveAs("comlabels_minDRB.png")

hNumMatched.Divide(hDenMatched)
c4 = rt.TCanvas("c4","",800,600)
hNumMatched.Draw("HIST");
hNumMatched.SetTitle("Label efficiency vs #matched daughters;#daughters ; Efficiency")
c4.SaveAs("comlabels_matchedeff.png")

c5 = rt.TCanvas("c5","",800,600)
h2DfractionNum.Draw("COLZ")
h2DfractionNum.SetTitle(";#matched daughters ; #matched / all costituents")
c5.SaveAs("comlabels_fraction.png")


c6 = rt.TCanvas("c5b","",800,600)
h2DdrNum.Draw("COLZ")
h2DdrNum.SetTitle(";#matched daughters ; DR(B,jet)")
c6.SaveAs("comlabels_2dDR.png")

c7 = rt.TCanvas("c7","",800,600)
h2DptNum.Draw("COLZ")
h2DptNum.SetTitle(";#matched daughters ; Pt(reco)")
c7.SaveAs("comlabels_pt.png")
