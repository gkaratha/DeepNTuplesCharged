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
#cc.Add("/eos/cms/store/cmst3/group/softJets/gkaratha/SoftMultiJet/DeepNtuples_v3/CRAB_UserFiles/PFC_Signal_chain_m70_dm20_13_04_25/250413_155615/0000/output_pfc_0_1.root")
cc.Add("output_pfc_0.root")

maxEvt=-1
CheckDRmatch=True
ReviewMatch=True




hall_compare = rt.TH2F("hall_compare","",3,0,3,3,0,3)
hdr_check = rt.TH2F("hdr_check","",100,0,0.05,100,0,0.05)
hidx_check = rt.TH2F("hidx_check","",50,0,50,50,0,50)
hmatch_check = rt.TH2F("hmatch_check","",5,0,5,5,0,5)
hMinDr_daughter = rt.TH1F("hMinDr_daughter","",100,0,0.5)
hPtRel_daughter = rt.TH1F("hPtRel_daughter","",100,-0.5,0.5)
hMinDr_daughter_sel = rt.TH1F("hMinDr_daughter_sel","",100,0,0.5)
hPtRel_daughter_sel = rt.TH1F("hPtRel_daughter_sel","",100,-0.5,0.5)
hMinDr_B = rt.TH1F("hMinDr_B","",100,0,2)
hPtRel_B = rt.TH1F("hPtRel_B","",100,-1,3)

hNumBoth = rt.TH1F("hNumBoth","",10,0,10)
hNumMine = rt.TH1F("hNumMine","",10,0,10)
hMinDrStd = rt.TH1F("hNumDrStd","",100,0,1)

hdrB_Both = rt.TH1F("hdrB_Both","",100,0,3)
hdrB_Mine = rt.TH1F("hdrB_Mine","",100,0,3)
hdrB_Std = rt.TH1F("hdrB_Std","",100,0,3)




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
   hall_compare.Fill(evt.isMatchedB+evt.isMatchedBB+evt.isMatchedBAndC,evt.isB+evt.isBB+evt.isGBB+evt.isLeptonicB+evt.isLeptonicB_C)
   minDR_B=100.
   for iB in range(evt.nGenBmeson):
     if minDR_B< deltaR(evt.jet_eta,evt.jet_phi,evt.genBmeson_eta[iB],evt.genBmeson_phi[iB]):
        continue
     minDR_B = deltaR(evt.jet_eta,evt.jet_phi,evt.genBmeson_eta[iB],evt.genBmeson_phi[iB])

   if CheckDRmatch: 
      nMatch=0.
      for irc in range(evt.nConstituent):
         minDR=100.
         minIdx=-1
         for ign in range(evt.nGenBmeson):
             if minDR< deltaR(evt.constituent_eta[irc],evt.constituent_phi[irc],evt.genBmeson_daughter_eta[ign],evt.genBmeson_daughter_phi[ign]):
                continue
             minDR = deltaR(evt.constituent_eta[irc],evt.constituent_phi[irc],evt.genBmeson_daughter_eta[ign],evt.genBmeson_daughter_phi[ign])
             minIdx=ign
         hdr_check.Fill(minDR,evt.constituent_BdaughterDr[irc])
         hMinDr_daughter.Fill(minDR)
         hPtRel_daughter.Fill( (evt.constituent_pt[irc]- evt.genBmeson_daughter_pt[int(minIdx)] )/ evt.genBmeson_daughter_pt[int(minIdx)] )

         if minDR<0.05: 
            nMatch+=1.0
            hidx_check.Fill(minIdx,evt.constituent_BdaughterIdx[irc])
      hmatch_check.Fill(nMatch,evt.nMatchedDaughtersB) 
      if ( nMatch>0 and evt.nMatchedDaughtersB==0 ) or (nMatch==0 and evt.nMatchedDaughtersB>0): 
          nfail+=1.0
 
   if ReviewMatch:
      FillB=True
      for irc in range(evt.nConstituent):
        if evt.constituent_BdaughterIdx[irc]>=0: 
           hMinDr_daughter_sel.Fill(evt.constituent_BdaughterDr[irc])
           hPtRel_daughter_sel.Fill( (evt.constituent_pt[irc]- evt.genBmeson_daughter_pt[int(evt.constituent_BdaughterIdx[irc])] )/ evt.genBmeson_daughter_pt[int(evt.constituent_BdaughterIdx[irc])] )
           if FillB:
              minDR_B2 = deltaR(evt.jet_eta,evt.jet_phi,evt.genBmeson_eta[ int(evt.constituent_BdaughterIdx[irc]) ],evt.genBmeson_phi[ int(evt.constituent_BdaughterIdx[irc]) ])
              relPt_B2 = (evt.jet_pt - evt.genBmeson_pt[ int(evt.constituent_BdaughterIdx[irc]) ] )/evt.genBmeson_pt[ int(evt.constituent_BdaughterIdx[irc]) ]
              hMinDr_B.Fill(minDR_B2)
              hPtRel_B.Fill(relPt_B2)
              FillB= False
  
   if evt.isB+evt.isBB+evt.isGBB+evt.isLeptonicB+evt.isLeptonicB_C>0 and evt.isMatchedB+evt.isMatchedBB+evt.isMatchedBAndC==0:
      nonly_int+=1.0
      #hMinDrStd.Fill( min(evt.daughter_dr) )
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
if CheckDRmatch: print("fail xcheck",nfail)

rt.gStyle.SetOptStat(0)
c1 = rt.TCanvas("c1","",800,600)
hall_compare.Draw("COLZ TEXT");
hall_compare.SetTitle("daughter-based vs standard match; track-based; standard")
c1.SetLogz()
c1.SaveAs("intmatch_tag.png")


c2 = rt.TCanvas("c2","",800,600)
hNumMine.Scale(1./hNumMine.Integral())
hNumBoth.Scale(1./hNumBoth.Integral())
hNumMine.Draw("HIST");
hNumBoth.SetLineColor(2)
hNumBoth.Draw("HIST sames")
hNumMine.SetTitle("# matched daughter; #matched; ")
c2.SaveAs("intmatch_nmatch.png")

c3 = rt.TCanvas("c3","",800,600)
hMinDrStd.Draw("HIST");
hMinDrStd.SetTitle("DR daughter in standard matching only; DR; ")
c3.SaveAs("intmatch_minelost.png")


c4 = rt.TCanvas("c4","",800,600)
hdrB_Both.Draw();
hdrB_Mine.SetLineColor(2)
hdrB_Mine.Draw("sames")
hdrB_Std.SetLineColor(3)
hdrB_Std.Draw("sames")
hdrB_Both.SetTitle("DR (jet,B); DR; ")
c4.SaveAs("intmatch_minDRB.png")


if CheckDRmatch:
  c2 = rt.TCanvas("c2","",800,600)
  hdr_check.Draw("COLZ");
  hdr_check.SetTitle("DR cross-check; Calculated DR; Saved DR")
  c2.SaveAs("intmatch_hdr.png")
  
  c3 = rt.TCanvas("c3","",800,600)
  hidx_check.Draw("COLZ");
  hdr_check.SetTitle("best gen Idx cross-check; Selected gen Idx; Saved gen Idx")
  c3.SaveAs("intmatch_hidx.png")
  
  c4 = rt.TCanvas("c4","",800,600)
  hmatch_check.Draw("COLZ");
  hmatch_check.SetTitle("# matched daughter cross-check; Calculated # matched daughters; Saved # matched daughters")
  c4.SetLogz()
  c4.SaveAs("intmatch_hmatch.png")

if ReviewMatch:
   c1b = rt.TCanvas("c1b","",800,600)
   hMinDr_daughter.Draw();
   hMinDr_daughter_sel.SetLineColor(2)
   hMinDr_daughter_sel.Draw("sames")
   hMinDr_daughter.SetTitle("minDR matched vs all; minDR;")
   c1b.SaveAs("intmatch_revdr.png")

   c2b = rt.TCanvas("c2b","",800,600)
   hPtRel_daughter.Draw();
   hPtRel_daughter_sel.SetLineColor(2)
   hPtRel_daughter_sel.Draw("sames")
   hPtRel_daughter.SetMinimum(1)
   hPtRel_daughter.SetTitle("rel pT matched vs all; rel pT;")
   c2b.SaveAs("intmatch_relpt.png")

   c3b = rt.TCanvas("c3b","",800,600)
   hMinDr_B.Draw();
   hMinDr_B.SetTitle("minDR matched B; minDR;")
   c3b.SaveAs("intmatch_revdrB.png")

   c4b = rt.TCanvas("c4b","",800,600)
   hPtRel_B.Draw();
   hPtRel_B.SetTitle("rel pT matched B; rel pT;")
   c4b.SaveAs("intmatch_relptB.png")

